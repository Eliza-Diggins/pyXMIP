.. _Astrometric_reduction:
===============================
Astrometric Reduction Process
===============================

The astrometric reduction process is designed to take into account the instrumental parameters of the observing instrument to
deal with error is astrometry which might lead to inaccurate matches. In effect, this reduction process boils down to determining how
confident one can be in a proposed match given their separation on the sky.

.. raw:: html

   <hr>

Overview
--------

Statistical inference regarding our confidence in the position of sources is a non-trivial problem. In ``pyXMIP``, we choose to treat
it in a `Bayesian framework <https://en.wikipedia.org/wiki/Bayesian_inference>`_. The basic premise of the method is to compare two hypotheses
for each proposed match:

1. (:math:`H_0`, the null-hypothesis) These two sources (the catalog source and the match candidate) are not the same source.
2. (:math:`H_1`, the test-hypothesis) These two sources are the same source.

For each hypothesis, the probability that the two sources would have the observed separation is calculated, :math:`P(\delta|H_i)`,
where :math:`\delta` is the separation between the two sources. Using these, we are able to calculate the `Bayes factor <https://en.wikipedia.org/wiki/Bayes_factor>`_
and from that we can obtain the posterior distribution for :math:`H_1`.

Because the model used is Bayesian, a `prior <https://en.wikipedia.org/wiki/Prior_probability>`_ is required and may be user specified or estimated
from the available data.

Mathematics Basis
'''''''''''''''''

A more complete picture of the mathematical ideas used in this process is presented below; however, all users should have a basic grasp of the core
ideas of this methodology.

.. dropdown:: Show the Mathematics!

    To begin, let's pose a basic question:

    Let a catalog object (from our survey) be detected at :math:`\textbf{x}` on the sky and let :math:`\textbf{y}` be the position of a proposed match
    from a database (as observed by some other instrument).

    .. hint::

        It's important to remember that neither position is the TRUE position, they are both sampled based on the parameters of the observing instrument.

    Let's assume that we know the **astrometric distributions** of the instruments, the probability distribution :math:`P(\textbf{r}|\textbf{R},M)`, where
    :math:`\textbf{r}` is the **observed** position, :math:`\textbf{R}` is the **true** position and :math:`M` is the set of parameters determined by the
    nature of the observatory. The natural question to ask is the following:

        Given two observations as described above and knowledge of the relevant astrometry, what is the probability that they
        are actually from the same source given their separation on the sky?


    This is actually a quite difficult question to answer; in this project, we follow Budavari et al. [BuSz08]_ in our approach. Consider the following precise statement of the
    problem

        Let :math:`\textbf{x}` and :math:`\textbf{y}` be two observed positions on the sky; one corresponding to the catalog object, the other to the
        proposed match to the catalog object. Let :math:`P(\textbf{x}|\textbf{R},M_{cat})` and :math:`P(\textbf{y}|\textbf{R},M_{db})` be the astrometric distributions
        given a known true position :math:`\textbf{R}` and instrumental parameters determined by models :math:`M_{cat}` and :math:`M_{db}`. What is the probability
        that the sources identified at :math:`\textbf{x}` and :math:`\textbf{y}` are the same object?

    Two hypotheses are considered:

    1. :math:`H`, These two sources are actually associated with the same underlying physical object.
    2. :math:`H_0`, These two sources are not associated with the same underlying physical object.

    To compare the hypotheses, we utilize the Bayes factor:

    .. math::

        B(H,H_0|D) = \frac{P(H|D)P(H_0)}{P(H_0|D)P(H)} = \frac{P(D|H)}{P(D|H_0)}.

    The probability of obtaining the observed data given the hypothesis is correct is simply

    .. math::

        P(D|H) = \int_{S^2} P(\mathbf{R}) P(\mathbf{x}|\mathbf{R},M_{cat})P(\mathbf{y}|\mathbf{R},M_{db}) d\mathbf{R}.

    Here :math:`P(R)` is the probability of the underlying physical source being positioned at a given position.

    Similarly, the probability assuming :math:`H_0` is

    .. math::

        P(D|H_0) = \int_{S^2}\int_{S^2} P(\mathbf{U})P(\mathbf{V}) P(\mathbf{x}|\mathbf{U},M_{cat})P(\mathbf{y}|\mathbf{V},M_{db}) d\textbf{U}d\textbf{V}.

    .. admonition:: Key Assumptions

        Assume (as is a reasonable thing to do) that :math:`P(\textbf{x}|\textbf{R},M) \sim f(|\textbf{x}-\textbf{r}|)` and is therefore symmetric about
        the displacement between the observation and the true position. In this case, the above integral simplifies becase

        .. math::

            P(D|H_0) = \int_{S^2}\int_{S^2} P(\mathbf{U})P(\mathbf{V}) P(\mathbf{x}|\mathbf{U},M_{cat})P(\mathbf{y}|\mathbf{V},M_{db}) d\textbf{U}d\textbf{V} = P(D|H_0) = \int_{S^2}\int_{S^2} P(\mathbf{U})P(\mathbf{V}) P(\mathbf{x}|\mathbf{U},M_{cat})P(\mathbf{y}|\mathbf{V},M_{db}) d\textbf{x}d\textbf{y}.
            = P(\mathbf{U})P(\mathbf{V}).

    Therefore, one arrives at the following Bayes factor:

    .. math::

        \boxed{B(H,H_0|D) = \frac{1}{P(\mathbf{R}|H)P(\mathbf{R}|H_0)} \int_{S^2} P(\mathbf{R}|H) P(\mathbf{x}|\mathbf{R},M_{cat})P(\mathbf{y}|\mathbf{R},M_{db}) d\mathbf{R}. }

    It is possible to compute this Bayes factor and convert it into a posterior distribution (see :ref:`mathematics` below).

Usage
-----

Like all reduction processes, the astrometric reduction is a :py:class:`structures.reduction.ReductionProcess` subclass (:py:class:`structures.reduction.AstrometricReductionProcess`).

.. hint::

    It's worth taking a look at the API documentation for :py:class:`structures.reduction.AstrometricReductionProcess`,
    there are a great many parameters which might be very helpful for your particular use case.

Setup and Parameters
''''''''''''''''''''

To get started, it's important to describe all of the relevant parameters and their behavior. The first consideration is how
to model the uncertainty in the positions of each source.

.. important::

    There are a lot of parameters here, but you usually won't need to specify them all by hand. If you run this process through
    a :py:class:`cross_reference.CrossMatchDatabase`, then the corresponding :py:class:`schema.CMDSchema` will be able to provide
    a lot of the necessary information.

Modeling Astrometric Uncertainty
++++++++++++++++++++++++++++++++

Every astrometric reduction process has two key parameters:

1. :py:attr:`~structures.reduction.AstrometricReductionProcess.CATALOG_ERR`: The position uncertainty in the catalog sources.
2. :py:attr:`~structures.reduction.AstrometricReductionProcess.DATABASE_ERR`: The position uncertainty in the database source.

Both of these parameters are :py:class:`utilities.types.CoordinateErrorSpecifier` instances.


.. raw:: html

   <hr>

Usage
-----

In this section, we present the details concerning how to set up one of these reduction processes.

Setup and Parameters
''''''''''''''''''''

The astrometric reduction process is highly customizable, meaning you can make it do whatever you need; however, that comes at the cost
of having **lots** of potential options! In this section, we'll go through the key considerations for how to setup the reduction process.

Astrometry Modes
++++++++++++++++

In some astronomical databases, the positions of objects may be specified by a single value (something like ``POSERR``). In others, you might
have errors provided in both RA and DEC. For some unfortunate cases, you might only have a general understanding of the instrumental PSF that has to
be provided to the reduction process by hand.

Each of these options corresponds to a given **astrometric-mode**.

- The ``CATALOG`` in your :py:class:`cross_reference.CrossMatchDatabase` has an astrometric mode.

  - Determined by the :py:attr:`~reduction.AstrometricReduction.astrometry_mode_cat` attribute of the reduction class.

- For **each** reference database, there is a different astrometry mode.

  - Determined by the :py:attr:`~reduction.AstrometricReduction.astrometry_mode_db` attribute of the reduction class.

There are 3 options for the astrometric mode:

.. tab-set::

    .. tab-item:: Circular

        .. hint::

            Enabled by setting the class astrometry mode to ``'circular'`` or specifying ``'circular'`` for the astrometry mode
            in your schema file (if using CLI).

        In the circular astrometry mode, **only one uncertainty is know** for each object. We therefore model the astrometry as symmetric about
        the source.

    .. tab-item:: Axial

        .. hint::

            Enabled by setting the class astrometry mode to ``'axial'`` or specifying ``'axial'`` for the astrometry mode
            in your schema file (if using CLI).

        In the axial astrometry mode, **2 uncertainties are know** (RA and DEC) for each object. We therefore model the astrometry as an ellipse about
        the source.

    .. tab-item:: None

        .. hint::

            Enabled by setting the class astrometry mode to ``None`` (or ``'None'``/``'none'``) or specifying ``null`` for the astrometry mode
            in your schema file (if using CLI).

        This assumes that we do not have any information about the precision of the astrometry. We assume that the positioning of this object is perfect.

Selecting a Prior
+++++++++++++++++

.. warning::

    This can have a **significant** influence on the outcome of your reduction process if altered. Unless you know what you're doing,
    we highly suggest that you simply use the naive but unbiased default estimate.

Our approach is Bayesian and therefore allows the user to specify a prior :math:`P(H)`. In many cases, the prior can be estimated without
input from the user, as described in :ref:`priors` and :ref:`self-consistent-priors`.

When preparing an astrometry reduction, there are 3 options for the prior setting (:py:attr:`~reduction.AstrometricReduction.prior`).

1. ``None``: This will use the naive estimate for the prior based on the assumed optimality of the catalog and databases.

.. hint::

    The default prior isn't the optimal choice, but it is unbiased. If you get good results without a prior, you can generally
    accept them unless you have further information not accounted for in this or any other reduction process. It's much easier to
    ruin your results by specifying a flawed prior than it is to ruin them by not specifying any prior.

2. ``Iterative``: This approach will use the iterative method described in :ref:`self-consistent-priors` to determine
   priors which are self-consistent with results. This can be a pricey operation, but it will improve the overall quality of the
   prior selection.

3. ``callable``: If the user provides a callable function, that function will be used as the prior. This must be a function with signature
   ``f(match_table,catalog_row)``.

   - The ``match_table`` parameter is the match table data (or a chunk of it) on which to evaluate the prior.
   - The ``catalog_table`` parameter is the section of ``CATALOG`` containing only those objects with candidates in the ``match_table``.

   How exactly one chooses to set up their prior is there own prerogative.

Initializing
''''''''''''

In this section, we will explain how to set up your astrometry reduction process. There are two options that are of interest, either
initializing the :py:class:`reduction.ReductionProcess` from within the code or using the CLI and a reduction schema instead.

.. tab-set::

    .. tab-item:: Python



    .. tab-item:: CLI



Running the Reduction
'''''''''''''''''''''

.. raw:: html

   <hr>

.. _mathematics:
Mathematics
-----------

Astrometry
''''''''''

The most important instrumental factor is astrometry. From a general standpoint, :math:`N` separate observations of sources at
sky positions :math:`\textbf{x}_i` may or may not come from the same physical source. The likelihood that they do come from the same
source depends directly on the resolution / astrometric accuracy of the observatories.

For a source with physical position :math:`\mathbf{\mu}`, let :math:`P(\textbf{x}_i|\mathbf{\mu},I)` be the probability that
instrument I would actually detect the source at position :math:`\textbf{x}_i`. Given a set of :math:`N` points, the probability that
they are all from the same source (hypothesis :math:`H`) is

.. math::

    P(D|H) = \int_{S^2} d\mathbf{\mu}\;\; p(\mathbf{\mu}|H) \prod_{i=0}^N p(\textbf{x}_i|\mathbf{\mu},I,H).

If the observations actually originate from different sources,

.. math::

    P(D|K) = \prod_{i=0}^N \int_{S^2} d\mathbf{\mu}_i\;\; p(\mathbf{\mu}_i|K)  p(\textbf{x}_i|\mathbf{\mu}_i,I,H).

Assuming uninformative priors :math:`p(\mathbf{\mu}_i|K)` and :math:`p(\mathbf{\mu}|H)` and requiring symmetry in the distribution models,
we find

.. math::

    P(D|K) = \prod_{i=0}^N p(\mathbf{\mu}_i|K) = \frac{1}{(4\pi)^N}.

Thus, the Bayes factor is

.. math::

    \boxed{B(H,K|D) = (4\pi)^N \int_{S^2} d\mathbf{\mu}\;\; p(\mathbf{\mu}|H) \prod_{i=0}^N p(\textbf{x}_i|\mathbf{\mu},I,H) }.

Without further refinement, this is as far as we can go analytically.

Normally Distributed Error
++++++++++++++++++++++++++

Let the instrumental model take the form

.. math::

    p(\textbf{x}|\mathbf{\mu},I) = p(\textbf{x}|\mathbf{\mu},\mathbf{\Sigma}) = \frac{1}{2\pi \sqrt{|\mathbf{\Sigma}|}} \exp\left[-\frac{1}{2}(\textbf{x}-\mathbf{\mu})^T\mathbf{\Sigma}^{-1}(\textbf{x}-\mathbf{\mu})\right].

Then

.. math::

    \boxed{B(H,K|D) =  \int_{S^2} d\mathbf{\mu}\;\; \prod_{i=0}^N p(\textbf{x}_i|\mathbf{\mu},\mathbf{\Sigma}_i,H) }.

Consider the case where we only have two detections. Then

.. math::

    B(H,K|D) =  \int_{S^2} d\mathbf{\mu}\;\; \mathcal{N}(\textbf{x}_1,\mathbf{\mu},\mathbf{\Sigma}_1)\mathcal{N}(\textbf{x}_2,\mathbf{\mu},\mathbf{\Sigma}_2) = \int_{S^2} d\mathbf{\mu}\;\; \mathcal{N}(0,\mathbf{\mu},\mathbf{\Sigma}_1)\mathcal{N}(\textbf{r},\mathbf{\mu},\mathbf{\Sigma}_2),

where :math:`\textbf{r}` is the displacement between the two detections.

.. dropdown:: Mathmatics Hint: Product of Gaussians

    The integral above,

    .. math::

        \int_{S^2} d\mathbf{\mu}\;\; \mathcal{N}(0,\mathbf{\mu},\mathbf{\Sigma}_1)\mathcal{N}(\textbf{r},\mathbf{\mu},\mathbf{\Sigma}_2)

    can be solved analytically. Consider the term in the exponent:

    .. math::

        \mu^T\Sigma_1^{-1}\mu + (r-\mu)^T\Sigma_2^{-1}(r-\mu) = \mu^T\Sigma_1^{-1}\mu + r^T\sigma_2^{-1}r - 2\mu^T\Sigma_2^{-1}r + \mu^T\Sigma_2^{-1}\mu.

    Letting :math:`\Pi^{-1} = \Sigma_1^{-1} + \Sigma_2^{-1}`, one gets

    .. math::

        \mu^T\Pi^{-1}\mu - 2\mu^T\Sigma_2^{-1}r + \underbrace{r^T\sigma_2^{-1}r}_{\text{independent of $\mu$}}.

    Transforming the second term so that

    .. math::

        -2\mu^T\Sigma^{-1}_2r = -2\mu^T\Pi^{-1}\Pi\Sigma^{-1}_2r = -2\mu^T\Pi^{-1}\nu,

    where :math:`\nu = \Pi\Sigma^{-1}_2 r`.

    It now follows that

    .. math::

        \mu^T\Pi^{-1}\mu - 2\mu^T\Sigma_2^{-1}r + \underbrace{r^T\sigma_2^{-1}r}_{\text{independent of $\mu$}} = \underbrace{\mu^T\Pi^{-1}\mu - 2\mu^T\Pi^{-1}\nu + \nu^T\Pi^{-1}\nu}_{(\mu-\nu)^T\Pi^{-1}(\mu-\nu)} -  \underbrace{\nu^T\Pi^{-1}\nu + r^T\sigma_2^{-1}r}_{\text{independent of $\mu$}}

    As such, the integral

    .. math::

        \frac{1}{2\pi \sqrt{\left|\prod_i \Sigma_i\right|}} \int_{S^2} d\mathbf{\mu} \exp\left(\frac{-1}{2}\left[(x_1-\mu)^T\Sigma_{1}^{-1}(x_1-\mu) + (x_2-\mu)^T\Sigma_2^{-1}(x_2-\mu) \right]\right)\;\;

    is simply

    .. math::

        \boxed{\frac{1}{2\pi}\sqrt{\frac{|\Pi|}{\left|\prod_i \Sigma_i\right|}} \exp\left(-\frac{1}{2}\left[r^T\Sigma_2^{-1}r-\nu^T\Pi^{-1}\nu\right]\right) \;\;}.

    Furthermore,

    .. math::

        \nu^T\Pi^{-1}\nu = (\Pi\Sigma_2^{-1}r)^T\Pi^{-1}(\Pi\Sigma_2^{-1}r),

    which, recalling that :math:`\Sigma_i` are symmetric, yields

    .. math::

        r^T\Sigma_2^{-1}\Pi\Pi^{-1}\Pi\Sigma_2^{-1}r = r^T\Sigma_2^{-1}\Pi\Sigma_2^{-1}r.

    As such,

    .. math::

        -\frac{1}{2}\left[r^T\Sigma_2^{-1}r-\nu^T\Pi^{-1}\nu\right] = -\frac{1}{2}\left[r^T\Sigma_2^{-1}\left(\Pi\Pi^{-1} - \Pi\Sigma_2^{-1}\right)r\right] = -\frac{1}{2}\left[r^T\Sigma_2^{-1}\Pi\left(\Pi^{-1} - \Sigma_2^{-1}\right)r\right] =  -\frac{1}{2}\left[r^T\Sigma_2^{-1}\left[\Sigma_2^{-1}+\Sigma_1^{-1}\right]^{-1}\Sigma_1^{-1}r\right].

    Thus,

    .. math::

        \boxed{\int_{S^2} d\mathbf{\mu}\;\; \mathcal{N}(0,\mathbf{\mu},\mathbf{\Sigma}_1)\mathcal{N}(\textbf{r},\mathbf{\mu},\mathbf{\Sigma}_2) = \frac{1}{2\pi}\frac{1}{\sqrt{|\Sigma_1||\Sigma_2||\Sigma_1^{-1}+\Sigma_2^{-1}|}} \exp\left[-\frac{1}{2}\left(r^T\Sigma_1^{-1}\left(\Sigma_1^{-1}+\Sigma_2^{-1}\right)^{-1}\Sigma_2^{-1}r\right)\right]}

One finds that

.. math::

    B(H,K|D) = \boxed{\int_{S^2} d\mathbf{\mu}\;\; \mathcal{N}(0,\mathbf{\mu},\mathbf{\Sigma}_1)\mathcal{N}(\textbf{r},\mathbf{\mu},\mathbf{\Sigma}_2) = \frac{1}{2\pi}\frac{1}{\sqrt{|\Sigma_1||\Sigma_2||\Sigma_1^{-1}+\Sigma_2^{-1}|}} \exp\left[-\frac{1}{2}\left(r^T\Sigma_1^{-1}\left(\Sigma_1^{-1}+\Sigma_2^{-1}\right)^{-1}\Sigma_2^{-1}r\right)\right]}

.. important::

    We have assumed that we can make a euclidean estimate of the sphere of the sky when performing these calculations. Obviously, we integrated over points in :math:`\mathbb{R}^2`; however,
    given that points at significantly different positions on the sky cannot be associated with the same object, this estimate is sound.

.. _priors:
From Bayes Factor to Posteriori Estimates
++++++++++++++++++++++++++++++++++++++++++

Recall that the Bayes factor is defined as

.. math::

    B(H,K|D) = \frac{P(H|D)/P(H)}{P(K|D)/P(K)}

Thus, we can write the posterior as

.. math::

    P(H|D) = \frac{P(K|D)P(H)B(H,K|D)}{P(K)}.

If :math:`P(K|D) = 1-P(H|D)`, and :math:`P(K)=1-P(H)` then

.. math::

    P(H|D) = \frac{(1-P(H|D))P(H)B(H,K|D)}{1-P(H)} = \left[1+\frac{1-P(H)}{P(H)B}\right]^{-1}.

The choice of prior here is intricate and can be informed by a plethora of different pieces of information. Naively, if the catalog you are matching
against has :math:`N` objects within the search window, then the probability that a given object is the true match is just :math:`1/N`. This fails to take into account
any of the instrumental information provided and represents the most naive possible prior.

.. hint::

    Regardless of its apparent naivety, this is the default prior using in ``pyXMIP``. Users can take further control of the prior in a variety
    of ways based on their use case, including specifying it entirely.

.. _self-consistent-priors:
Self-Consistent Priors By Iterative Methods
+++++++++++++++++++++++++++++++++++++++++++

The determination of a prior can be something of a headache. In an ideal sense, where the database being matched against contains **all objects**,
then the correct prior is, very simply, :math:`1/N`; where :math:`N` is the number of objects in the search frame. Unfortunately, we know that no database
contains all objects in the universe. Furthermore, based on instrumentation, it is possible that a given source has :math:`P(H) = 0` when matched against a given catalog.

Consider a search frame around 1 source object. A database (SIMBAD for example) returns :math:`N` potential matches within a search frame area :math:`A`. Now, the expected number of matches
from the prior estimate would be

.. math::

    \left<N\right> = \sum_{i} P(H),

which should also be self-consistent with the informed estimate

.. math::

    \left<N\right> = \sum_i P_i(H|D).

Thus, one can use an iterative process to get a self-consistent (if entirely naive) estimate of the prior:

1. Start with the initial guess of 1 match in the search area: :math:`P(H) = 1/N`.
2. Calculate :math:`P_i(H|D)` for each of the potential matches using the prior you have set.
3. Calculate the effective number of matches: :math:`\left<N\right> = \sum_i P_i(H|D)`.
4. Recalculate :math:`P(H) = \left<N\right>/N = \left(\sum_i P_i(H|D)\right)/N`.
5. Continue the iteration process until a stability threshold has been accomplished.

References
----------

.. [BuSz08] Budavári, T. and Szalay, A.S., 2008. Probabilistic cross-identification of astronomical sources. The Astrophysical Journal, 679(1), p.301.
