Density Mapping
===============

One of the core objectives of pyXMIP is to cross match novel sources in astronomical surveys against databases of known
objects. One of the challenges presented is that, for any given sky position, there is always a chance that a database contains an entirely
unrelated source at or near the same sky position. This will preferentially bias source types which are more common in a given database and
bias databases with more sources over those with fewer sources. To improve our methodology, we use statistical methods to model the probability
of a spurious match in some region of sky. This process is roughly broken into 4 steps:

1. Sample sources from the matching database.
2. Construct a :py:class:`stats.density.DensityModel` and :py:class:`stats.density.PopulationDensity` map of the different source populations in the matching database.
3. Cross match a local catalog (unknown sources) against the matching database.
4. Model the probability of a given match given the density of database sources from the map.

On this page, we present the basic ideas behind this process.

Statistical Methods
===================

Consider a collection of sources (preferably a self-consistent population of sources) from the matching database.

.. admonition:: Example

    We could be looking to map the density of QSO objects in NED.

We begin by randomly sampling from that population and constructing a :py:class:`stats.density.PopulationDensity` object with a collection of source positions on the sky and
the number of counts within a given solid angle around that point. Each :py:class:`stats.density.PopulationDensity` instance is effectively a collection of sky positions
:math:`\textbf{r} = (\phi_i,\theta_i)`, sampling areas :math:`A_i` and source counts :math:`k_i`.

.. note::

    TODO: this should link to info about the random sampling methodology.

Once a :py:class:`stats.density.PopulationDensity` object has been produced for the given database or catalog, we then seek to fit a given
density model to the data. A density model (encapsulated by the :py:class:`stats.density.DensityModel` class) is a function

.. math::

    f(\phi,\theta;\mathbf{\Omega}) \;\text{such that}\; \left<k_{\mathrm{matches}}\right>_{\phi,\theta} = f(\phi,\theta; \mathbf{\Omega})

for some set of fit parameters :math:`\mathbf{\Omega}`. As such, :math:`f(\phi,\theta;\mathbf{\Omega})` effectively encapsulates the density
of sources at a given position on the sky.

This information can then be used in the cross matching process. If a source at :math:`\textbf{r}` on the sky is matched to a database source at :math:`\textbf{r}'` which
has some angular distance :math:`\delta` from the true source, the probability of finding such a source randomly is

.. math::

    1-P(k=0|\phi,\theta,\delta) = 1 - \lambda(\phi,\theta,\delta) \exp\left(-\lambda(\phi,\theta,\delta)\right)

where

.. math::

    \lambda(\phi,\theta,\delta) = f(\phi,\theta;\mathbf{\Omega}_{\mathrm{fit}}) \pi^2 \delta.

From this, we may determine our confidence in the legitimacy of a particular match.

.. note::

    There are two cruxes in this process:

    1. Computing a reasonable density fit.
    2. Constructing the density maps.

    For common databases, the density maps for each of the object types are updated regularly and available directly from
    the repository.

.. raw:: html

   <hr style="color:black">


Fitting Density Models
----------------------

For a given population density sample consisting of :math:`N` points :math:`\mathbf{\phi},\mathbf{\theta}`, their associated
source counts :math:`\textbf{k}`, and the sampling area :math:`\textbf{A}`, the liklihood of obtaining that data from a Poisson process is

.. math::

  p(\textbf{k}|\mathbf{\phi},\mathbf{\theta},\mathbf{A}) = \prod_{i=1}^N \mathrm{Poisson}\left(k_i,\Lambda_i)\right)

The (negative) log-likelihood associated with the distribution is

.. math::

  \mathcal{L} = - \log p(\textbf{k}|\mathbf{\phi},\mathbf{\theta},\mathbf{A}) = \sum_{i=1}^N \Lambda_i - k_i\log \Lambda_i - \log k_i!,

where

.. math::

    \Lambda_i = A_if(\phi_i,\theta_i,\mathbf{\Omega}).

To this point, all of the various density estimation methods are homogeneous. From here, there are several possible classes of fit that
the user may choose to employ:

.. tab-set::

    .. tab-item:: Bayesian Estimation

       If the user elects to use a Bayesian framework for the model fit, then


    .. tab-item:: Maximum Likelihood Estimation

      If the user elects to use maximum likelihood estimation (MLE) to determine the fit, the backend will attempt to find a minima of :math:`\mathcal{L}` corresponding
      to the effective best fit parameters. Generically, this problem takes the form

      .. math::

        \partial_\mu \mathcal{L} =  \Lambda_{i,\mu}\left(1 - \frac{k_i}{\Lambda_i}\right) = A_if_{i,\mu}\left(1 - \frac{k_i}{A_i f_i}\right) = 0

      Letting

      .. math::

        \textbf{F}_i^j(\mathbf{\Omega}) = f_{i,\Omega_j}(\Omega), \;\text{and}\; \textbf{U}_i = 1 - k_i A_i^{-1} f_{i}^{-1}(\mathbf{\Omega}),

      The problem reduces to a non-linear optimization problem of the form

      .. math::

        \textbf{A} \cdot \textbf{F}(\mathbf{\Omega})U(\mathbf{\Omega}) = 0

      which contains :math:`|\mathbf{\Omega}|` variables in :math:`N` equations. In almost all cases, these equations should be tractable. For the most part,
      we accomplish the minimization by numerical method; however, some special cases have analytically tractable solutions.

      The covariance matrix for the model fit is (by the `Cramer-Rao Bound <https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Rao_bound>`_)

      .. math::

        \mathrm{Cov}\mathbf{\Omega}_{ij} = \left(- \mathcal{L}_{,ij}(\mathbf{\Omega}_0)\right)^{-1/2}


Built-in Models
================

.. tab-set::

    .. tab-item:: Constant Density

        The constant density model assumes (for a given probability model :math:`P(N_E|A,f)`) a distribution function :math:`f` of the form

        .. math::

            f(\phi,\theta;\alpha) = \alpha.

        This is, in effect, the simplest possible density model. For a given population, the log-likelihood goes as

        .. math::

            -\mathcal{L}(\Phi,\Theta;\alpha) = -\sum_{i}^N \log P(N_i|A_i,\alpha).

        **Maximum Likelihood Estimation**:

        In the case of a Poisson probability model (the most typical option), the log-likelihood becomes

        .. math::

            -\mathcal{L}(\Phi,\Theta;\alpha) = \sum_{i}^N A_i\alpha - N_i \log A_i \alpha + \log N_i!.

        This takes on a maxima for

        .. math::

            -\partial_\alpha \mathcal{L} = \sum_i^N \left[ A_i - \frac{N_i}{\alpha} \right] = 0.

        Thus, we find the MLE estimator for :math:`\alpha` is

        .. math::

            \boxed{\hat{\alpha}_{\mathrm{MLE}} = \frac{\sum N_i}{\sum A_i}}

        Additionally, the second derivative is

        .. math::

            \left.-\partial_\alpha^2 \mathcal{L}\right|_{\alpha = \hat{\alpha}} = \sum_i^N \frac{N_i}{\hat{\alpha}^2}.

        The variance is then

        .. math::

            \boxed{ \hat{\sigma}^2_\alpha = \frac{\sum_i^N N_i}{\left(\sum_i^N A_i\right)^2}.}
