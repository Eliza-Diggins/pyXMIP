.. _cross_Referencing:
=========================
Cross Referencing
=========================

Cross Matching
--------------

The first step in the cross identification process is **cross-matching**. During this stage, our goal is to find all
of the plausible object matches to the sources in your catalog. After the cross-matching process is completed, we then
undertake the **reduction** process to determine which of the identified matches is actually the **best match**.

Cross matching in ``pyXMIP`` is super simple and can be done either via the command line or via python directly!

.. hint::

    Under the hood, ``pyXMIP`` will load your catalog into an ``astropy`` table. Thus, most common formats are are acceptable including
    ``.fits``, ``.csv``, ``.txt``, etc.

    If you experience difficulties with loading, your first stop should be assuring that the filetype you've stored your catalog in is
    actually valid.

To cross match a catalog stored in (for example) ``example.fits``, we can do either of the following:

.. tab-set::

    .. tab-item:: CLI

        To use the command line interface (CLI), we can simply use the 1-line command

        .. code-block:: bash

            >>> pyxmip xmatch run "path/to/example.fits" "output_path.db"

        There are also two optional flags: ``-db`` and ``-f``. ``-f`` will allow you to overwrite existing ``.db`` output
        files with the same path. ``-db db1 db2 ... dbN`` allows you to overwrite the default set of databases to use for
        cross matching and instead use your own.

        .. hint::

            For this to work, the specified databases must be included in ``DEFAULT_DATABASE_REGISTRY``.

    .. tab-item:: Python

        To accomplish this task from within a python script, we need only do the following:

        .. code-block:: python

            from pyXMIP.cross_reference import cross_match

            cross_match("input_path","output_path")

        Full details can be found at :py:func:`cross_reference.cross_match`.

Exploring Cross Referencing Outputs
'''''''''''''''''''''''''''''''''''

Once you've obtained a cross-referencing output, a lot of information can be obtained from the resulting database.



Match Reduction
---------------

The source matching process is relatively naive; it simply samples from the specified databases subject to your chosen
parameters but doesn't pay any attention to further "common-sense" decisions that could be used to improve the fidelity
of the matching process. Furthermore, during the search, all sources within the search radius are added, not just the best match.

As such, the next step in the reduction process is to run the match reduction algorithm, which takes the large database
of identified sources and their match in your catalog and determines which matches are legitimate and which are spurious.
This is a complex process, and can be controlled significantly by the user. In this section, we will provide an overview of the
user's options for this process.

Mathematical Overview
'''''''''''''''''''''

The reduction process is based on a cost minimization framework determined by each of the user's selected sub-processes.
Effectively, each subprocess run will determine a "cost" for each possible match to a given source. Exactly how that cost is calculated is specific
to the particular sub-process; however, the value is always in the interval :math:`[0,1]`, where 0 indicates a perfect match.

In the reduction schema file, the user may assign weights to each of the sub-processes to fine-tune the minimization process. In general,
each potential match is assigned a **collective cost** :math:`C(x|y)`, representing the total cost across all sub-processes of matching source
:math:`y` to database object :math:`x`. Each sub-process has a weight :math:`c_i(x|y)` and

.. math::

    C(x|y) = \sum_{i} \alpha_i c_i(x|y),

where :math:`\alpha_i` is the user assigned weight for the :math:`i` th sub-process. Generically, the :math:`\alpha_i` should be reflective of the user's
confidence in the success of the given sub-processes' model given the available information.

Sub-Process Overview
''''''''''''''''''''


.. tab-set::

    .. tab-item:: Instrumental Matching

        .. note::

            Instrumental matching is implemented using the ``RUN_PARAMS`` option ``INSTRUMENTAL=true`` and then adding the ``INSTRUMENT_PARAMS``
            heading to the schema file for your reduction run.

        .. caution::

            This is currently an in-development feature of the pyXMIP package. Ongoing development is underway to improve the quality and physical sophistication
            of the PSF modeling used here.

        Instrument specific physics can help to inform how reasonable a given source is as a match to the given catalog object. Effectively, this comes down
        to determining whether or not a given source falls within the PSF of the telescope that was used to create the user's catalog. As such, we model
        the cost function of a given match :math:`x` to a catalog object :math:`y` as a 0 covariance Gaussian with a user-specified :math:`\sigma`. Thus,

        .. math::

            c_{\mathrm{instrument}}(x|y) = 1 - \sqrt{2\pi \sigma^2} N(\delta_{xy},\sigma),

        where :math:`\delta_{xy}` is the angular distance between the catalog object and the potential match object.


    .. tab-item:: Object Type Matching

        .. note::

            Object type matching is implemented using the ``RUN_PARAMS`` option ``OTYPES=true`` and then adding ``OBJECT_PARAMS`` to
            the schema file for your reduction run.

        The intention of object matching is to constrain the types of objects that can be permitted and those which cannot be. The pyXMIP algorithms
        all use the standardized SIMBAD convention for object names / definitions. If the database being referenced has a well defined schema with
        object type mapping enabled, then these are automatically converted to the correct SIMBAD type. In the reduction schema file, the user may
        then include a mapping of object types to cost values from 0 to 1 depending on the relative feasability of that object type.

        .. hint::

            SIMBAD type objects are defined hierarchically. If only a parent classification is found in the schema, then all of the subtypes
            inherit the same cost.

    .. tab-item:: Poisson Modeling

        .. note::

            Poisson modeling is implemented using the ``RUN_PARAMS`` option ``POISSON=true`` and then adding the ``POISSON_PARAMS``
            heading to the schema file for your reduction run.

        Depending on the nature of the database against which you are cross-referencing, there may be an abundance of spurious matches from
        a particularly prevalent type of object such as stars or background galaxies. This is partially managed by the Object Type matching protocol; however,
        additional robustness can be enforced using Poisson modeling, which comes with the added benefit of permitting matches to potentially poorly described
        database objects.

        .. hint::

            A particularly common example of this is in the NED database, which contains sources classified as ``IrS`` corresponding to 2MASS, WISE, and SDSS
            observations. These objects may be stars, galaxies, or other objects. As such, it may be too restrictive to simply dismiss the entire class of
            potential matches, but one still wants to avoid spuriously including them.

        Poisson mapping utilizes the pre-constructed object poisson-maps for a given database against which the user is cross-referencing. If you are not already
        familiar with the pyXMIP poisson mapping system, see :ref:`poisson-mapping` for a more comprehensive description. In essence, Poisson maps are generated prior
        to the cross-matching process by randomly sampling the database at points on the sky and performing a density estimation algorithm to determine the
        frequency at which a particular type of object in the specific database appears on the sky. The probability of then finding such an object within an angular circle
        of radius :math:`r` is

        .. math::

            P(k>0|r,\lambda) = 1 - \exp\left(-\lambda(\phi,\theta) \pi r^2 \right).

        This can be used to determine the probability that a match candidate :math:`x` is a spurious match (occuring due to random chance) as

        .. math::

            P_{\mathrm{spurious}}(x|y) = 1 - \exp\left(-\lambda(\phi_y,\theta_y) \pi \delta_{xy}^2 \right)

        This is directly a proxy for the cost function of this procedure.

Reduction Schema
''''''''''''''''

The cornerstone of the reduction process is the schema file, which dictates the runtime behavior of the reduction algorithm. In this section, we will describe the
layout of these files and how one can go about constructing them.

General Formatting
++++++++++++++++++

Like all schema files in the pyXMIP ecosystem, the reduction schema is a ``.yaml`` file with the following required headers:

- ``RUN_PARAMS``: Core settings for the reduction run. This includes selecting what additional data should be considered,
  which subprocesses should be enabled, etc.
- ``IO_PARAMS``: IO related parameters, including file paths and other information.

In addition to these two required sets of parameters, the following may also be specified / enabled

- ``POISSON_PARAMS``: Parameters for the Poisson mapping sub-process.
- ``OBJECT_PARAMS``: Parameters for the object type sub-process.
- ``INSTRUMENT_PARAMS``: Parameters for the instrument specific sub-process.

.. tab-set::

    .. tab-item:: RUN_PARAMS

        .. csv-table:: RUN_PARAM table
            :class: longtable
            :align: center
            :width: 100%
            :widths: 3,3,10,1
            :file: _tables/run_params.csv
            :header-rows: 1

    .. tab-item:: IO_PARAMS

        .. csv-table:: Table Title
            :file: _tables/run_params.csv
            :header-rows: 1
