.. _source_data:
===============================================
Loading Source Data
===============================================

The first step to successfully using ``pyXMIP`` for your scientific purpose is actually loading data into the environment.
For the most part, this task should be quite simple; however, there are a few things to keep in mind when doing this that
should be known. In this reference page, we provide a comprehensive description of loading data into the ``pyXMIP`` environment
and the possible interactions that can be had therein.

Loading a Catalog from Disk
---------------------------

``pyXMIP`` is designed from cross-referencing / identification of astronomical objects within catalogs. To load a catalog into
the library, simply use :py:func:`~structures.table.load`. In this reference, we're going to use the eROSITA eRASS 1
all sky survey data located `here <https://erosita.mpe.mpg.de/dr1/erodat/>`_.

.. code-block:: python

    >>> import pyXMIP as pyxmip
    >>> catalog = pyxmip.load("eRASS1_Hard.v1.0.fits")
    >>> print(catalog)
    <SourceTable length=5466>
            IAUNAME                      DETUID              ... FLAG_OPT
                                                             ...
            bytes23                     bytes32              ...  int16
    ----------------------- -------------------------------- ... --------
    1eRASS J002937.2-310209 eb01_009120_020_ML00004_002_c010 ...        0
    1eRASS J004207.0-283154 eb01_012120_020_ML00002_002_c010 ...        0
    1eRASS J004922.3-293108 eb01_012120_020_ML00003_002_c010 ...        0
    1eRASS J005448.9-311230 eb01_012120_020_ML00004_002_c010 ...        0
    1eRASS J012910.8-214156 eb01_021111_020_ML00001_002_c010 ...        0
    1eRASS J012338.1-231059 eb01_021114_020_ML00001_002_c010 ...        0
    1eRASS J013729.2-195637 eb01_024111_020_ML00005_003_c010 ...        0
    1eRASS J015219.0-183235 eb01_027108_020_ML00003_002_c010 ...        0
    1eRASS J015721.4-161420 eb01_029105_020_ML00003_002_c010 ...        0
                        ...                              ... ...      ...
    1eRASS J061245.9-520045 em01_095141_020_ML01024_003_c010 ...        0
    1eRASS J130424.5+132131 em01_197078_020_ML00270_003_c010 ...        0
    1eRASS J130833.7+032916 em01_197087_020_ML00223_003_c010 ...        0
    1eRASS J131142.7-303349 em01_197120_020_ML00303_002_c010 ...        0
    1eRASS J060719.0-715042 em01_095162_020_ML01380_002_c010 ...        0
    1eRASS J061241.0-705749 em01_095162_020_ML02377_002_c010 ...        0
    1eRASS J062748.9-230934 em01_096114_020_ML00478_002_c010 ...        0
    1eRASS J062429.5-351724 em01_096126_020_ML00282_002_c010 ...        0
    1eRASS J131252.8-161006 em01_198105_020_ML00374_002_c010 ...        0
    1eRASS J130904.4-260824 em01_198117_020_ML00487_002_c010 ...        0

Loading a catalog results in a class called :py:class:`~structures.table.SourceTable`. This is a *slightly modified* version of
the typical :py:class:`astropy.table.table.Table` class with some special purpose components. We can do all of the classic modifications
to the table that would be possible in :py:mod:`astropy`.

Just like any other table of data, we can plot the sources in this catalog:

.. code-block::

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from astropy.units import Quantity
    >>> ra,dec = Quantity(catalog['RA']).to_value('rad')-np.pi,Quantity(catalog["DEC"]).to_value('rad')

    >>> fig = plt.figure(figsize=(12,8))
    >>> ax = fig.add_subplot(111,projection='aitoff')
    >>> ax.scatter(ra,dec,c=np.log10(catalog['ML_RATE_0']))
    >>> plt.show()

.. image:: ../images/plots/source_data_1.png

Working With Catalogs
---------------------

Before moving on to other, more complex, aspects of ``pyXMIP``, it's worth describing some additional details of the
:py:class:`structures.table.SourceTable` class.

Table Schema
''''''''''''

The first thing to familiarize yourself with is the schema associated with the table. You can find a comprehensive guide to
schemas in ``pyXMIP`` on the :ref:`schema` page as well as on the API documentation of the :py:mod:`schema` module. For now,
we'll introduce only the broader aspects of the schema.

The :py:class:`schema.SourceTableSchema` class is associated with every table and effectively tells ``pyXMIP`` what all of the
different column names and associated data actually mean. These schema located the coordinates, object types, redshifts, and other
key information about the catalog. They also determine the base coordinate system and are necessary for much of the useful operations
provided by this software.

In many cases, ``pyXMIP`` can deduce some parts of the schema based on the table itself. Let's try this with the eROSITA data.
To access the schema, we simply use the :py:attr:`~structures.table.SourceTable.schema` attribute:

.. code-block:: python

    >>> schema = catalog.schema
    pyXMIP : [DEBUG    ] 2024-04-27 21:13:12,308  [SourceTableSchema] Constructing SourceTableSchema from fits table.
    pyXMIP : [DEBUG    ] 2024-04-27 21:13:12,318  [SourceTableSchema] Failed to identify automatic match for special column Z.
    pyXMIP : [DEBUG    ] 2024-04-27 21:13:12,318  [SourceTableSchema] Failed to identify automatic match for special column TYPE.
    pyXMIP : [DEBUG    ] 2024-04-27 21:13:12,318  [SourceTableSchema] Identified special key NAME with column IAUNAME of the table.
    pyXMIP : [DEBUG    ] 2024-04-27 21:13:12,318  [SourceTableSchema] Identified special key RA with column RA of the table.
    pyXMIP : [DEBUG    ] 2024-04-27 21:13:12,318  [SourceTableSchema] Identified special key DEC with column DEC of the table.
    pyXMIP : [DEBUG    ] 2024-04-27 21:13:12,318  [SourceTableSchema] Identified special key L with column LII of the table.
    pyXMIP : [DEBUG    ] 2024-04-27 21:13:12,318  [SourceTableSchema] Identified special key B with column BII of the table.
    pyXMIP : [DEBUG    ] 2024-04-27 21:13:12,320  [SourceTableSchema] Located 2 possible coordinate frames. Selected ICRS as default.

Whoa! That's a lot of outputs, what do they actually tell us?

The :py:class:`schema.SourceTableSchema` class will first try to identify special columns (things like redshift, coordinates, object name, etc.) from
the available column names. According to the log, we have failed to find a redshift or object type (eRASS data doesn't have these), but we did find
the ``IAUNAME`` column as the ``NAME`` of each object, along with ``RA``, ``DEC``, ``L``, ``B``; each of which is a coordinate column. Finally,
we note that two viable coordinate systems were identified and that ICRS (RA/DEC) was chosen as the default.

How does this actually help us? **Now pyXMIP can identify important information automatically**!

The schema allows the table to locate all of the critical information for us. Let's see an example of how this works: Let's extract the
source names from the table:

.. code-block:: python

    >>> print(catalog.NAME)
            IAUNAME
    -----------------------
    1eRASS J002937.2-310209
    1eRASS J004207.0-283154
    1eRASS J004922.3-293108
    1eRASS J005448.9-311230
    1eRASS J012910.8-214156
    1eRASS J012338.1-231059
    1eRASS J013729.2-195637
    1eRASS J015219.0-183235
    1eRASS J015721.4-161420
    1eRASS J020037.2-164905
                        ...
    1eRASS J062733.8-513408
    1eRASS J061245.9-520045
    1eRASS J130424.5+132131
    1eRASS J130833.7+032916
    1eRASS J131142.7-303349
    1eRASS J060719.0-715042
    1eRASS J061241.0-705749
    1eRASS J062748.9-230934
    1eRASS J062429.5-351724
    1eRASS J131252.8-161006
    1eRASS J130904.4-260824
    Length = 5466 rows

You can also easily access the coordinates, redshift, object type, etc. through the same attributes that would represent
those columns in the schema class. It is highly recommended that new users familiarize themselves with the available attributes
as described in the :py:class:`structures.table.SourceTable` documentation.
