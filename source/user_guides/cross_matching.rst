.. _cross_Referencing:
=========================
Cross Referencing
=========================

The primary use of ``pyXMIP`` is cross identification of catalog sources from known databases. On its surface, this is a relatively
simple process: query databases, identify matches, return matched catalog. In practice, there are several complications that arise:

- Is a given match a sensible match given the available information about my instrumentation?
- How likely is a given match given its distance from my detection?
- Could this match have occurred by chance?
- What is the "best" match?

``pyXMIP`` provides a framework for managing all of these different parts of the problem. In this guide, we'll discuss all of the
operational abilities of the package and describe how to fit them to your use-case.

.. raw:: html

   <hr>


Cross Matching
--------------

The first step in the cross identification process is **cross-matching**. During this stage, our goal is to find all
of the plausible object matches to the sources in your catalog. After the cross-matching process is completed, we then
undertake the **reduction** process to determine which of the identified matches is actually the **best match**.

Cross matching in ``pyXMIP`` is super simple and can be done either via the command line or via python directly!

.. hint::

    Under the hood, ``pyXMIP`` will load your catalog into a :py:class:`structures.table.SourceTable`. Thus, most common formats are are acceptable including
    ``.fits``, ``.csv``, ``.txt``, etc.

    If you experience difficulties with loading, your first stop should be assuring that the filetype you've stored your catalog in is
    actually valid.

To cross match a catalog stored in (for example) ``example.fits``, we can do either of the following:

.. tab-set::

    .. tab-item:: CLI

        For convenience, ``pyXMIP`` provides a command (automatically added to path on install) for generating cross-matching
        databases: ``pyxmip xmatch run``. The help information is as follows:

        .. code-block:: bash

            >>> pyxmip xmatch run -h
            usage: pyxmip xmatch run [-h] [-db DATABASES [DATABASES ...]] [-f] input_path output_path

            positional arguments:
              input_path            The path to the catalog you want to cross-reference
              output_path           The path to output the data to (SQL database).

            options:
              -h, --help            show this help message and exit
              -db DATABASES [DATABASES ...], --databases DATABASES [DATABASES ...]
                                    Databases to include
              -f, --force           Allow overwriting of existing data.

        As such, to compile a given catalog into a cross matching database, you need only run a command like the following.

        .. code-block:: bash

            >>> pyxmip xmatch run docs/source/examples/data/eRASS1_Hard.v1.0.fits docs/source/examples/data/cross_matched.db -f

        There are also two optional flags: ``-db`` and ``-f``. ``-f`` will allow you to overwrite existing ``.db`` output
        files with the same path. ``-db db1 db2 ... dbN`` allows you to overwrite the default set of databases to use for
        cross matching and instead use your own.

        .. hint::

            For this to work, the specified databases must be included in :py:data:`structures.databases.DEFAULT_DATABASE_REGISTRY`. If they are not in the
            :py:data:`structures.databases.DEFAULT_DATABASE_REGISTRY`, you should use a proper script to add your custom database to a new registry and then
            cross match on that registry. See :py:func:`cross_reference.cross_match` for details on using alternative registries.

        A successful run will look something like the following:

        .. code-block:: bash

            >>> pyxmip xmatch run docs/source/examples/data/eRASS1_Hard.v1.0.fits docs/source/examples/data/cross_matched.db -f
            pyXMIP : [INFO     ] 2024-04-28 13:34:12,111 X-Matching docs/source/examples/data/eRASS1_Hard.v1.0.fits into docs/source/examples/data/cross_matched.db.
            pyXMIP : [INFO     ] 2024-04-28 13:34:12,265 Cross matching with 2 databases: ['NED', 'SIMBAD'].
            pyXMIP : [WARNING  ] 2024-04-28 13:34:12,322 Table NED_MATCH exists in docs/source/examples/data/cross_matched.db. Overwrite = True -> deleting.
            pyXMIP : [WARNING  ] 2024-04-28 13:34:12,699 Table SIMBAD_MATCH exists in docs/source/examples/data/cross_matched.db. Overwrite = True -> deleting.
            pyXMIP : [INFO     ] 2024-04-28 13:34:12,750 Source matching 5466 against NED.
            pyXMIP : [DEBUG    ] 2024-04-28 13:34:12,751  [SourceTableSchema] Constructing SourceTableSchema from fits table.
            pyXMIP : [DEBUG    ] 2024-04-28 13:34:12,752  [SourceTableSchema] Failed to identify automatic match for special column Z.
            pyXMIP : [DEBUG    ] 2024-04-28 13:34:12,752  [SourceTableSchema] Failed to identify automatic match for special column TYPE.
            pyXMIP : [DEBUG    ] 2024-04-28 13:34:12,752  [SourceTableSchema] Identified special key NAME with column IAUNAME of the table.
            pyXMIP : [DEBUG    ] 2024-04-28 13:34:12,752  [SourceTableSchema] Identified special key RA with column RA of the table.
            pyXMIP : [DEBUG    ] 2024-04-28 13:34:12,752  [SourceTableSchema] Identified special key DEC with column DEC of the table.
            pyXMIP : [DEBUG    ] 2024-04-28 13:34:12,753  [SourceTableSchema] Identified special key L with column LII of the table.
            pyXMIP : [DEBUG    ] 2024-04-28 13:34:12,753  [SourceTableSchema] Identified special key B with column BII of the table.
            pyXMIP : [DEBUG    ] 2024-04-28 13:34:12,753  [SourceTableSchema] Located 2 possible coordinate frames. Selected ICRS as default.
            pyXMIP : [DEBUG    ] 2024-04-28 13:34:12,763 Querying with threading.
            pyXMIP : [INFO     ] 2024-04-28 13:37:38,003 Source matching 5466 against SIMBAD.
            pyXMIP : [DEBUG    ] 2024-04-28 13:37:38,006 Querying with threading.
            Matching from 2 databases: 100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 2/2 [05:15<00:00, 157.57s/it]


        .. hint::

            Generally databases can be queried ~ 50 times / second. Thus, you should expect runtime to scale linearly with the size of your catalog.
            A 5000 entry database takes on the order of minutes, 1,000,000 entries equates to around 2 hours.

    .. tab-item:: Python

        To accomplish this task from within a python script, we need only do the following:

        .. code-block:: python

            from pyXMIP.cross_reference import cross_match

            cross_match("input_path","output_path")

        Full details can be found at :py:func:`cross_reference.cross_match`.

What Exactly is Cross Matching?
'''''''''''''''''''''''''''''''

When you cross match a catalog in ``pyXMIP``, you're converting a :py:class:`structures.table.SourceTable` (the catalog) into
a :py:class:`cross_reference.CrossMatchDatabase` (a sql database) which contains sources matching the objects in your catalog.

Inside of the sql database, you'll find a bunch of tables (called match tables), named something like ``<database_name>_MATCH``.
These tables are just conjoined outputs from querying the specified database (using :py:meth:`structures.databases.SourceDatabase.source_match`) for
each of the objects in your catalog. Generically, the cross-matching process goes like this:

.. code-block:: python

    for database in cross_matching_databases:
        # We iterate through all of our provided databases.

        # create the blank cross-matching table.
        create_blank_cross_matching_table()

        for object in catalog:
            # search the database for matches to the catalog object.
            table = search_database(object)

            # add some extra info.
            table = add_extra_info(table)
            # write to the sql database.
            cross_matching_table += table

Obviously, its not quite that simple to actually implement, but these are the important steps to remember.


Match Tables
------------

For each database included in the cross-matching process, there will be a ``<name>_MATCH`` table in the output ``.sql`` file.
Generically, these tables can look a lot different for different databases and different query outputs.

.. hint::

    Under the hood, we query the database and join the output table to information about the object we're matching against. Thus
    you may see any number of columns from the database's query output table.

While tables can look different from one another, there are a few **standard columns** that you'll see in all cross-matching databases:

+-----------------------+----------------------------+-----------+--------------------------------------------------------+
| Column Name           | Query Schema Equivalent    | Required? | Description                                            |
+=======================+============================+===========+========================================================+
| ``NAME``              | ``schema.NAME``            | ``True``  | The name of the match candidate.                       |
+-----------------------+----------------------------+-----------+--------------------------------------------------------+
| ``RA``                | ``schema.RA``              | ``True``  | The RA of the match candidate (in degrees). Regardless |
|                       |                            |           | of the available coordinate columns, RA / DEC are      |
|                       |                            |           | always included as the base coordinate system for      |
|                       |                            |           | further manipulations.                                 |
+-----------------------+----------------------------+-----------+--------------------------------------------------------+
|``DEC``                | ``schema.DEC``             | ``True``  | The DEC of the match candidate (in degrees). See ``RA``|
+-----------------------+----------------------------+-----------+--------------------------------------------------------+
|``TYPE``               |``schema.TYPE``             | ``False`` | The object type for the candidate.                     |
+-----------------------+----------------------------+-----------+--------------------------------------------------------+
|``Z``                  |``schema.Z``                |``False``  | The redshift of the candidate.                         |
+-----------------------+----------------------------+-----------+--------------------------------------------------------+
|``CATOBJ``             |                            |``True``   | The catalog object that was matched to this candidate. |
+-----------------------+----------------------------+-----------+--------------------------------------------------------+
|``CATRA``              |                            |``True``   | The catalog object's RA (degrees).                     |
+-----------------------+----------------------------+-----------+--------------------------------------------------------+
|``CATDEC``             |                            |``True``   | The catalog declination (degrees).                     |
+-----------------------+----------------------------+-----------+--------------------------------------------------------+

The additional columns that you might find are additional (potentially useful) columns provided by your :py:class:`structures.databases.SourceDatabase` instances
used for the matching process.

.. important::

    Everytime the underlying :py:class:`structures.databases.SourceDatabase` instance is queried, it returns a :py:class:`structures.table.SourceTable`.
    That table is then put through the "cleaning process" defined by :py:attr:`structures.databases.SourceDatabase.correct_query_output`. The cleaning
    process is designed to do the **bare-minimum** to make the query writable to a ``sql`` format. After the cleaning has occured, the table is written
    to your ``sql`` database and we proceed to the next query.

    If your output table (:py:class:`structures.table.SourceTable`) is not formatted in a way which can be written to file (problematic columns,
    bad object types, etc.) then you will need to alter :py:attr:`structures.databases.SourceDatabase.correct_query_output` to get the table ready
    to write to disk.

    For built-in :py:class:`structures.databases.SourceDatabase`, we have already written these :py:attr:`structures.databases.SourceDatabase.correct_query_output`;
    thus, this should only be an issue if you are writing custom databases.

Exploring Cross Matching Outputs
'''''''''''''''''''''''''''''''''''


Once you've obtained a cross-referencing output, a lot of information can be obtained from the resulting database. Let's look at the following example:

.. code-block:: python

    >>> import pyXMIP as pyxm
    >>> catalog_table = pyxm.load("/home/ediggins/pyROSITA_test/eRASS1_Hard.v1.0.fits")

    >>> pyxm.cross_match_table(q[:10],"test.db",overwrite=True)

If we now access the ``test.db`` SQL file, we find the following:

.. code-block:: shell

    > sqlite3 test.db
    SQLite version 3.31.1 2020-01-27 19:55:54
    Enter ".help" for usage hints.

    sqlite> .schema
    CREATE TABLE IF NOT EXISTS "CATALOG" (
            "CATALOG_OBJECT" TEXT,
            "DETUID" TEXT,
            "SKYTILE" INTEGER,
            "ID_SRC" INTEGER,
            "UID" BIGINT,
            ...
            "FLAG_SP_LGA" SMALLINT,
            "FLAG_SP_GC_CONS" SMALLINT,
            "FLAG_NO_RADEC_ERR" SMALLINT,
            "FLAG_NO_EXT_ERR" SMALLINT,
            "FLAG_NO_CTS_ERR" SMALLINT,
            "FLAG_OPT" SMALLINT
    );
    CREATE TABLE IF NOT EXISTS "META" (
            "PROCESS" TEXT,
            "TABLE" TEXT,
            "DATE_RUN" TEXT
    );
    CREATE TABLE IF NOT EXISTS "NED_STD_MATCH" (
            "No." BIGINT,
            "NAME" TEXT,
            "RA" FLOAT,
            "DEC" FLOAT,
            "TYPE" TEXT,
            "Velocity" FLOAT,
            "Redshift" FLOAT,
            "Redshift Flag" TEXT,
            "Magnitude and Filter" TEXT,
            "Separation" FLOAT,
            "References" BIGINT,
            "Notes" BIGINT,
            "Photometry Points" BIGINT,
            "Positions" BIGINT,
            "Redshift Points" BIGINT,
            "Diameter Points" BIGINT,
            "Associations" BIGINT,
            "CATOBJ" TEXT,
            "CATRA" FLOAT,
            "CATDEC" FLOAT
    );
    CREATE TABLE IF NOT EXISTS "SIMBAD_STD_MATCH" (
            "NAME" TEXT,
            "RA" TEXT,
            "DEC" TEXT,
            "RA_PREC" BIGINT,
            "DEC_PREC" BIGINT,
            "COO_ERR_MAJA" FLOAT,
            "COO_ERR_MINA" FLOAT,
            "COO_ERR_ANGLE" BIGINT,
            "COO_QUAL" TEXT,
            "COO_WAVELENGTH" TEXT,
            "COO_BIBCODE" TEXT,
            "TYPE" TEXT,
            "RA_d_A" FLOAT,
            "DEC_d_D" FLOAT,
            "SCRIPT_NUMBER_ID" BIGINT,
            "CATOBJ" TEXT,
            "CATRA" FLOAT,
            "CATDEC" FLOAT
    );

As you can see, we now have the ``CATALOG`` and the ``META`` table, and the standard columns are all placed in the tables.

.. raw:: html

   <hr>

Cross Match Database
--------------------

The :py:class:`cross_reference.CrossMatchDatabase` class is one of the key classes in ``pyXMIP``. It provides a number
of powerful methods and attributes for successfully interacting with / reducing cross-matching data directly from external databases.

In this document, we will describe the use of these class instances and demonstrate how to make the most out of their
backend.

What is a ``CrossMatchDatabase``?
''''''''''''''''''''''''''''''''''

Once a catalog has been referenced against a set of specified databases, the result is a ``SQL`` database containing all of the
cross matching data for each of the databases in the set and for each of the sources in the catalog. These ``SQL`` databases
can be difficult to effectively interact with, particularly for those without great experience using ``SQL``. Fortunately, in
``pyXMIP``, the :py:class:`cross_reference.CrossMatchDatabase` plays the critical role of accessing, organizing, updating, and
analyzing these ``SQL`` databases all from within the python environment.

Making ``CrossMatchingDatabase`` Instances
++++++++++++++++++++++++++++++++++++++++++

As was mentioned above, :py:class:`cross_reference.CrossMatchDatabase` instances are effectively ``SQL`` representations in ``python``. They are made
when the user cross-matches against databases. In many cases, such as when using functions like :py:func:`cross_reference.cross_match` or :py:func:`cross_reference.cross_match_table`
the returned object is already a :py:class:`cross_reference.CrossMatchDatabase` instance. In other cases, a cross-matching database can be loaded directly from it's filepath:

.. code-block:: python

    >>> import pyXMIP as pyxmip
    >>> cmd = pyxmip.CrossMatchDatabase("docs/source/examples/data/cross_matched.db")
    >>> print(cmd)
    <CrossMatchDatabase @ docs/examples/data/cross_matched.db>

In this case, we've just opened a cross-matching database created by cross-matching eROSITA eRASS1 data against NED and SIMBAD.

.. raw:: html

   <hr>

Database Alterations
'''''''''''''''''''''

When cross-matching is performed, ``pyXMIP`` will return a :py:class:`cross_reference.CrossMatchDatabase` which has been corrected
slightly to standardize things like column names, object types, etc. Each of these processes, along with a variety of other protocols are
examples of "database alterations." There are a number of contexts in which these can be useful for your analysis, and so we will present them here.

The ``META`` table
++++++++++++++++++

Let's begin by calling the ``cmd.tables`` attribute and see what comes up:

.. code-block:: python

    >>> cmd.tables
    ['CATALOG', 'META', 'NED_MATCH', 'SIMBAD_MATCH']

As one might guess, these are the tables inside of the underlying ``SQL`` database.

- ``NED_MATCH`` and ``SIMBAD_MATCH`` are the **raw match tables** for each of the databases.

  - These contain all of the match candidates for each of the catalog sources.

- ``CATALOG`` is precisely what it sounds like: the original catalog being used for the cross-matching.

- ``META`` is special: it contains a log of all of the things that have been done to this cross-matching database.

Taking a look at the ``META`` table, we can see all of the various "processes" which have been run:

.. code-block:: python

    >>> cmd.meta
                PROCESS         TABLE                  DATE_RUN
    0  CATALOG_INCLUDED           all  Sun Apr 28 18:50:23 2024
    1    OBJECT_CORRECT     NED_MATCH  Sun Apr 28 18:50:27 2024
    2    OBJECT_CORRECT  SIMBAD_MATCH  Sun Apr 28 18:50:27 2024
    3    COLUMN_CORRECT     NED_MATCH  Sun Apr 28 18:50:27 2024
    4    COLUMN_CORRECT  SIMBAD_MATCH  Sun Apr 28 18:50:27 2024

Each line represents a different action that was taken on a specific table.

- ``CATALOG_INCLUDED`` indicates that the original catalog has been written to the sql database.
- ``OBJECT_CORRECT`` indicates that the ``TYPE`` column of each table has been converted to the SIMBAD type conventions.
- ``COLUMN_CORRECT`` indicates that the columns have been renamed to meet our standard schema.

Under the hood, each of these is just a method of the :py:class:`cross_reference.CrossMatchDatabase` database. There are a number
of these processes, which can be read about in detail below.

.. rubric:: Available CMD Processes

+----------------------------------+--------------------------------------------------------------------+---------------------------------------------------------------------------------------+
| META TAG                         | Method                                                             | Description                                                                           |
+==================================+====================================================================+=======================================================================================+
| ``CATALOG_INCLUDED``             | :py:meth:`cross_reference.CrossMatchDatabase.add_catalog`          | Add a catalog (any readable source table) to the database. In most cases, this is     |
|                                  |                                                                    | done automatically during the cross-matching process.                                 |
+----------------------------------+--------------------------------------------------------------------+---------------------------------------------------------------------------------------+
| ``CORRECT_OBJECT_TYPES``         | :py:meth:`cross_reference.CrossMatchDatabase.correct_object_types` | Correct the object types so that they match the ``SIMBAD`` object types.              |
+----------------------------------+--------------------------------------------------------------------+---------------------------------------------------------------------------------------+
| ``CORRECT_COLUMN_NAMES``         | :py:meth:`cross_reference.CrossMatchDatabase.correct_column_names` | Correct the names of columns to fit the ``pyXMIP`` standard.                          |
+----------------------------------+--------------------------------------------------------------------+---------------------------------------------------------------------------------------+
| ``CORRECT_COORDINATE_COLUMNS``   | :py:meth:`cross_reference.CrossMatchDatabase.correct_coordinates`  | Correct the available coordinate columns.                                             |
+----------------------------------+--------------------------------------------------------------------+---------------------------------------------------------------------------------------+

.. important::

    Later in this article, we discuss reduction processes, which are simply a special subset of these processes acting in
    :py:class:`cross_reference.CrossMatchDatabase` instances.

.. raw:: html

   <hr>


Interacting with Raw Data
'''''''''''''''''''''''''

.. hint::

    For a comprehensive list on all methods and attributes associated with cross match databases, check out the API documentation (:py:class:`cross_reference.CrossMatchDatabase`).


Now that we've demonstrated how to obtain / open :py:class:`cross_reference.CrossMatchDatabase` instances, we can start exploring what
they're capable of. It's important to remember that there are 3 key features of :py:class:`cross_reference.CrossMatchDatabase` classes:

1. Ability to interact directly with the cross matching data.
2. Ability to organize, sort, and standardize the available data.
3. Ability to perform statistical reduction processes on the data available.

In this section, we will introduce the methods that allow :py:class:`cross_reference.CrossMatchDatabase` to accomplish these tasks.


Pull Matches from Tables
++++++++++++++++++++++++

Let's start with the most basic information you can pull from :py:class:`cross_reference.CrossMatchDatabase`. Because these classes wrap
underlying SQL data, we are often interested in what match tables are in the file:

.. code-block:: python

    >>> print(cmd.match_tables)
    ['NED_MATCH', 'SIMBAD_MATCH']

You can use :py:attr:`cross_reference.CrossMatchDatabase.match_tables` to obtain only ``MATCH`` data instead of any other tables.

Let's try to access the ``MATCH`` data. To pull a table, you can simply index into the database like it's a dictionary:

.. code-block:: python

    >>> ned_match_table = cmd['NED_MATCH']
    >>> ned_match_table
                     CATALOG_OBJECT  CATALOG_RA  ...                Object Name   Type
    0       1eRASS J013729.2-195637   24.371790  ...  WISEA J013725.22-195648.6  |IrS|
    1       1eRASS J013729.2-195637   24.371790  ...  WISEA J013726.30-195650.1  |IrS|
    2       1eRASS J013729.2-195637   24.371790  ...  WISEA J013726.61-195716.5  |IrS|
    3       1eRASS J013729.2-195637   24.371790  ...  WISEA J013726.73-195620.2    |G|
    4       1eRASS J013729.2-195637   24.371790  ...  WISEA J013726.91-195636.1  |IrS|
                             ...         ...  ...                        ...    ...
    183741  1eRASS J062052.1-284050   95.217384  ...  WISEA J062052.79-284119.1  |IrS|
    183742  1eRASS J062052.1-284050   95.217384  ...  WISEA J062053.30-284146.9  |IrS|
    183743  1eRASS J062052.1-284050   95.217384  ...  WISEA J062053.43-284054.5  |IrS|
    183744  1eRASS J062052.1-284050   95.217384  ...  WISEA J062053.89-284012.8  |IrS|
    183745  1eRASS J062052.1-284050   95.217384  ...  WISEA J062055.65-284021.4  |IrS|
    [183746 rows x 7 columns]

There are any number of things that can now be done with this data.

Plot A Match
++++++++++++

One useful utility provided by the :py:class:`cross_reference.CrossMatchDatabase` class is the ability to visualize the specific
matches for a given object. This is done with the :py:meth:`cross_reference.CrossMatchDatabase.plot_matches` method. As an example, let's look at
matches in SIMBAD to the Bullet Cluster.

Let's start by looking to see if there are any matching detections to the JWST deep field cluster SMACS J0723.3-7327:

.. code-block:: python

    >>> cmd.query("SELECT * FROM SIMBAD_MATCH WHERE NAME == 'SMACS J0723.3-7327'")
                CATALOG_OBJECT  ...                                               TYPE
    0  1eRASS J072316.9-732718  ...  |gLe|gLe|ClG|ClG|ClG|ClG|ClG|ClG|C?G|C?G|C?G|X|X|
    [1 rows x 7 columns]

We see that the object is matched only to the eRASS1 catalog object 1eRASS J072316.9-732718.

In this case, we've used DESI and HST images as the background. We can do similar things from within ``pyXMIP`` without any external software!

Let's go ahead and plot the matches like this:

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> from matplotlib.colors import SymLogNorm
    >>> q.plot_matches("1eRASS J072316.9-732718", "NED_MATCH",
                   norm=SymLogNorm(linthresh=1e-1),
                   cmap=plt.cm.gnuplot,
                   scatter_kwargs=dict(c='w',marker='o',s=5),
                   fov='2 arcmin',
                   resolution=1000,
                   hips_kwargs=dict(hips_path='CDS/P/JWST/F444W'))
    >>> plt.show()

.. image:: ../images/plots/cross_match_database_2.png

We can compare this output with the corresponding JWST composite image as well:

.. image:: https://blogs-images.forbes.com/startswithabang/files/2017/01/abell_lens-1.jpg

The white points in our image are the cross-matched reference points detected for this cluster!


.. raw:: html

   <hr>

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

Below, you can find information about each of the available reduction processes.

- :ref:`instrumental_reduction`
- :ref:`object_type_reduction`

Further Information
'''''''''''''''''''

- :ref:`reduction_schema`
- :ref:`custom_reductions`
