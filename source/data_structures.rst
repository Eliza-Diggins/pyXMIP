.. _data_structures:
Data Structures
===============

DATABASES
---------

SCHEMAS
-------

Schemas (subclasses of the :py:class:`schema.Schema` abstract base class) are objects which are attached to various
collection classes (databases, catalogs, etc.) to provide those collection classes with the information they need to
access the data they are given. There are a variety of types of schema that are used throughout the pyXMIP package which the user
should be aware of.

Sky Collection Schema
+++++++++++++++++++++

Collection schema (:py:class:`schema.SourceTableSchema`) are the most important schema class in pyXMIP. These schema connect the
raw data in a table (underlying a catalog or a database) to the corresponding physical properties of interest. In practice, the
most important job of :py:class:`schema.SourceTableSchema` instances is to connect the available ``columns`` of the underlying data table
to standardized conventions used by the class for cross matching or other tasks.

Formatting
``````````

- ``column_map`` is the most important key of the :py:class:`schema.SourceTableSchema`. This specifies the names of special columns in your catalog.

  - You **DO NOT** need to include every column in the ``column_map``. Only special columns have to be specified.
  - Coordinate systems are determined from the ``column_map``. You must have at least one complete coordinate system in your ``column_map``.
  - The following are the available special columns:

  .. table::

    =================== ============================================================
    Special Column Name Description
    =================== ============================================================
    ``Z``               The column representing the redshift. (optional)
    ``RA``              The column containing RA coordinates. (optional)
    ``DEC``             The column containing DEC coordinates. (optional)
    ``b``               The column containing galactic latitudes (optional)
    ``l``               The column containing galactic longitudes (optional)
    =================== ============================================================

- ``default_coord_system`` is the default coordinate system to apply to your data. In general, catalogs may contain more than 1 complete
  coordinate system and so this column is used to determine which frame to use by default. It can be changed on the fly by setting :py:attr:`schema.SourceTableSchema.coordinate_system` .

The following is the standard format for the :py:class:`schema.SourceTableSchema` class.

.. code-block:: yaml

    #-----------------------------------------#
    # MAPPINGS
    #-----------------------------------------#
    column_map: # The column map tells the schema what different columns in your catalog actually are.
        RA: "<name of RA column in catalog>"
        Z: "<name of the redshift column in catalog>"

    #-----------------------------------------#
    # DEFAULTS
    #-----------------------------------------#
    default_coord_system: "ICRS" # The default coordinate system to use.
