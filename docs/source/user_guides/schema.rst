.. _schema:
========================
pyXMIP Schema Classes
========================

``pyXMIP`` uses classes (called ``Schema`` classes) to translate the various conventions and datatypes in the user's data and
the standardized versions that can be used in the backend. This includes doing everything from identifying the coordinate
columns in a table to determining what types of objects are present in a given database and everything in between.

In this guide, the :py:mod:`schema` module is described in detail with all of the critical knowledge to utilize ``Schema``
classes successfully.

.. raw:: html

   <hr>

What is a Schema
----------------

There are a variety of different types of schema in ``pyXMIP``, the most common of which is :py:class:`schema.SourceTableSchema`.
Schema have a variety of purposes, but they all do effectively the same thing: they allow ``pyXMIP`` to convert your data
to a form it understands.

Under-the-hood, schemas are **dataclasses** built on top of the ``pydantic`` ``BaseModel`` class. Effectively, this means that
schema are just classes with a bunch of attributes you can access as you would any other class in python. What's different is that,
from the developer's perspective, there's quite a bit of additional logic to assure that datatypes, values, etc. are actually valid.
``pydantic`` provides structure (very powerful structure) for these sorts of data validation tasks and we thus use it here.

Let's begin by listing some of the most important schema:

- :py:class:`schema.SourceTableSchema` - By far the most important of the schemas.

  - Everytime the user loads a table in ``pyXMIP``, it get an associated :py:class:`schema.SourceTableSchema` instance.

    - Sometimes, this can be deduced from the available column names (i.e. RA and DEC columns); however, in other cases
      it cannot.
    - ``pyXMIP`` will attempt to generate a schema with access to all of the correct parameters; however, if this fails it may
      be necessary to write your own.

  - The corresponding :py:class:`schema.SourceTableSchema` allows ``pyXMIP`` to convert column names to standardized forms, understand
    different object types, and manage a variety of settings.

- :py:class:`schema.ReductionSchema`

  - The reduction process converts raw match data to a statistically optimized cross-reference catalogs. To do this, the user generally
    writes a **reduction file** (``yaml``) containing all of the settings for that reduction process. This is read into ``pyXMIP`` and
    becomes an instance of the :py:class:`schema.ReductionSchema` class.

.. hint::

    Generally, if you're writing a ``.yaml`` file to feed into ``pyXMIP``, **you're writing a schema file**. The conventions for how to
    specify different data types are standardized between different types of schema, so you only need to get familiar with one syntax.

Writing Schema
--------------

The first step to working with ``Schema`` classes is understanding what they contain and how they work. In this section, we'll begin
by describing how to write a schema file before moving on to interacting with these objects within ``pyXMIP``.

.. important::

    **Every** schema is also a ``BaseModel`` from ``pydantic``. As such, you can initialize it with keyword arguments or provide
    a dictionary of the necessary data and, while loading the model, ``pyXMIP`` will validate that all of the data you provided
    makes sense.

The exact format, including what entries are allowed and what headers should contain, is determined by the type of schema
file you're writing, but there are some conventions which are maintained across all schema.


Schema Types
++++++++++++++

There are a variety of different types of schema with different formats. The tabs below contain the available schema types
and their conventions for formatting.


.. tab-set::

    .. tab-item:: Source Table Schema

        .. hint::

            For a comprehensive listing of the available options / settings recognized by these schema, look at the
            API documentation: :py:class:`schema.SourceTableSchema`. It may also be worthwhile to look at various examples
            in the documentation where these schema are interacted with.

        Source table schema are the most important schema in ``pyXMIP``. They are necessary whenever you need ``pyXMIP`` to
        recognize table structures from your own data.

        A schema can be written to / read from disk in a number of formats; in this document, we use ``.yaml``, but ``.json``, and ``.toml`` are
        also recognized. To write a :py:class:`schema.SourceTableSchema`, you need 3 sections of your ``.yaml`` files as follows:

        .. code-block:: yaml

            # Example yaml SourceTableSchema
            column_map:
                # The column map encodes special columns in your catalog table.
                TYPE:
                  name: "my_type_column"
                Z:
                  name: "my_redshift_column"
                NAME:
                  name: "my_object_id"
                RA:
                  name: "my_RA"
                  unit: 'deg'
                DEC:
                  name: "my_DEC"
                  unit: 'deg'

            # ... More special column definitions...
            object_map:
                # The object map converts your object types to SIMBAD / pyXMIP object types
                star: "*"
                g_cluster: "GClstr"

            # ... All of your object types ...

            settings:
                default_coord_system: "ICRS"

        For comprehensive details, see :py:class:`schema.SourceTableSchema`.

        Some special schemas are already built-in. In particular, **all** of the built-in remote databases already have
        schema in the ``pyXMIP`` environment. These default schema are stored in registries (:py:class:`schema.SchemaRegistry`).
        In particular, the ``schema.DEFAULT_SOURCE_SCHEMA_REGISTRY`` contains the list of built-in table schema available.

        For example, we can get a list of these registries as follows:

        .. code-block:: python

            >>> from pyXMIP.schema import DEFAULT_SOURCE_SCHEMA_REGISTRY
            >>> print(f"The available schema are {DEFAULT_SOURCE_SCHEMA_REGISTRY.as_list()}")
            The available schema are ['pyXMIP/bin/builtin_schema/source_table/SIMBAD.yaml', 'pyXMIP/bin/builtin_schema/source_table/NED.yaml']

            >>> ned_schema = DEFAULT_SOURCE_SCHEMA_REGISTRY['NED']
            >>> print(ned_schema.column_map.dict())
            {'RA': {'name': 'RA', 'unit': 'deg'},
             'DEC': {'name': 'DEC', 'unit': 'deg'},
             'RA_ERR': None,
             'DEC_ERR': None,
             'GAL_L': None,
             'GAL_B': None,
             'GAL_L_ERR': None,
             'GAL_B_ERR': None,
             'NAME': {'name': 'Object Name', 'unit': 'None'},
             'Z': {'name': 'Redshift', 'unit': 'None'},
             'TYPE': {'name': 'Type', 'unit': 'None'}}


        The user can also create these registries and / or store alter the existing ones to load their own set of default schema.




    .. tab-item:: Reduction Schema

        .. hint::

            For a comprehensive listing of the available options / settings recognized by these schema, look at the
            API documentation: :py:class:`schema.ReductionSchema`.

        Reduction schema tell ``pyXMIP`` how to perform the statistical analyses of interest on an existing cross-matching database.
        The structure of the schema will largely depend on the exact nature of the reductions being performed. Users can write their own
        reductions and register the parameters for those reductions so that they can be entered into the schema.

        Because reduction schema are an involved topic, users will find formatting options at


Interaction With Schema
-----------------------

Now that we've introduce the :py:class:`Schema` class, we can now dive into the functionality of these objects. For the most
part, schema behave like a regular data class. You can interact with their attributes, alter their attributes, etc.
For example, if we want to see all of the attributes in the schema, we need only use

.. code-block:: python

    >>> from pyXMIP.schema import DEFAULT_SOURCE_SCHEMA_REGISTRY
    >>> ned_schema = DEFAULT_SOURCE_SCHEMA_REGISTRY['NED']
    >>> print(ned_schema.__fields__.keys())
    dict_keys(['column_map', 'object_map', 'settings'])

We can also iterate through the column map just like we would anything else:

.. code-block:: python

    >>> for k,v in (NED.default_query_schema.column_map.__fields__.items()):
            print(k,getattr(NED.default_query_schema.column_map,k))
    RA name='RA' unit=Unit("deg")
    DEC name='DEC' unit=Unit("deg")
    RA_ERR None
    DEC_ERR None
    GAL_L None
    GAL_B None
    GAL_L_ERR None
    GAL_B_ERR None
    NAME name='Object Name' unit=None
    Z name='Redshift' unit=None
    TYPE name='Type' unit=None

In addition to the typical dictionary-like behaviors, schema also have type-specific methods and properties to
manage various aspects of their use-case.

.. tab-set::

    .. tab-item:: Source Table Schema

        The key purpose of :py:class:`schema.SourceTableSchema` instances is to allow ``pyXMIP`` to translate between
        your table format and the ``pyXMIP`` backend. As such, there are a lot of methods which are of little use to
        the user, but are utilized within the :py:class:`structures.table.SourceTable` backend.

        More usefully, all of the different settings and ``column_map`` special columns are immediately available as
        properties.

        .. code-block:: python

            >>> print(f"The NED database TYPE column is {ned_schema.TYPE}.")
            'The NED database TYPE column is Type.'

        In many cases, some special columns don't actually appear in the ``column_map``. These can still be accessed and
        set from within the ``pyXMIP`` environment:

        .. code-block:: python

            >>> print(ned_schema.column_map)
            {'RA': 'RA',
             'DEC': 'DEC',
             'NAME': 'Object Name',
             'TYPE': 'Type',
             'Z': 'Redshift'}

            >>> print(ned_schema.L)
            None

            >>> ned_schema.L = "Galactic Longitude"

            >>> print(ned_schema.L)
            'Galactic Longitude'

            >>> print(ned_schema.column_map)
            {'RA': 'RA',
            'DEC': 'DEC',
            'NAME': 'Object Name',
            'TYPE': 'Type',
            'Z': 'Redshift',
            'L': 'Galactic Longitude'}

        One of the most important functionalities if the :py:class:`schema.SourceTableSchema` class is managing coordinate systems.
        It is critical that ``pyXMIP`` is able to identify the locations of the catalog items and convert them accurately too and
        from the correct systems. To do this, the :py:class:`schema.SourceTableSchema` has a several convenient methods.

        The first of these methods is :py:meth:`schema.SourceTableSchema.coordinate_system`, which returns the columns and class
        representing the current coordinate system:

        .. code-block:: python

            >>> ned_schema.coordinate_system
            (astropy.coordinates.builtin_frames.icrs.ICRS, ['RA', 'DEC'])

        This can also be changed to utilize a different coordinate system!

        .. hint::

            In order for a coordinate system change to be valid, you must actually have corresponding columns in your
            table.

        .. code-block:: python

            >>> ned_schema.available_coordinate_frames()
                {astropy.coordinates.builtin_frames.icrs.ICRS: ['RA', 'DEC']}

            >>> ned_schema.coordinate_system = 'Galactic'
                Traceback (most recent call last):
                  File "<ipython-input-16-073ab00b423e>", line 1, in <module>
                    ned_schema.coordinate_system = 'Galactic'
                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                  File "pyXMIP/schema.py", line 455, in coordinate_system
                    getattr(astro_coords, value) in self.available_coordinate_frames()
                AssertionError: The coordinate frame Galactic is not a valid coordinate frame for this schema.

    .. tab-item:: Reduction Schema

        .. important::

            Section in progress...




Automated Schema
----------------

.. important::

    Section in progress...
