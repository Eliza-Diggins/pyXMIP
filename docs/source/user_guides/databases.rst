.. _databases:
==========
Databases
==========

pyXMIP's core functionality is the cross-mapping / cross-referencing of source data against known databases of astronomical objects.
The :py:class:`structures.databases.SourceDatabase` class is the corner-stone of this functionality. On this page, we'll describe
how these objects work and how you can utilize them to meet your scientific needs.

The :py:class:`~structures.databases.SourceDatabase` class and it's subclasses are designed (in a general sense), to allow for easy
interaction between the user and underlying data sources in the process of performing cross-matching operations. As such, the goal of these
classes is to

1. Simplify the process of querying databases for objects; particularly in the case of remote databases.
2. Provide an interface for generating critical data for reduction and statistical inference during the cross-matching process.
3. Provide intuitive methods for cross-matching catalogs against these databases.

To do that, a variety of different databases are implemented as described in this document.

Databases At A Glance
---------------------

As was stated above, databases are all about interacting with and standardizing catalog data (in various forms). To do this, databases are represented
by **class-instances**. For example, a user can access the IPAC-Caltech NED database simply by creating an instance of the corresponding remote database class:

.. code-block:: python

    >>> from pyXMIP import NED
    >>> ned_instance = NED()
    >>> print(ned_instance.name)
    NED_STD

These instances of database classes have a couple of universally shared properties:

- Each instance has a :py:meth:`structures.databases.SourceDatabase.query_radius` method which allows the user to select all of the
  sources within a given radius of a particular coordinate.
- Each instance has a source-matching method (specific to the type of database) which matches sources from *another* catalog against
  the given database.
- Each instance has a couple of important parameters:

  - :py:attr:`structures.databases.SourceDatabase.name`: The name of the instance.
  - :py:attr:`structures.databases.SourceDatabase.query_schema`: The :py:class:`schema.SourceTableSchema` corresponding to the tables returned when querying the database.
  - :py:attr:`structures.databases.SourceDatabase.query_config`: Configuration settings for the query's behavior.

Querying Databases
++++++++++++++++++

The simplest operation that can be performed by a database is **querying**.

.. admonition:: definition

    In the context of ``pyXMIP``, **querying** a database the operation which takes the sky position and angular radius of a
    search area as input and returns all of the objects within the search area in the database.

To query a database, the :py:meth:`structures.databases.SourceDatabase.query_radius` method is used:

.. code-block:: python

    >>> from pyXMIP import NED
    >>> from astropy import units
    >>> from astropy.coordinates import SkyCoord

    # Load the NED database class (IPAC/NED extragalactic database)
    >>> ned_instance = NED()

    # Create query position and radius
    >>> gc = SkyCoord(0,0,unit='deg',frame='galactic') # The galactic center.
    >>> radius = 0.1*units.arcmin # the search radius.

    # Query
    >>> query_output = ned_instance.query_radius(gc,radius)
    pyXMIP : [INFO     ] 2024-05-05 13:54:59,485 Querying <SkyCoord (Galactic): (l, b) in deg
        (0., 0.)> in NED_STD...
    >>> for object in query_output['Object Name']:
            print(f"Object {object} is in the query region!")
    Object 2MASS J17453717-2856090 is in the query region!
    Object [TSK2008] 0223 is in the query region!


Poisson Atlases
+++++++++++++++

.. note::

    Poisson Atlases are an involved topic in and of themselves; for a more detailed description of the processes involved in these
    structures, see :ref:`poisson_mapping`.

One of the useful tools provided by ``pyXMIP`` when cross-referencing / cross-matching catalogs is the ability to model the likelihood of
a spurious match given the type of astronomical object. To achieve this, ``pyXMIP`` utilizes database-linked structures called :py:class:`structures.map.PoissonAtlas` to
model the density of object types on the sky.

.. hint::

    An example of this would be the likelihood that there is a star within distance 0.21 arcmin from our source. If this
    probability is high, then it's not surprising of we find a matching star when cross-referencing.

While providing comprehensive information about Poisson atlases is outside of the scope of this document, it is worth being aware that
**all databases** can have an attached Poisson atlas, which is used to perform these computations. This is stored in the :py:attr:`structures.databases.SourceDatabase.poisson_atlas`
attribute:

.. code-block:: python

    >>> psn_map = ned_instance.poisson_atlas
    >>> print(psn_map.path)
    pyXMIP/bin/psn_maps/NED.poisson.fits

As you can see, there is a built-in Poisson atlas for NED. You can also add points directly:

.. code-block:: python

    >>> print(len(psn_map.COUNTS))
    1790
    >>> ned_instance.add_sources_to_poisson(10,1*units.arcmin)
    100%|█████████████████████████████████████████| 10/10 [00:08<00:00,  1.11it/s]
    pyXMIP : [INFO     ] 2024-05-05 14:31:07,043 Adding data to the Poisson map at /home/ediggins/pyXs/pyXMIP/bin/psn_maps/NED.poisson.fits.
    >>> print(len(psn_map.COUNTS))
    1800

Cross Matching
++++++++++++++

.. important::

    Cross matching is a very important functionality; so much so that there is an entirely separate page for it: :ref:`cross_Referencing`.

Types of Databases
------------------

To begin, it should be understood that there are **2 types** of databases in ``pyXMIP``: the :py:class:`structures.databases.LocalDatabase` and the
:py:class:`structures.databases.RemoteDatabase` classes.

.. note::

    Unlike local databases, which are generated directly from the :py:class:`structures.databases.LocalDatabase`, most
    remote databases are actually subclasses of the :py:class:`structures.databases.RemoteDatabase` class.

To a certain extent, the differences between these classes are obvious from the name: local databases allow the user to convert a catalog of sources into a searchable
database, remote databases allow access to resources like NED, SIMBAD, and Vizier. Because the necessary structures to interact with each of these two types of data
are vastly different, the two classes operate very distinctly.

.. important::

    The differences between these two classes is important to understand when performing analyses with ``pyXMIP``. Despite our best
    efforts to create a very uniform experience between the two, the vastly different necessary operations has led to some potential gotchas
    if the user doesn't understand the difference.

Local Databases
+++++++++++++++

A local database is a database with is loaded directly from a dataset on the user's local computer. Generally, these are smaller and contained in some
standard table format (``.fits``, ``.csv``, etc.) In this document, we're going to use the eROSITA Hard Band data as an example.
You can find it `here <https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Hard.tar.gz>`_.

.. important::

    :py:class:`structures.databases.LocalDatabase` is really just additional structure on top of the :py:class:`structures.table.SourceTable` class.

In order to load data as a local database, you can simply use the following lines of code:

.. code-block:: python

    >>> tb = SourceTable.read("/home/ediggins/pyROSITA_test/eRASS1_Hard.v1.0.fits")
    >>> db = LocalDatabase(tb,"example_database")
    >>> print(db)
    pyXMIP : [DEBUG    ] 2024-05-05 17:17:00,721  [SourceTableSchema] Constructing SourceTableSchema from fits table.
    pyXMIP : [DEBUG    ] 2024-05-05 17:17:00,722  [SourceTableSchema] Failed to identify automatic match for special column Z.
    pyXMIP : [DEBUG    ] 2024-05-05 17:17:00,722  [SourceTableSchema] Failed to identify automatic match for special column TYPE.
    pyXMIP : [DEBUG    ] 2024-05-05 17:17:00,722  [SourceTableSchema] Identified special key NAME with column IAUNAME of the table.
    pyXMIP : [DEBUG    ] 2024-05-05 17:17:00,722  [SourceTableSchema] Identified special key RA with column RA of the table.
    pyXMIP : [DEBUG    ] 2024-05-05 17:17:00,722  [SourceTableSchema] Identified special key DEC with column DEC of the table.
    pyXMIP : [DEBUG    ] 2024-05-05 17:17:00,722  [SourceTableSchema] Identified special key L with column LII of the table.
    pyXMIP : [DEBUG    ] 2024-05-05 17:17:00,723  [SourceTableSchema] Identified special key B with column BII of the table.
    pyXMIP : [DEBUG    ] 2024-05-05 17:17:00,723  [SourceTableSchema] Located 2 possible coordinate frames. Selected ICRS as default.
    <LocalDatabase example_database>

.. note::

    The ``DEBUG`` statements that appear here are ``pyXMIP`` constructing a :py:class:`schema.SourceTableSchema` for the :py:class:`structures.table.SourceTable` instance
    that's being used.

There are some **very important** things to know about local databases:

- Local databases typically don't have default attribute values. Thus,

  - :py:attr:`~structures.databases.LocalDatabase.name` must be set when the instance is initialized.
  - :py:attr:`~structures.databases.LocalDatabase.query_schema` is **automatically** set to the schema of the initializing table.
  - :py:attr:`~structures.databases.LocalDatabase.correct_query_output` is the same standard function as :py:attr:`~structures.databases.SourceDatabase.correct_query_output`.
  - :py:attr:`~structures.databases.LocalDatabase.query_config` is ``None``. The query is trivial and thus doesn't have a special configuration to worry about.

- Local databases don't need to utilize remote queries, which means that **parallelism is not implemented**.

  - Note that many methods still take a ``parallel_kwargs`` keyword, but they are only present for generality of the code.

Remote Databases
++++++++++++++++

Complementary to the :py:class:`~structures.databases.LocalDatabase` class, it is also possible to configure remote databases like NED, SIMBAD, and Vizier.
These classes tend to be marginally more complex due to the API interactions involved.

Just like :py:class:`~structures.databases.LocalDatabase`, :py:class:`~structures.databases.RemoteDatabase` can be queried using :py:meth:`structures.databases.RemoteDatabase.query_radius`, and
they share many of the same parameters.

.. hint::

    Many of the most popular databases are built-in as subclasses of :py:class:`~structures.databases.RemoteDatabase`. If a desired database is not already
    written as a class, then the user many need to write a new database class.

Remote databases utilize all of the various defaults:

- :py:attr:`~structures.databases.LocalDatabase.name`: If the database is being used in default configuration, it will be ``<DB_NAME>_STD``, otherwise, it will be required like it is for local databases.
- :py:attr:`~structures.databases.LocalDatabase.query_schema` is set at the class level as a default. It almost never needs to be overridden.
- :py:attr:`~structures.databases.LocalDatabase.correct_query_output` is a purpose built method to standardize format.
- :py:attr:`~structures.databases.LocalDatabase.query_config` is ``None`` is a very useful tool for configuring the details of a given databases outputs.

Just like local databases, remote databases are **instances of a given database class**; for example, you can access the NED database via the built-in class
:py:class:`structures.databases.NED`.

.. hint::

    **What's the difference between the class and the instance?**

    For :py:class:`structures.databases.LocalDatabase`, the underlying class isn't very useful;
    you have to provide a table (during instantiation) and interact with the instance, not the class. For :py:class:`structures.databases.RemoteDatabase`, the picture
    isn't quite as simple; the underlying class represents a **default configuration** of the database while instances can represent various modifications of the underlying settings.

    Nonetheless, **you always need to work with instances, not classes**!

Let's look at the example shown at the top of this page:

.. code-block:: python

    >>> from pyXMIP import NED
    >>> from astropy import units
    >>> from astropy.coordinates import SkyCoord

    # Load the NED database class (IPAC/NED extragalactic database)
    >>> ned_instance = NED()

    # Create query position and radius
    >>> gc = SkyCoord(0,0,unit='deg',frame='galactic') # The galactic center.
    >>> radius = 0.1*units.arcmin # the search radius.

    # Query
    >>> query_output = ned_instance.query_radius(gc,radius)
    pyXMIP : [INFO     ] 2024-05-05 13:54:59,485 Querying <SkyCoord (Galactic): (l, b) in deg
        (0., 0.)> in NED_STD...
    >>> for object in query_output['Object Name']:
            print(f"Object {object} is in the query region!")
    Object 2MASS J17453717-2856090 is in the query region!
    Object [TSK2008] 0223 is in the query region!

Here, we see that the ``ned_instance`` object is a **blank instance** of the :py:class:`structures.databases.NED` database class. This is
the **default** instantiation of NED; but many things might be worth changing!

Unlike local databases, all of the built-in remote databases have **built-in Poisson Atlases**! These can be accessed directly by
calling the :py:attr:`structures.databases.RemoteDatabase.poisson_atlas`:

.. code-block:: python

    >>> import matplotlib.pyplot as plt

    >>> db = NED()
    >>> poisson_atlas = db.poisson_atlas

    >>> map = ps.get_map("IrS") # Fetch the IR sources PSN map.
    >>> map.plot(cmap='gnuplot')
    >>> plt.show()

.. image:: ../images/plots/databases_1.png

Differences to Know
+++++++++++++++++++

Now that we've gone through the two critical types of databases, it's worth summarizing the similarities
and differences between the two.

+---------------+--------------------------------------------------+--------------------------------------------------+
| Difference    | :py:class:`structures.databases.LocalDatabase`   | :py:class:`structures.databases.RemoteDatabase`  |
+===============+==================================================+==================================================+
| Poisson Maps  | Generated by individual instances, not built-in  | Classes have built-in Poisson Atlases.           |
+---------------+--------------------------------------------------+--------------------------------------------------+
| Schema        | Deduced from the source table used to build the  | Schema is either built-in or provided by hand.   |
|               | database.                                        |                                                  |
+---------------+--------------------------------------------------+--------------------------------------------------+
| Correction    | Specified by custom prescription if necessary.   | Generally built-in, can be custom built if needed|
+---------------+--------------------------------------------------+--------------------------------------------------+
| Query Config  | Not needed                                       | Optional, needed to make alterations to query.   |
+---------------+--------------------------------------------------+--------------------------------------------------+


Creating Custom Databases
-------------------------

In many cases, a combination of :py:class:`structures.databases.LocalDatabase` and the built-in databases are sufficient
for most cross matching needs; however, it may be the case that a new :py:class:`structures.databases.RemoteDatabase` needs to
be implemented for your use case or a specialized subclass of :py:class:`structures.databases.LocalDatabase` is needed.

In general, writing such custom classes is a non-trivial activity; however, we have put together the resources below to help
users achieve this goal with minimal stress.

Local Databases
+++++++++++++++

Because :py:class:`structures.databases.LocalDatabases` are rather simple classes, we do not provide a
template for generating subclasses. In all cases, the ``__init__()`` method needs to have the following form

.. code-block:: python

    def __init__(self, table, db_name, **kwargs):

        # Default behavior from standard LocalDatabase
        self.table = table
        super().__init__(db_name, query_schema=self.table.schema, **kwargs)

        # Additional custom code can go here.

There are only 3 methods in :py:class:`structures.databases.LocalDatabase`:

- ``source_match`` and ``source_match_memory`` are the core methods for cross-matching and, if edited, should
  follow the same conventions as the base class.
- ``_query_radius`` must also be present and, if altered, must maintain its call signature.

Beyond these 3 methods, you have total freedom over custom class implementations. You can simply
subclass :py:class:`structures.databases.LocalDatabase` and add new methods as need be.


Remote Databases
++++++++++++++++

Remote databases are far more complex to implement than their local counterparts. For remote databases, we encourage you to
work off of the following class template:

.. code-block:: python

    class _RemoteTemplate(RemoteDatabase):
    """
    Template of a custom RemoteDatabase class.
    """

    default_poisson_atlas_path = os.path.join(poisson_map_directory, "<class_name>.poisson.fits")
    default_query_config = {} #--> THIS SHOULD NEVER CHANGE. ADD _DatabaseConfigSetting's instead.
    default_query_schema = DEFAULT_SOURCE_SCHEMA_REGISTRY["<class_name>"]

    # -- Settings -- #
    query_config_setting_1 = _databaseConfigSetting(default='default')
    query_config_setting_2 = _databaseConfigSetting(default='default')

    def __init__(self, name="<class_name>_STD", **kwargs):
        super().__init__(name, **kwargs)

        self.config() # This can be renamed.

    def config_ned(self):

        # This is where you should take all of the _databaseConfigSetting's you've got implemented
        # and use them to actually modify the behavior of your querying system.

        # When you're using an Astro-query based system, you can usually just change the class parameters in that package.

    def default_correct_query_output(self, table):

        # Make alterations as necessary to standardize format #

        table = super().default_correct_query_output(table)

        return table

    def _query_radius(self, position, radius):
        """
        Query the remote database at the specified position and pull all sources within a given radius.

        Parameters
        ----------
        position: :py:class:`astropy.coordinates.SkyCoord`
            The position at which to query.
        radius: :py:class:`astropy.units.Quantity`
            The angular area about which to query.

        Returns
        -------
        :py:class:`astropy.table.Table`
        """
        # -- Attempt the query -- #
        try:
            output = SourceTable(your_custom_query_code(position,radius))
        except requests.exceptions.ConnectionError:
            raise DatabaseError(
                f"Failed to complete query [{position},{radius}] to NED due to timeout."
            )

        # -- return data if valid -- #
        output.schema = self.query_schema
        return output

Additional methods can be added as you see fit. We also encourage you to read the Developer Information section of this
page.

---

Developer Information
---------------------

For developers working on this project, there is some additional, in-depth information
worth keeping in mind when writing code.

Query Schema
++++++++++++

For remote databases, there should **always** be a custom :py:class:`schema.SourceTableSchema` associated with the output
from the default configuration of that database. If users want to alter defaults in such a way as to break that functionality,
they will need to provide their own version of the schema. Built-in schema should be included in the ``/bin/builtin_schema` directory
and follow the standard naming convention used there.

.. hint::

    While you can get away without adding all of the possible detail to a given schema, we highly encourage
    creating comprehensive / complete schema.

Query Configurations
++++++++++++++++++++

In the template shown above, we include the ``self.config()`` method. In cases where you are writing a custom class, this should
be used to modify the underlying query backend to account for the configuration settings. Each configuration setting is a
descriptor instance (``_databaseConfigSetting``) which looks up its value in the ``self.query_config`` dictionary.

Output Standardization
++++++++++++++++++++++

One of the typical issues faced when developing database classes is that the output table is not generally in a state
which meets our formatting requirements for use with other databases.
