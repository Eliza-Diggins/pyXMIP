"""
Module for interactions with external and local databases.

Notes
-----

For more information on databases, see :ref:`databases`.

"""
import os
import threading
import time
import warnings
from abc import ABC, abstractmethod
from itertools import repeat
from typing import Any, Callable, Generic, Type, TypeVar

import numpy as np
import requests.exceptions
import sqlalchemy as sql
from astropy import units
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table, vstack
from astroquery.ipac.ned import Ned
from astroquery.simbad import Simbad
from tqdm.auto import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from pyXMIP.schema import DEFAULT_SOURCE_SCHEMA_REGISTRY, SourceTableSchema
from pyXMIP.structures.map import PoissonAtlas
from pyXMIP.structures.table import SourceTable, correct_column_types
from pyXMIP.utilities.core import bin_directory, enforce_units
from pyXMIP.utilities.logging import mainlog
from pyXMIP.utilities.types import Registry, convert_np_type_to_sql

poisson_map_directory: str = os.path.join(bin_directory, "psn_maps")

# -- Typing Variables -- #
Instance = TypeVar("Instance")
Value = TypeVar("Value")
Attribute = TypeVar("Attribute")


class DatabaseError(Exception):
    """
    Collective error type for issues with database queries.
    """

    def __init__(self, message: str = None):
        self.message = message
        super().__init__(self.message)


class _DatabaseConfigSetting(Generic[Instance, Attribute, Value]):
    # Descriptor class for query configuration settings.
    # Stored in query_config.
    def __init__(self, default: Any = None):
        self.default = default

    def __set_name__(self, owner: Type[Instance], name: str) -> None:
        self._name = name

    def __get__(self, instance: Instance, owner: Type[Instance]) -> Any:
        if self._name in instance.query_config:
            return instance.query_config[self._name]
        else:
            return self.default

    def __set__(self, instance: Instance, value: Any) -> None:
        instance.query_config[self._name] = value


class SourceDatabase(ABC):
    """
    Abstract class representation of a database class. All other database types are subclasses of this class.


    Notes
    -----

    .. warning::

        This is an abstract class; there is no reason to instantiate this class.

    """

    # -- class default variables and methods -- #
    default_poisson_atlas_path: str = os.path.join(poisson_map_directory, "NONE")
    """str: The path to the default poisson atlas file.

    For remote databases, this is generally non-trivial. For local databases,this is typically None.
    """
    default_query_config: dict[str, Any] = None
    """dict: The configuration information for performing the query.

    Depending on the particular database, this may vary substantially."""
    default_query_schema: SourceTableSchema = None
    """:py:class:`schema.SourceTableSchema`: The schema for the returned table after a query.

    For local databases, this is automatically set to the base catalog. For remotes, it must generally be constructed.
    """
    _thread_lock = threading.Lock()

    def __init__(self, db_name: str, *args, **kwargs):
        """
        Initialize the :py:class:`SourceDatabase`.

        Parameters
        ----------
        db_name: str
            The name of the database instance.
        poisson_path: str
            The path to the Poisson atlas file.
        query_config: dict
            Configuration settings for the querying process.
        query_schema: :py:class:`schema.SourceTableSchema`
            The schema to associated with query outputs.
        correct_query_output: callable
            The callable for correcting formatting of the query output.
        """
        self.name: str = db_name
        """str: The name of this database instance.
        """
        self._poisson_atlas_path: str = kwargs.pop(
            "poisson_path", self.__class__.default_poisson_atlas_path
        )
        """str: The path to the Poisson Atlas file.

        By default, it is ``None``; but remote subclasses should have a built-in poisson path.
        """
        self.query_config: dict = kwargs.pop(
            "query_config", self.__class__.default_query_config
        )
        """dict: The query configuration options.

        By default, it is ``None``; this varies for subclasses.
        """
        self.query_schema: SourceTableSchema = kwargs.pop(
            "query_schema", self.__class__.default_query_schema
        )
        """:py:class:`schema.SourceTableSchema`: The schema for the returned table after a query.

        For local databases, this is automatically set to the base catalog. For remotes, it must generally be constructed.
        """
        self.correct_query_output: Callable[[SourceTable], SourceTable] = kwargs.pop(
            "correct_query_output", self.default_correct_query_output
        )
        """callable: The method for altering returned source tables to be writable to sql.

        By default, these are pre-set. See the relevant documentation.
        """

    def __str__(self):
        return f"<SourceDatabase {self.name}>"

    def __repr__(self):
        return self.__str__()

    def register(self, registry: Registry):
        """Add this instance to a DBRegistry."""
        if self.name not in registry:
            registry[self.name] = self

    # -- Querying -- #

    @abstractmethod
    def _query_radius(self, position: SkyCoord, radius: units.Quantity) -> SourceTable:
        """
        DEVELOPERS: This is the query radius that should be altered in sub-classes.
        """
        pass

    def query_radius(self, position: SkyCoord, radius: units.Quantity) -> SourceTable:
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
        :py:class:`structures.table.SourceTable`

        Notes
        -----

        .. admonition:: Dev-Note

            This method doesn't need to be altered when writing a custom database class. It simply provides a standardized
            wrapper for the query. The ``_query_radius`` method contains the code which varies from class to class.

        Examples
        --------

        **Query for NED objects within 1 arcmin of the galactic center**:

        .. code-block:: python

            >>> from pyXMIP import NED
            >>> from astropy.coordinates import SkyCoord
            >>> from astropy import units

            >>> database_instance = NED()
            >>> query = database_instance.query_radius(SkyCoord(0,0,unit='deg',frame='galactic'),radius=1*units.arcmin)

            >>> print(f"The first matching object is {query[0]['Object Name']}, Separation={query[0]['Separation']} arcmin.")
            The first matching object is 2MASS J17453276-2856166, Separation=0.981 arcmin.

        """

        radius = enforce_units(radius, units.arcmin)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mainlog.info(
                f"Querying lat={np.round(position.frame.spherical.lat.deg,decimals=3)}, lon={np.round(position.frame.spherical.lon.deg,decimals=3)} in {self.name}..."
            )
            output_table = self._query_radius(
                position,
                radius,
            )

        return output_table

    @classmethod
    def _default_correct_query_output(
        cls, table: SourceTable | Table, schema: SourceTableSchema = None
    ) -> SourceTable | Table:
        """
        The default query output correction callable.

        This is only used if :py:attr:`SourceTable.correct_query_output` is not specified during initialization.

        Parameters
        ----------
        table: :py:class:`astropy.table.table.Table`
            The table being corrected for formatting.
        schema: :py:class:`schema.SourceTableSchema`
            The relevant query schema.

        Returns
        -------
        :py:class:`astropy.table.table.Table`
            The output table.

        Notes
        -----
        By default, the only thing that this function does is correct column types.

        """
        # -- remove all table meta-data -- #
        table.meta = None
        table.schema = schema if schema is not None else cls.default_query_schema

        # -- standardize type columns -- #
        table = correct_column_types(table)

        return table

    def default_correct_query_output(
        self, table: SourceTable | Table
    ) -> SourceTable | Table:
        """
        The default query output correction callable.

        This is only used if :py:attr:`SourceTable.correct_query_output` is not specified during initialization.

        Parameters
        ----------
        table: :py:class:`astropy.table.table.Table`
            The table being corrected for formatting.

        Returns
        -------
        :py:class:`astropy.table.table.Table`
            The output table.

        Notes
        -----
        By default, the only thing that this function does is correct column types.

        """
        return self.__class__._default_correct_query_output(
            table, schema=self.query_schema
        )

    def count(
        self, positions: SkyCoord, radii: units.Quantity, parallel_kwargs: dict = None
    ):
        """
        Count the number of each object type in a set of position queries.

        Parameters
        ----------
        positions: list of SkyCoord
            A list of sky coordinates to query at.
        radii: units.Quantity
            An array of the same size as ``positions`` containing the sample radii.
        parallel_kwargs: dict
            Parameters for parallelization. For more information see :ref:`parallelization`.

            .. note::

                This method implements multi-threading in some subclasses.

        Returns
        -------
        Table
            A table containing counts for each of the object types and the positions of queries.
        """
        from pyXMIP.utilities.optimize import map_to_threads

        # -------------------------------------------#
        # Managing args / kwargs and paths
        # -------------------------------------------#
        radii = enforce_units(radii, units.arcmin)
        if radii.isscalar:
            radii = radii * np.ones(len(positions))
        # ------------------------------------------------ #
        # Running the queries through the database         #
        # ------------------------------------------------ #
        with logging_redirect_tqdm(loggers=[mainlog]), mainlog.fixed_verb(
            info=0
        ), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pbar = tqdm(zip(positions, radii), total=len(radii))
            results = map_to_threads(
                self._thread_pooled_count,
                positions,
                radii,
                repeat(pbar),
                threading_kw=parallel_kwargs,
            )

            output_table = correct_column_types(
                vstack([k for k in results if k is not None])
            )

        output_table["RA"] = [p.icrs.ra for p in positions]
        output_table["DEC"] = [p.icrs.dec for p in positions]
        output_table["RAD"] = radii
        output_table["TIME"] = time.asctime()
        return output_table

    def _thread_pooled_count(
        self, position: SkyCoord, radius: units.Quantity, progress_bar: Any
    ) -> Table:
        try:
            _output_var = self.query_radius(position, radius).count_types()
        except Exception as exception:
            mainlog.error(exception.__repr__())
            _output_var = None
        progress_bar.update()
        return _output_var

    def random_sample_count(
        self, points: int, radii: units.Quantity, parallel_kwargs: Any = None
    ) -> Table:
        """
        Count the number of instances of each object type at each of a number of random positions on the sky.

        Parameters
        ----------
        points: int
            The number of randomly sampled points to query.
        radii: units.Quantity
            The radii of the search area for each query.
        parallel_kwargs: dict
            Parameters for parallelization. For more information see :ref:`parallelization`.

            .. note::

                This method implements multi-threading in some subclasses.

        Returns
        -------
        Table
            Table of counts for each of the object types.
        """
        from pyXMIP.stats.utilities import uniform_sample_spherical

        mainlog.info(f"Querying for {points} random counts on {self.name}.")

        # -------------------------------------------#
        # Managing args / kwargs and paths
        # -------------------------------------------#
        radii = enforce_units(radii, units.arcmin)
        if radii.isscalar:
            radii = radii * np.ones(points)

        # -------------------------------------------#
        # Pull the random samples
        # -------------------------------------------#
        phi, theta = uniform_sample_spherical(points)
        theta = np.pi / 2 - theta
        positions = SkyCoord(phi, theta, frame="galactic", unit="rad")

        return self.count(positions, radii, parallel_kwargs=parallel_kwargs)

    # -- Poisson Atlas related methods -- #

    @property
    def poisson_atlas(self) -> PoissonAtlas:
        """
        The Poisson-atlas corresponding to this database.

        Returns
        -------
        :py:class:`structures.map.PoissonAtlas`
            The corresponding atlas.

        Notes
        -----
        If the user wishes to set the Poisson-atlas, they should set ``instance.poisson_atlas = path``, where the path
        directs to the file in which the atlas is stored.

        .. warning ::

            Generally, these are ``None`` by default; particularly for :py:class:`LocalDatabase` instances. For built-in
            :py:class:`RemoteDatabase` classes, there are pre-built Poisson atlases by default.
        """
        try:
            return PoissonAtlas(self._poisson_atlas_path)
        except FileNotFoundError:
            mainlog.warning(
                f"There is no Poisson atlas at {self._poisson_atlas_path}. Reverting to default."
            )
            self._poisson_atlas_path = self.__class__.default_poisson_atlas_path
            return PoissonAtlas(self._poisson_atlas_path)

    @poisson_atlas.setter
    def poisson_atlas(self, value: str):
        self._poisson_atlas_path = value

    def add_sources_to_poisson(
        self, points: int, radii: units.Quantity, parallel_kwargs: dict = None
    ):
        """
        Add randomly sampled sources to the Poisson-atlas of this database instance.

        Parameters
        ----------
        points: int
            The number of randomly sampled points to add to the :py:class:`structures.map.PoissonAtlas` instance's ``COUNTS`` table.
        radii: :py:class:`astropy.units.Quantity`
            The search radii for each of the points. Can be either an ``(n,1)`` array where ``n`` is the number of points, or a scalar value.
        parallel_kwargs: dict
            Parameters for parallelization. For more information see :ref:`parallelization`.

            .. note::

                This method implements multi-threading in some subclasses.
        """
        mainlog.info(f"Generating random sample of {points} counts.")
        point_data = self.random_sample_count(
            points, radii, parallel_kwargs=parallel_kwargs
        )

        mainlog.info(f"Adding data to the Poisson map at {self._poisson_atlas_path}.")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.poisson_atlas.append_to_fits(point_data, "COUNTS")


class LocalDatabase(SourceDatabase):
    """
    Generic representation of a local database.
    """

    def __init__(self, table, db_name, **kwargs):
        self.table = table
        super().__init__(db_name, query_schema=self.table.schema, **kwargs)

    def __str__(self):
        return f"<LocalDatabase {self.name}>"

    def __repr__(self):
        return f"<LocalDatabase {self.name}, N={len(self.table)}>"

    def source_match_memory(
        self, source_table, search_radius=1 * units.arcmin, parallel_kwargs=None
    ):
        """
        Match a :py:class:`SourceTable` against this database.

        Parameters
        ----------
        source_table: :py:class:`SourceTable`
            The table to cross match against.
        search_radius: :py:class:`astropy.units.Quantity`
            The search radii for each of the source points.
        parallel_kwargs: dict
            Kwargs for parallelism.

            .. warning::

                This kwarg has no effect for local databases.

        Returns
        -------
        :py:class:`astropy.table.table.Table`
        """
        mainlog.info(
            f"Source matching table with {len(source_table)} entries against local database {self.name}."
        )

        # -- pull the coordinates from both tables -- #
        db_positions = self.table.get_coordinates()
        other_positions = source_table.get_coordinates()

        # ---------------------------------------------------#
        # Managing args and kwargs
        # ---------------------------------------------------#

        if not isinstance(search_radius, units.Quantity):
            mainlog.warning(
                f"Search radii is a data type without standard units ({type(search_radius)}). Defaulting to arcmin."
            )
            search_radius = search_radius * units.arcmin

        if not search_radius.isscalar:
            raise ValueError("The search radius for a local database must be a scalar.")

        _ = parallel_kwargs  # --> trick IDEs
        # ---------------------------------------------------#
        # Matching
        # ---------------------------------------------------#
        idxother, idxself, separation, _ = db_positions.search_around_sky(
            other_positions, search_radius
        )

        # ---------------------------------------------------#
        # Constructing the Table
        # ---------------------------------------------------#
        # Pull the base data provided by the local database.
        matched_table = self.table[idxself]  # --> this is the matched data.

        # Merge catalog and local database information.
        matched_table["CATOBJ"] = source_table[source_table.schema.NAME][idxother]
        matched_table["CATRA"], matched_table["CATDEC"] = (
            other_positions[idxother].ra.deg,
            other_positions[idxother].dec.deg,
        )
        matched_table["CATNMATCH"] = [
            len(list(set(idxself[idxother == idxo]))) for idxo in idxother
        ]

        return matched_table

    def source_match(
        self, path, source_table, search_radius=1 * units.arcmin, parallel_kwargs=None
    ):
        """
        Match a :py:class:`SourceTable` against this database.

        Parameters
        ----------
        path: str
            The path at which to write the match data.
        source_table: :py:class:`SourceTable`
            The table to cross match against.
        search_radius: :py:class:`astropy.units.Quantity`
            The search radii for each of the source points.
        parallel_kwargs: dict
            Kwargs for parallelism.

            .. warning::

                This kwarg has no effect for local databases.

        Returns
        -------
        None
        """
        import sqlalchemy as sql

        engine = sql.create_engine(f"sqlite:///{path}")
        matched_data = self.source_match_memory(
            source_table, search_radius=search_radius, parallel_kwargs=parallel_kwargs
        )
        matched_data.append_to_sql(f"{self.name}_MATCH", engine)

    def _query_radius(self, position, radius):
        # -- Pull the matches -- #
        positions = self.table.get_coordinates()
        return self.table[positions.separation(position) < radius]


class RemoteDatabase(SourceDatabase, ABC):
    """
    Generic class representation of a source database.
    """

    def __init__(self, db_name, **kwargs):
        super().__init__(db_name, **kwargs)

    def __str__(self):
        return f"<RemoteDatabase {self.name}>"

    def source_match(
        self, path, source_table, search_radii=1 * units.arcmin, parallel_kwargs=None
    ):
        """
        Match a :py:class:`SourceTable` against this database.

        Parameters
        ----------
        path: str
            The path at which to write the match data.
        source_table: :py:class:`SourceTable`
            The table to cross match against.
        search_radii: :py:class:`astropy.units.Quantity`
            The search radii for each of the source points.
        parallel_kwargs: dict
            Kwargs for parallelism.

        Returns
        -------
        None
        """
        import sqlalchemy as sql

        from pyXMIP.utilities.optimize import map_to_threads

        mainlog.info(f"Source matching {len(source_table)} against {self.name}.")

        # ---------------------------------------------------#
        # Managing args and kwargs
        # ---------------------------------------------------#
        engine = sql.create_engine(f"sqlite:///{path}")

        if not isinstance(search_radii, units.Quantity):
            mainlog.warning(
                f"Search radii is a data type without standard units ({type(search_radii)}). Defaulting to arcmin."
            )
            search_radii = np.array(search_radii) * units.arcmin

        if search_radii.isscalar:
            search_radii = search_radii * np.ones(len(source_table))

        # ---------------------------------------------------#
        # Running queries
        # ---------------------------------------------------#
        positions = source_table.get_coordinates()

        with logging_redirect_tqdm(loggers=[mainlog]):
            pbar = tqdm(desc=f"Querying {self.name}", total=len(search_radii))
            # -- Run once without pass to threads -- #
            result = map_to_threads(
                self._thread_pooled_source_match,
                positions,
                source_table,
                repeat(source_table.schema),
                search_radii,
                repeat(pbar),
                repeat(engine),
                threading_kw=parallel_kwargs,
            )

            for _ in result:
                pass

    def _thread_pooled_source_match(
        self,
        position,
        source_table_entry,
        source_table_schema,
        search_radius,
        pbar,
        engine,
    ):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                query = self.correct_query_output(
                    self._query_radius(position, search_radius)
                )
                if not len(query):
                    pbar.update()
                    return None

                query["CATOBJ"] = source_table_entry[source_table_schema.NAME]
                query["CATRA"] = position.ra.deg
                query["CATDEC"] = position.dec.deg
                query["CATNMATCH"] = len(query)
                query.meta = {}
                pbar.update()
            except Exception as exception:
                mainlog.error(exception.__repr__())
                return None

            with self._thread_lock:
                if not sql.inspect(engine).has_table(f"{self.name}_MATCH"):
                    mainlog.info(
                        f"[{threading.current_thread().name}] Creating table {self.name}_MATCH schema."
                    )
                    metadata = sql.MetaData()
                    _ = sql.Table(
                        f"{self.name}_MATCH",
                        metadata,
                        *[
                            sql.Column(k, convert_np_type_to_sql(v))
                            for k, v in dict(query.to_pandas().dtypes).items()
                        ],
                    )
                    metadata.create_all(engine)

                query.append_to_sql(f"{self.name}_MATCH", engine, method="multi")

        return None


class NED(RemoteDatabase):
    """
    Built-in remote database access to the `IPAC / CALTECH NED <https://ned.ipac.caltech.edu/>`_ database.

    Notes
    -----

    *from NED website:*

    NED is a comprehensive database of multiwavelength data for extragalactic objects, providing a systematic, ongoing
    fusion of information integrated from hundreds of large sky surveys and tens of thousands of research publications.
    The contents and services span the entire observed spectrum from gamma rays through radio frequencies.
    As new observations are published, they are cross-identified or statistically associated with previous data and
    integrated into a unified database to simplify queries and retrieval.
    Seamless connectivity is also provided to data in NASA astrophysics mission archives (IRSA, HEASARC, MAST),
    to the astrophysics literature via ADS, and to other data centers around the world.

    """

    default_poisson_atlas_path = os.path.join(poisson_map_directory, "NED.poisson.fits")
    default_query_config = {}
    default_query_schema = DEFAULT_SOURCE_SCHEMA_REGISTRY["NED"]

    # -- NED Settings -- #
    TIMEOUT = _DatabaseConfigSetting(default=60)
    """
    float: The maximum query response time to allow without raising an error.

    Default is 60 seconds.
    """

    def __init__(self, name="NED_STD", **kwargs):
        super().__init__(name, **kwargs)

        self.config_ned()

    def config_ned(self):
        Ned.TIMEOUT = self.TIMEOUT

    @classmethod
    def _default_correct_query_output(cls, table, schema=None):
        for col in table.columns:
            table[col].format = None
            if table[col].unit == "degrees":
                table[col].unit = "deg"

        table = super()._default_correct_query_output(table, schema=schema)

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
            output = SourceTable(Ned.query_region(position, radius))
        except requests.exceptions.ConnectionError:
            raise DatabaseError(
                f"Failed to complete query [{position},{radius}] to NED due to timeout."
            )

        # -- return data if valid -- #
        output.schema = self.query_schema
        return output

    def query_object(self, object_name):
        """
        Query NED for data related to a particular object.

        Parameters
        ----------
        object_name: str
            The name of the object to pull data for.

        Returns
        -------
        :py:class:`structures.table.SourceTable`
        """
        try:
            output = SourceTable(Ned.query_object(object_name))
        except requests.exceptions.ConnectionError:
            raise DatabaseError(
                f"Failed to complete query [{object_name}] to NED due to timeout."
            )

        # -- return data if valid -- #
        output.schema = self.query_schema
        return output

    @staticmethod
    def get_image_list(object_name, item="image"):
        """
        Query NED image or spectra urls.

        Parameters
        ----------
        object_name: str
            The name of the object to pull data for.

        Returns
        -------
        list of str
            The resulting URLs available for the object.
        """
        return Ned.get_image_list(object_name, item=item)

    @staticmethod
    def get_table(object_name, table):
        """
        Query NED for a particular data table related to an object.

        Parameters
        ----------
        object_name: str
            The name of the object to pull data for.
        table: str
            The table to query for.

        Returns
        -------
        :py:class:`astropy.table.table.Table`
            The resulting output table.
        """
        return Ned.get_table(object_name, table=table)


class SIMBAD(RemoteDatabase):
    """
    Built-in remote database access to the `SIMBAD <https://simbad.cds.unistra.fr/simbad/>`_ database.

    Notes
    -----

    *from SIMBAD website:*

    The SIMBAD database is managed by the Centre de Données astronomiques de Strasbourg (CDS).

    SIMBAD is the acronym for: **S** et of **I** dentifications, **M** easurments and **B** ibliography for
    **A** stronomical **D** ata.

    The SIMBAD software is developed by Marc Wenger (retired), Anaïs Oberto, Grégory Mantelet (CDS, Strasbourg).
    with contributions of students during trainings.

    The updating of SIMBAD is a continuous process.

    """

    default_poisson_atlas_path = os.path.join(
        poisson_map_directory, "SIMBAD.poisson.fits"
    )
    default_query_schema = DEFAULT_SOURCE_SCHEMA_REGISTRY["SIMBAD"]
    default_query_config = {}
    TIMEOUT = _DatabaseConfigSetting(default=120)
    """
    float: The maximum query response time to allow without raising an error.

    Default is 120 seconds.
    """
    EXTRA_COLUMNS = _DatabaseConfigSetting(default=["otypes", "ra(d;A)", "dec(d;D)"])
    """
    list: Additional VOTable columns to add to the returned query tables.

    Default is ``['otypes','ra(d;A)','dec(d;D)']``.
    """
    REMOVED_COLUMNS = _DatabaseConfigSetting(default=[])
    """
    list: VOTable columns to remove from the default list.

    Default is ``[]``.
    """
    ROW_LIMIT = _DatabaseConfigSetting(default=0)
    """
    int: The maximum allowed rows to return in a query.

    Default is ``0``, which allows for all results.
    """

    def __init__(self, name="SIMBAD_STD", **kwargs):
        super().__init__(name, **kwargs)

        self.config_simbad()

    @property
    def available_tables(self):
        """
        list: The tables accessible via SIMBAD.
        """
        return Simbad.list_tables()

    @staticmethod
    def clear_cache():
        """Clear the cache associated with SIMBAD queries."""
        Simbad.clear_cache()

    @staticmethod
    def list_columns(*args, **kwargs):
        """
        List the available columns.
        """
        return Simbad.list_columns(*args, **kwargs)

    @staticmethod
    def list_linked_tables(*args, **kwargs):
        """List the linked tables."""
        return Simbad.list_linked_tables(*args, **kwargs)

    def config_simbad(self):
        Simbad.TIMEOUT = self.TIMEOUT
        Simbad.ROW_LIMIT = self.ROW_LIMIT
        Simbad.add_votable_fields(*self.EXTRA_COLUMNS)
        Simbad.remove_votable_fields(*self.REMOVED_COLUMNS)

    @classmethod
    def _default_correct_query_output(cls, table: SourceTable, schema=None):
        if "RA" in table.columns:
            table["RA"] = Angle(table["RA"], unit="h").deg

        if "DEC" in table.columns:
            table["DEC"] = Angle(table["DEC"], unit="deg").deg

        table = super()._default_correct_query_output(table, schema=schema)

        return table

    def _query_radius(self, position: SkyCoord, radius: units.Quantity) -> SourceTable:
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
            output = SourceTable(Simbad.query_region(position, radius))
        except requests.exceptions.ConnectionError:
            raise DatabaseError(
                f"Failed to complete query [{position},{radius}] to Simbad due to timeout."
            )

        # -- return data if valid -- #
        output.schema = self.query_schema
        return output

    def query_object(self, object_name: str):
        """
        Query SIMBAD for data related to a particular object.

        Parameters
        ----------
        object_name: str
            The name of the object to pull data for.

        Returns
        -------
        :py:class:`structures.table.SourceTable`
        """
        try:
            output = SourceTable(Simbad.query_object(object_name))
        except requests.exceptions.ConnectionError:
            raise DatabaseError(
                f"Failed to complete query [{object_name}] to Simbad due to timeout."
            )

        # -- return data if valid -- #
        output.schema = self.query_schema
        return output

    @staticmethod
    def query_tap(query: str):
        """
        Query Simbad via ASQL query.

        Parameters
        ----------
        query: str
            The ASQL query to run.

        Returns
        -------
        :py:class:`astropy.table.table.Table`
            The resulting output table.
        """
        return Simbad.query_tap(query)


class DBRegistry(Registry):
    """
    A :py:class:`DBRegistry` instance is a collection of identifiable database classes.

    Notes
    -----

    The primary purpose of database registries is to allow the user and the backend to keep a "list" of the available
    databases.
    """

    def __init__(self, databases: list[SourceDatabase]):
        """
        Initialize a database registry.

        Parameters
        ----------
        databases: list of SourceDatabase
            The list of databases to load into the registry.

        """
        super().__init__({database.name: database for database in databases})

        # -- check that the classes are actually valid -- #
        assert all(isinstance(v, SourceDatabase) for v in self.values())

    @property
    def names(self) -> list[str]:
        return list(self.keys())

    @property
    def locals(self) -> list[str]:
        return [k for k, v in self.items() if isinstance(v, LocalDatabase)]

    @property
    def remotes(self) -> list[str]:
        return [k for k, v in self.items() if isinstance(v, RemoteDatabase)]

    @classmethod
    def _default_registry(cls):
        return cls([NED(), SIMBAD()])


DEFAULT_DATABASE_REGISTRY = DBRegistry._default_registry()

# ================================================================ #
# Additional Utility Functions                                     #
# ================================================================ #


def builtin_databases():
    """List the available built-in databases."""
    return DEFAULT_DATABASE_REGISTRY.names


def get_database(database_name, registry=None):
    """
    Fetch the initialized database from a given registry.
    Parameters
    ----------
    database_name: str
        The database name to fetch.
    registry: DBRegistry, optional
        The registry to fetch from. If not specified, will search through the default registry.

    Returns
    -------
    SourceDatabase
        The output database.
    """
    if not registry:
        return DEFAULT_DATABASE_REGISTRY.get(database_name)
    else:
        return registry.get(database_name)


def get_poisson_path(database_name, registry=None):
    """
    Fetch the path to the Poisson Atlas for a given database.
    Parameters
    ----------
    database_name: str
        The name of the database to fetch.
    registry: DBRegistry, optional
        The registry to search.

    Returns
    -------
    str
    """
    return get_database(database_name, registry=registry)._poisson_atlas_path


def get_poisson_atlas(database_name, registry=None):
    """
    Fetch the path to the Poisson Atlas for a given database.
    Parameters
    ----------
    database_name: str
        The name of the database to fetch.
    registry: DBRegistry, optional
        The registry to search.

    Returns
    -------
    str
    """
    return get_database(database_name, registry=registry).poisson_atlas


def add_points_to_poisson_map(database_name, n, radii, registry=None, **kwargs):
    """
    Add new points to an existing :py:class:`structures.map.PoissonAtlas` attached to a database.

    Parameters
    ----------
    database_name: str
        The database to load.
    n: int
        The number of new points to add.
    radii: Quantity
        Value (scalar or vector) for the radius of search for each sample.
    registry: DBRegistry, optional
        The registry to search. If not provided, searches the default registry.
    **kwargs:
        Additional kwargs to pass to the random counting process.

    Returns
    -------
    None
    """
    if isinstance(radii, str):
        radii = float(radii.split(" ")[0]) * units.Unit(radii.split(" ")[1])
    elif isinstance(radii, (float, int)):
        radii = radii * units.arcmin

    database = get_database(database_name, registry)

    database.add_sources_to_poisson(n, radii, **kwargs)
