"""
Database module for X-matching and querying databases.
"""
import os
import threading
import time
import warnings
from abc import ABC
from itertools import repeat

import numpy as np
import requests.exceptions
from astropy import units
from astropy.coordinates import Angle, SkyCoord
from astropy.table import vstack
from astroquery.ipac.ned import Ned
from astroquery.simbad import Simbad
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from pyXMIP.schema import SourceTableSchema
from pyXMIP.structures.map import PoissonAtlas
from pyXMIP.structures.table import SourceTable
from pyXMIP.utilities._registries import _Registry
from pyXMIP.utilities.core import _bin_directory, enforce_units, mainlog

poisson_map_directory = os.path.join(_bin_directory, "psn_maps")


class DatabaseError(Exception):
    """
    Collective error type for issues with database queries.
    """

    def __init__(self, message=None):
        self.message = message
        super().__init__(self.message)


class SourceDatabase(ABC):
    """
    Generic class representation of a source database.
    """

    class_table_schema = None
    default_poisson_map = os.path.join(poisson_map_directory, "NONE")
    poisson_map = default_poisson_map

    @classmethod
    def _fix_output_table_formatting(cls, table):
        """This is a backend function for fixing query tables with non-standard outputs."""

        # -- remove all table meta-data -- #
        table.meta = None
        table.schema = cls.class_table_schema
        # -- Fix object column types -- #
        for col in table.columns:
            if table[col].dtype == "object":
                table[col] = np.array([str(j) for j in table[col]], dtype="<U8")

        if table.schema.TYPE in table.columns:
            # fix the types
            table[table.schema.TYPE] = [
                f"|{k}" if k[0] != "|" else k for k in table[table.schema.TYPE]
            ]
            table[table.schema.TYPE] = [
                f"{k}|" if k[-1] != "|" else k for k in table[table.schema.TYPE]
            ]

        return table

    @classmethod
    def _query_radius_uncorrected(cls, position, radius):
        raise NotImplementedError

    @classmethod
    def query_radius(cls, position, radius):
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
        """
        output_table = cls._query_radius_uncorrected(position, radius)
        if len(output_table):
            return cls._fix_output_table_formatting(output_table)
        else:
            return output_table

    @classmethod
    def random_sample_count(cls, points, radii, threading=True, thread_kw=None):
        """
        Count the number of instances of each object type at each of a number of random positions on the sky.

        Parameters
        ----------
        points: int
            The number of randomly sampled points to query.
        radii: units.Quantity
            The radii of the search area for each query.
        threading: bool, optional
            If ``True``, the queries will be threaded.
        thread_kw: dict, optional
            Additional kwargs to dictate the threading behavior.

        Returns
        -------
        Table
            Table of counts for each of the object types.
        """
        from pyXMIP.stats.utilities import uniform_sample_spherical

        mainlog.info(f"Querying for {points} random counts on {cls.__name__}.")

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
        # -------------------------------------------#
        # Managing Threads
        # -------------------------------------------#
        if thread_kw is None:
            thread_kw = {}

        if threading:
            from concurrent.futures import ThreadPoolExecutor

            from pyXMIP.utilities._mp_utils import split

            mainlog.debug("Querying with threading.")
            max_workers = thread_kw.pop("max_workers", 5)
            chunk_size = thread_kw.pop("chunk_size", 5)
            c_radii = split(radii, len(radii) // chunk_size)
            c_positions = split(positions, len(positions) // chunk_size)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                with logging_redirect_tqdm(loggers=[mainlog]):
                    pbar = tqdm("Querying", total=len(radii))
                    with ThreadPoolExecutor(max_workers=max_workers) as executor:
                        results = executor.map(
                            cls._count_multithreaded, c_positions, c_radii, repeat(pbar)
                        )

            output_tables = [r for r in results]
            output_tables = vstack(output_tables)
            output_tables = cls._fix_output_table_formatting(output_tables)
            output_tables["RA"] = [p.icrs.ra for p in positions]
            output_tables["DEC"] = [p.icrs.dec for p in positions]
            output_tables["RAD"] = radii
            output_tables["TIME"] = time.asctime()
            return output_tables

        else:
            mainlog.debug("Querying without threading.")
            return cls.count(positions, radii)

    @classmethod
    def _count_multithreaded(cls, positions, radii, progress_bar):
        output_tables = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for position, radius in zip(positions, radii):
                # Iterate over all of the queries and grab them all.
                try:
                    output_tables.append(
                        cls._query_radius_uncorrected(position, radius).count_types()
                    )
                except Exception as exception:
                    mainlog.error(exception.__repr__())

                progress_bar.update()

            return vstack(output_tables)

    @classmethod
    def count(cls, positions, radii):
        """
        Count the number of each object type in a set of position queries.

        Parameters
        ----------
        positions: list of SkyCoord
            A list of sky coordinates to query at.
        radii: units.Quantity
            An array of the same size as ``positions`` containing the sample radii.

        Returns
        -------
        Table
            A table containing counts for each of the object types and the positions of queries.
        """
        output_tables = []
        with logging_redirect_tqdm(loggers=[mainlog]):
            for position, radius in tqdm(zip(positions, radii), total=len(radii)):
                # Iterate over all of the queries and grab them all.
                try:
                    output_tables.append(
                        cls._query_radius_uncorrected(position, radius).count_types()
                    )
                except Exception as exception:
                    mainlog.error(exception.__repr__())

        output_table = vstack(output_tables)
        output_table = cls._fix_output_table_formatting(output_table)
        output_table["RA"] = [p.icrs.ra for p in positions]
        output_table["DEC"] = [p.icrs.dec for p in positions]
        output_table["RAD"] = radii
        output_table["TIME"] = time.asctime()
        return output_table

    @classmethod
    def source_match(
        cls,
        path,
        source_table,
        search_radii=1 * units.arcmin,
        threading=True,
        thread_kw=None,
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
        threading: bool, optional
            If ``True``, multithreading will be used.
        thread_kw: dict, optional
            Additional dictionary with threading keywords.

        Returns
        -------
        None
        """
        import sqlalchemy as sql

        mainlog.info(f"Source matching {len(source_table)} against {cls.__name__}.")

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

        if thread_kw is None:
            thread_kw = {}  # --> Keyword arguments for threading.

        # ---------------------------------------------------#
        # Running queries
        # ---------------------------------------------------#
        positions = source_table.get_coordinates()
        if threading:
            from concurrent.futures import ThreadPoolExecutor

            from pyXMIP.utilities._mp_utils import split

            mainlog.debug("Querying with threading.")
            max_workers = thread_kw.pop("max_workers", 5)
            chunk_size = thread_kw.pop("chunk_size", 5)
            c_radii = split(search_radii, len(search_radii) // chunk_size)
            c_positions = split(positions, len(positions) // chunk_size)
            c_entries = split(source_table, len(source_table) // chunk_size)

            with logging_redirect_tqdm(loggers=[mainlog]):
                pbar = tqdm(
                    desc=f"Querying {cls.__name__}",
                    total=len(search_radii),
                    leave=False,
                )
                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    executor.map(
                        cls._source_match_multithreaded,
                        c_positions,
                        c_entries,
                        c_radii,
                        repeat(pbar),
                        repeat(engine),
                    )

        else:
            mainlog.info(f"Querying without threading. [Calls={len(source_table)}]")
            with logging_redirect_tqdm(loggers=[mainlog]):
                for _position, _st_entry, _radius in tqdm(
                    zip(positions, source_table, search_radii),
                    desc=f"Querying {cls.__name__}",
                    total=len(source_table),
                    leave=False,
                ):
                    try:
                        query = cls._query_radius_uncorrected(_position, _radius)
                    except Exception as exception:
                        mainlog.error(exception.__repr__())

                    query["CATALOG_OBJECT"] = _st_entry[source_table.schema.NAME]
                    query["CATALOG_RA"] = _position.ra
                    query["CATALOG_DEC"] = _position.dec

                    # -- Writing the query to the path -- #
                    query = query[
                        [
                            "CATALOG_OBJECT",
                            "CATALOG_RA",
                            "CATALOG_DEC",
                            cls.class_table_schema.coordinate_columns[0],
                            cls.class_table_schema.coordinate_columns[1],
                            cls.class_table_schema.NAME,
                            cls.class_table_schema.TYPE,
                        ]
                    ]
                    query = cls._fix_output_table_formatting(query)
                    query.append_to_sql(f"{cls.__name__}_MATCH", engine)

    @classmethod
    def _source_match_multithreaded(
        cls, positions, source_table, search_radii, pbar, engine
    ):
        output_tables = []
        for _position, _st_entry, _radius in zip(positions, source_table, search_radii):
            try:
                query = cls._query_radius_uncorrected(_position, _radius)
            except Exception as exception:
                mainlog.error(exception.__repr__())
                continue

            query["CATALOG_OBJECT"] = _st_entry[source_table.schema.NAME]
            query["CATALOG_RA"] = _position.ra
            query["CATALOG_DEC"] = _position.dec
            query.meta = {}
            output_tables.append(query)

            pbar.update()

        # -- Writing the query to the path -- #
        output_tables = vstack(output_tables)
        output_tables = output_tables[
            [
                "CATALOG_OBJECT",
                "CATALOG_RA",
                "CATALOG_DEC",
                cls.class_table_schema.coordinate_columns[0],
                cls.class_table_schema.coordinate_columns[1],
                cls.class_table_schema.NAME,
                cls.class_table_schema.TYPE,
            ]
        ]
        output_tables = cls._fix_output_table_formatting(output_tables)
        with threading.Lock():
            output_tables.append_to_sql(f"{cls.__name__}_MATCH", engine)

    @classmethod
    def get_poisson_atlas(cls):
        """Get the currently set PoissonMap instance."""
        try:
            return PoissonAtlas(cls.poisson_map)
        except FileNotFoundError:
            mainlog.warning(
                f"Class {cls.__name__} has no map at {cls.poisson_map}. Defaulting."
            )
            return cls.get_default_poisson_atlas()

    @classmethod
    def get_default_poisson_atlas(cls):
        try:
            return PoissonAtlas(cls.default_poisson_map)
        except FileNotFoundError:
            raise ValueError(
                f"It appears there is no existing PoissonAtlas at {cls.default_poisson_map}. You will have to generate one."
            )

    @classmethod
    def set_poisson_atlas(cls, path):
        cls.poisson_map = path

    @classmethod
    def add_sources_to_poisson(cls, points, radii, threading=True, thread_kw=None):
        mainlog.info(f"Generating random sample of {points} counts.")
        point_data = cls.random_sample_count(
            points, radii, threading=threading, thread_kw=thread_kw
        )

        mainlog.info(f"Adding data to the Poisson map at {cls.poisson_map}.")
        psmap = cls.get_poisson_atlas()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            psmap.append_to_fits(point_data, "COUNTS")

    @classmethod
    def create_poisson_atlas(cls, *args, **kwargs):
        kwargs["database"] = cls.__name__
        path = kwargs.pop("path", cls.default_poisson_map)
        psn = PoissonAtlas.generate(path, *args, **kwargs)

        cls.set_poisson_atlas(psn.path)
        return psn


class RemoteDatabase(SourceDatabase):
    """
    Class representation of general remotely hosted databases which can be queried for source objects.
    """


class LocalDatabase(SourceDatabase):
    """
    Generic representation of a local database.
    """


class NED(RemoteDatabase):
    """
    Class instantiation of the IPAC / CALTECH NED database.
    """

    default_poisson_map = os.path.join(poisson_map_directory, "NED.poisson.fits")
    poisson_map = default_poisson_map
    class_table_schema = SourceTableSchema.from_file(
        os.path.join(_bin_directory, "builtin_schema", "NED.yaml")
    )

    @classmethod
    def _fix_output_table_formatting(cls, table):
        for col in table.columns:
            table[col].format = None
            if table[col].unit == "degrees":
                table[col].unit = "deg"

        table = super()._fix_output_table_formatting(table)

        return table

    @classmethod
    def _query_radius_uncorrected(cls, position, radius):
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
        output.schema = cls.class_table_schema
        return output


class SIMBAD(RemoteDatabase):
    """
    Class instantiation of the SIMBAD database.
    """

    Simbad.add_votable_fields("otypes")
    Simbad.TIMEOUT = 120
    default_poisson_map = os.path.join(poisson_map_directory, "SIMBAD.poisson.fits")
    poisson_map = default_poisson_map
    class_table_schema = SourceTableSchema.from_file(
        os.path.join(_bin_directory, "builtin_schema", "SIMBAD.yaml")
    )

    @classmethod
    def _fix_output_table_formatting(cls, table):
        if "RA" in table.columns:
            table["RA"] = Angle(table["RA"], unit="h").deg

        if "DEC" in table.columns:
            table["DEC"] = Angle(table["DEC"], unit="deg").deg

        table = super()._fix_output_table_formatting(table)

        return table

    @classmethod
    def _query_radius_uncorrected(cls, position, radius):
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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                output = SourceTable(Simbad.query_region(position, radius))
            except requests.exceptions.ConnectionError:
                raise DatabaseError(
                    f"Failed to complete query [{position},{radius}] to SIMBAD due to timeout."
                )

        # -- return data if valid -- #
        output.schema = cls.class_table_schema
        return output


class DBRegistry(_Registry):
    def __init__(self, mapping):
        super().__init__(mapping)

        assert all(issubclass(k, SourceDatabase) for k in list(mapping.values()))

    @classmethod
    def _default_registry(cls):
        _mapping = {"NED": NED, "SIMBAD": SIMBAD}
        return cls(_mapping)


DEFAULT_DATABASE_REGISTRY = DBRegistry._default_registry()


def add_points_to_poisson_map(database_name, n, radii, **kwargs):
    if isinstance(radii, str):
        radii = float(radii.split(" ")[0]) * units.Unit(radii.split(" ")[1])
    elif isinstance(radii, (float, int)):
        radii = radii * units.arcmin

    database = DEFAULT_DATABASE_REGISTRY[database_name]

    database.add_sources_to_poisson(n, radii, **kwargs)


if __name__ == "__main__":
    from astropy import units

    p = SIMBAD.query_radius(SkyCoord(ra=1, dec=1, unit="deg"), 10 * units.arcmin)
    print(list(p["OTYPES"]))
    print(p.count_types()["LM*"])
    # print(p['TYPE'])
    # print(p.count_types())
