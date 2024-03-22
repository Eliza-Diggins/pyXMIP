"""
Database module for X-matching and querying databases.
"""
import threading
import time
from abc import ABC
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
import requests.exceptions
from astroquery.ipac.ned import Ned
from pyXMIP.structures.table import SourceTable
from astropy.table import vstack
from pyXMIP.schema import SourceTableSchema
from pyXMIP.utilities.core import _bin_directory, mainlog, enforce_units
import os
import astropy.coordinates as coords
from astropy import units
import pathlib as pt
import numpy as np
from astropy.io import fits
from itertools import repeat
from time import asctime
poisson_map_directory = os.path.join(_bin_directory,'psn_maps')


class DatabaseError(Exception):
    """
    Collective error type for issues with database queries.
    """

    def __init__(self,message=None):
        self.message = message
        super().__init__(self.message)

class SourceDatabase(ABC):
    """
    Generic class representation of a source database.
    """

class RemoteDatabase(SourceDatabase):
    """
    Class representation of general remotely hosted databases which can be queried for source objects.
    """

    class_table_schema = None

    @classmethod
    def query_radius(cls,position,radius):
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
        raise NotImplementedError

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

        #-------------------------------------------#
        # Managing args / kwargs and paths
        #-------------------------------------------#
        radii = enforce_units(radii,units.arcmin)
        if radii.isscalar:
            radii = radii * np.ones(points)

        #-------------------------------------------#
        # Pull the random samples
        #-------------------------------------------#
        phi,theta = uniform_sample_spherical(points)
        theta = np.pi/2 - theta
        positions = SkyCoord(phi, theta, frame='galactic', unit='rad')
        #-------------------------------------------#
        # Managing Threads
        #-------------------------------------------#
        if thread_kw is None:
            thread_kw = {}

        if threading:
            from pyXMIP.utilities.mp_utils import split
            from concurrent.futures import ThreadPoolExecutor
            mainlog.debug("Querying with threading.")
            max_workers = thread_kw.pop('max_workers', 5)
            chunk_size = thread_kw.pop('chunk_size', 5)
            c_radii = split(radii,len(radii)//chunk_size)
            c_positions = split(positions,len(positions)//chunk_size)

            with logging_redirect_tqdm(loggers=[mainlog]):
                pbar = tqdm("Querying",total=len(radii))
                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    results = executor.map(cls._count_multithreaded,c_positions,c_radii,repeat(pbar))

            output_tables = [r for r in results]
            output_tables = vstack(output_tables)
            output_tables['RA'] = [p.icrs.ra for p in positions]
            output_tables['DEC'] = [p.icrs.dec for p in positions]
            output_tables['RAD'] = radii
            output_tables['TIME'] = time.asctime()
            return output_tables

        else:
            mainlog.debug("Querying without threading.")
            return cls.count(positions,radii)

    @classmethod
    def _count_multithreaded(cls,positions,radii,progress_bar):
        output_tables = []

        for position,radius in zip(positions,radii):
            # Iterate over all of the queries and grab them all.
            try:
                output_tables.append(cls.query_radius(position,radius).count_types())
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
            for position, radius in tqdm(zip(positions, radii),total=len(radii)):
                # Iterate over all of the queries and grab them all.
                try:
                    output_tables.append(cls.query_radius(position, radius).count_types())
                except Exception as exception:
                    mainlog.error(exception.__repr__())

        output_table = vstack(output_tables)
        output_table['RA'] = [p.icrs.ra for p in positions]
        output_table['DEC'] = [p.icrs.dec for p in positions]
        output_table['RAD'] = radii
        output_table['TIME'] = time.asctime()
        return output_table

    @classmethod
    def source_match(cls,path,source_table,search_radii=1*units.arcmin,threading=True,thread_kw=None):
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
        mainlog.info(f"Source matching {len(source_table)} against {cls.__name__}.")

        #---------------------------------------------------#
        # Managing args and kwargs
        #---------------------------------------------------#
        path = pt.Path(path)
        if not path.exists():
            hdr = fits.Header()
            hdr['DATE'] = asctime()
            fits.HDUList([fits.PrimaryHDU(header=hdr)]).writeto(path)

        if not isinstance(search_radii,units.Quantity):
            mainlog.warning(f"Search radii is a data type without standard units ({type(search_radii)}). Defaulting to arcmin.")
            search_radii = np.array(search_radii) * units.arcmin

        if search_radii.isscalar:
            search_radii = search_radii * np.ones(len(source_table))

        if thread_kw is None: thread_kw = {} # --> Keyword arguments for threading.

        #---------------------------------------------------#
        # Running queries
        #---------------------------------------------------#
        positions = source_table.get_coordinates()
        if threading:
            from pyXMIP.utilities.mp_utils import split
            from concurrent.futures import ThreadPoolExecutor
            mainlog.debug("Querying with threading.")
            max_workers = thread_kw.pop('max_workers', 5)
            chunk_size = thread_kw.pop('chunk_size', 5)
            c_radii = split(search_radii,len(search_radii)//chunk_size)
            c_positions = split(positions,len(positions)//chunk_size)
            c_entries = split(source_table,len(source_table)//chunk_size)

            with logging_redirect_tqdm(loggers=[mainlog]):
                pbar = tqdm("Querying",total=len(search_radii))
                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    executor.map(cls._source_match_multithreaded,c_positions,c_entries,c_radii,repeat(pbar),repeat(path))

        else:
            mainlog.info(f"Querying without threading. [Calls={len(source_table)}]")
            with logging_redirect_tqdm(loggers=[mainlog]):
                for _position, _st_entry, _radius in tqdm(zip(positions,source_table,search_radii),total=len(source_table)):
                    try:
                        query = cls.query_radius(_position,
                                                 _radius)
                    except Exception as exception:
                        mainlog.error(exception.__repr__())

                    query['CATALOG_OBJECT'] = _st_entry[source_table.schema.NAME]
                    query['CATALOG_RA'] = _position.ra
                    query['CATALOG_DEC'] = _position.dec


                    # -- Writing the query to the path -- #
                    query = query[
                        ['CATALOG_OBJECT',
                         'CATALOG_RA',
                         'CATALOG_DEC',
                         cls.class_table_schema.coordinate_columns[0],
                         cls.class_table_schema.coordinate_columns[1],
                         cls.class_table_schema.NAME,
                         cls.class_table_schema.TYPE]
                    ]
                    query.append_to_fits(path,f'{cls.__name__}_MATCH')

    @classmethod
    def _source_match_multithreaded(cls,positions,source_table,search_radii,pbar,path):
        output_tables = []
        for _position, _st_entry, _radius in zip(positions, source_table, search_radii):
            try:
                query = cls.query_radius(_position,
                                         _radius)
            except Exception as exception:
                mainlog.error(exception.__repr__())

            query['CATALOG_OBJECT'] = _st_entry[source_table.schema.NAME]
            query['CATALOG_RA'] = _position.ra
            query['CATALOG_DEC'] = _position.dec
            query.meta = {}
            output_tables.append(query)

            pbar.update()

        # -- Writing the query to the path -- #
        output_tables = vstack(output_tables)
        output_tables = output_tables[
            ['CATALOG_OBJECT',
             'CATALOG_RA',
             'CATALOG_DEC',
             cls.class_table_schema.coordinate_columns[0],
             cls.class_table_schema.coordinate_columns[1],
             cls.class_table_schema.NAME,
             cls.class_table_schema.TYPE]
        ]
        with threading.Lock():
            output_tables.append_to_fits(path, f'{cls.__name__}_MATCH')


class LocalDatabase(SourceDatabase):
    """
    Generic representation of a local database.
    """
    def __init__(self,filepath,schema=None,tablename=None):
        pass

class NED(RemoteDatabase):
    """
    Class instantiation of the IPAC / CALTECH NED database.
    """
    class_table_schema = SourceTableSchema.from_file(os.path.join(_bin_directory, 'builtin_schema', 'NED.yaml'))

    @classmethod
    def query_radius(cls,position,radius):
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
            output = SourceTable(Ned.query_region(position,radius))
        except requests.exceptions.ConnectionError as error:
            raise DatabaseError(f"Failed to complete query [{position},{radius}] to NED due to timeout.")

        # -- return data if valid -- #
        output.schema = cls.class_table_schema
        return output



if __name__ == '__main__':
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import SkyCoord
    import matplotlib.pyplot as plt
    #fits_table = SourceTable.read('/home/ediggins/pyROSITA_test/eRASS1_Hard.v1.0.fits',format='fits')

    #NED.source_match('test.fits',fits_table,threading=True)

    with fits.open('test.fits','update') as hudl:
        print(hudl.info())