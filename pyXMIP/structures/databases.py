"""
pyXMIP database implementations for easy querying of databases.
"""
from abc import ABC, abstractmethod
from astroquery.ipac.ned import Ned
from table import SourceTable
from pyXMIP.schema import SkyCollectionSchema
from pyXMIP.utilities.core import _bin_directory
import os

class SourceDatabase(ABC):
    """
    Generic class representation of a source database.
    """

class RemoteDatabase(SourceDatabase):
    """
    Generic representation of remote hosted source databases.
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
        :py:class:`astropy.table.Table`
        """
        raise NotImplementedError

class LocalDatabase(SourceDatabase):
    """
    Generic representation of a local database.
    """
    def __init__(self,filepath,schema=None,tablename=None):
        pass


class NED(RemoteDatabase):
    class_table_schema = SkyCollectionSchema.from_file(os.path.join(_bin_directory,'builtin_schema','NED.yaml'))

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
        output = SourceTable(Ned.query_region(position,radius))
        output.schema = cls.class_table_schema
        return output

if __name__ == '__main__':
    from astropy.coordinates import SkyCoord

    u = SkyCoord(ra=1,dec=1,unit='deg')
    from astropy.units import arcmin
    print(NED.query_radius(u,5*arcmin))