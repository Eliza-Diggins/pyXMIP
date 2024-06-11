"""
Testing suite for the :py:mod:`pyXMIP.databases` module.
"""
import os

import pytest
from astropy import units
from astropy.coordinates import SkyCoord

from pyXMIP.structures.databases import DEFAULT_DATABASE_REGISTRY, RemoteDatabase
from pyXMIP.tests.utils import check_astropy_table

database_answer_subdir = "database_answers"
query_radius_check_position = SkyCoord(ra=1, dec=3, unit="deg")
radius = 2 * units.arcmin


@pytest.mark.parametrize(
    "database",
    [DEFAULT_DATABASE_REGISTRY[j] for j in DEFAULT_DATABASE_REGISTRY.remotes],
)
def test_radius_query(database: RemoteDatabase, answer_store, answer_dir):
    """
    Queries each database for a specified radius and position. Checks against saved answer tables.
    """
    query = database.query_radius(query_radius_check_position, radius=radius)
    check_astropy_table(
        query,
        answer_dir,
        os.path.join(database_answer_subdir, f"{database.name}_query_table.fits"),
        answer_store,
    )


class RemoteDatabaseTest:
    """
    General class for setting up and accessing testing resources for the database tests.
    """

    remote_database_name = None  # The name of the remote database.
