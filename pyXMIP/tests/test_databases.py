import pytest
from astropy import units
from astropy.coordinates import SkyCoord

from pyXMIP.structures.databases import NED, SIMBAD
from pyXMIP.tests.utils import table_answer_testing


@pytest.mark.usefixtures("answer_dir", "answer_store")
class TestDatabase:
    query_position = SkyCoord(ra=1, dec=1, unit="deg")
    database = None

    def test_query(self, answer_dir, answer_store):
        if self.__class__.database is not None:
            data = self.__class__.database.query_radius(
                self.__class__.query_position, 1 * units.arcmin
            )
        else:
            return None

        table_answer_testing(
            data,
            f"{self.__class__.database.__class__.__name__}_query_test.fits",
            answer_store,
            answer_dir,
        )


@pytest.mark.usefixtures("answer_dir", "answer_store")
class TestNED(TestDatabase):
    database = NED


@pytest.mark.usefixtures("answer_dir", "answer_store")
class TestSIMBAD(TestDatabase):
    database = SIMBAD
