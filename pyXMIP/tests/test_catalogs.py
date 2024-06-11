"""
Testing module for pyXMIP's catalogs.
"""
import os
import pathlib as pt

import pytest

from pyXMIP import load
from pyXMIP.utilities.core import bin_directory

test_catalog_path = pt.Path(os.path.join(bin_directory, "testobj", "test_catalog.fits"))


class TestCatalog:
    """
    Test the base catalog.
    """

    @pytest.fixture(scope="class")
    def base_catalog(self):
        yield load(test_catalog_path)

    def test_catalog_exists(self, base_catalog):
        """
        Test that the sample catalog does exist.
        """
        assert 1, "The sample catalog loaded successfully."

    def test_schema(self, base_catalog):
        """
        Test that a schema exists and has the information we anticipate.
        """
        _expected_column_map = {
            "RA": {"name": "RA", "unit": "deg"},
            "DEC": {"name": "DEC", "unit": "deg"},
            "RA_ERR": None,
            "DEC_ERR": None,
            "GAL_L": {"name": "LII", "unit": "deg"},
            "GAL_B": {"name": "BII", "unit": "deg"},
            "GAL_L_ERR": None,
            "GAL_B_ERR": None,
            "SEP": None,
            "NAME": {"name": "IAUNAME", "unit": None},
            "Z": None,
            "TYPE": None,
        }

        assert (
            _expected_column_map == base_catalog.schema.column_map.model_dump()
        ), "The test catalog has a different column map in its schema than expected."

    def test_coordinates(self, base_catalog):
        """
        Simply check that RA and DEC are the same as we expect them to be.
        """
        from numpy.testing import assert_allclose

        lat, lon = base_catalog.lat.to_value("deg"), base_catalog.lon.to_value("deg")
        _elat = [
            -31.0361403,
            -28.53148995,
            -29.51964773,
            -31.20848839,
            -21.6990906,
            -23.18313657,
            -19.9442472,
            -18.54302228,
            -16.23896735,
            -16.81838409,
        ]
        _elon = [
            7.40528954,
            10.52964788,
            12.34260232,
            13.7041894,
            22.29515815,
            20.90908527,
            24.37178961,
            28.07930583,
            29.33960449,
            30.15507953,
        ]

        assert_allclose(_elat, lat, rtol=1e-3)
        assert_allclose(_elon, lon, rtol=1e-3)
