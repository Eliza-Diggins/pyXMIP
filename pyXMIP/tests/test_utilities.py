import numpy as np
from numpy.testing import assert_allclose

from pyXMIP.utilities.geo import convert_coordinates


def test_coord_conversion():
    from itertools import product

    coordinate_test_set = {
        "healpix": [(0, 0), (0, np.pi / 2), (0, np.pi), (np.pi, 0)],
        "latlon": [(0, np.pi / 2), (0, 0), (0, -np.pi / 2), (np.pi, np.pi / 2)],
    }

    for coordsystem1, coordsystem2 in list(
        product(list(coordinate_test_set.keys()), list(coordinate_test_set.keys()))
    ):
        for i in range(4):
            assert_allclose(
                np.array(coordinate_test_set[coordsystem2][i]),
                convert_coordinates(
                    *coordinate_test_set[coordsystem1][i],
                    from_system=coordsystem1,
                    to_system=coordsystem2,
                ),
                err_msg=f"failed for {coordsystem2},{coordsystem1}",
            )
