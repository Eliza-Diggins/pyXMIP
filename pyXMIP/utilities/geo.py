r"""Module for managing conversion between different coordinate conventions.

Notes
-----
The utilities in this module are largely focused on managing the interface between HEALPix, astronomical coordinates,
and mathematical convention.

In pyXMIP, coordinates are specified in 1 of 3 systems:

- **Mathematical Convention**: In the mathematical convention, spherical coordinates are denoted :math:`(r,\phi,\theta)`,
  where :math:`\phi` is the azimuthal angle measured from the :math:`x`-axis counter-clockwise around positive :math:`z`.
  :math:`\theta` is the elevation angle, measured from the :math:`xy`-plane to :math:`\pm \pi/2`.
- **HEALPix Convention**: The HEALPix coordinate convention is marginally different from the others in that it measures the
  elevation angle from :math:`0\to \pi`. In this convention, the angle :math:`\omega` is the *co-latitude* and measures from the
  north pole. The longitude :math:`\phi` is measured from :math:`0\to 2\pi` going eastward.
- **Coordinate Convention**: This is the convention for all the astronomical coordinate systems. These obey longitude and
  latitude as it is typically defined with longitude increasing to the east from 0 to 360 degrees.
"""
import numpy as np

# ================================================================================= #
# COORDINATE MANAGEMENT                                                             #
# ================================================================================= #

sky_coordinate_systems = {
    "latlon": {"phi": (-np.pi, "east"), "theta": (-np.pi / 2, "north")},
    "healpix": {"phi": (-np.pi, "east"), "theta": (np.pi, "south")},
    "math_convention": {"phi": (0, "east"), "theta": (-np.pi / 2, "north")},
}


def _get_system_offsets(to_system, from_system, dp, dt):
    po, to = sky_coordinate_systems[to_system]["phi"][0] - (
        dp * sky_coordinate_systems[from_system]["phi"][0]
    ), sky_coordinate_systems[to_system]["theta"][0] - (
        dt * sky_coordinate_systems[from_system]["theta"][0]
    )
    return po, to


def _get_system_direction(to_system, from_system):
    pd = (
        1
        if sky_coordinate_systems[to_system]["phi"][1]
        == sky_coordinate_systems[from_system]["phi"][1]
        else -1
    )
    td = (
        1
        if sky_coordinate_systems[to_system]["theta"][1]
        == sky_coordinate_systems[from_system]["theta"][1]
        else -1
    )
    return pd, td


def convert_coordinates(phi, theta, from_system="healpix", to_system="latlon"):
    # -- Sanity Checks -- #
    assert (
        from_system in sky_coordinate_systems
    ), f"Failed to recognize from_system {from_system}."
    assert (
        to_system in sky_coordinate_systems
    ), f"Failed to recognize to_system {to_system}."

    # -- calculating offset differences -- #
    phi_dir, theta_dir = _get_system_direction(to_system, from_system)
    phi_offset, theta_offset = _get_system_offsets(
        to_system, from_system, phi_dir, theta_dir
    )
    return (phi_dir * phi) + phi_offset, (theta_dir * theta) + theta_offset


def convert_skycoord(skycoord, to_system):
    # skycoords always go straight to latlon
    phi, theta = (
        skycoord.frame.spherical.lon.wrap_at("180 deg").rad,
        skycoord.frame.spherical.lat.rad,
    )

    return convert_coordinates(phi, theta, from_system="latlon", to_system=to_system)
