"""
Plotting utilities for use in the pyXMIP backend.
"""
import functools

import numpy as np
from astropy.coordinates import ICRS, Angle
from astropy.units import deg, rad
from matplotlib import pyplot as plt

from pyXMIP.utilities.core import xsparams


# ======================================================================================================================#
# STYLE FUNCTIONS                                                                                                      #
# ======================================================================================================================#
def _enforce_style(func):
    """Enforces the mpl style."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        _rcp_copy = plt.rcParams.copy()

        for _k, _v in xsparams["plotting"]["defaults"].items():
            plt.rcParams[_k] = _v

        out = func(*args, **kwargs)

        plt.rcParams = _rcp_copy
        del _rcp_copy

        return out

    return wrapper


def set_style():
    for _k, _v in xsparams["plotting"]["defaults"].items():
        plt.rcParams[_k] = _v


def _format_hours_minutes_ra_left_longitude(x, pos):
    """Format as H:M:S but with longitude going left."""
    val = Angle(-x * rad).wrap_at(360 * deg)
    return val.to_string(unit="hour", sep=":")


def _format_hours_minutes_ra_right_longitude(x, pos):
    """Format as H:M:S but with longitude going left."""
    val = Angle(x * rad).wrap_at(360 * deg)
    return val.to_string(unit="hour", sep=":")


def _format_deg_minutes_dec_latitude(x, pos):
    """Format as H:M:S but with longitude going left."""
    val = Angle(x * rad).wrap_at(90 * deg)
    return val.to_string(unit="deg", sep=("d", ":"))


tickers = {
    "icrs": {
        "phi_l": _format_hours_minutes_ra_left_longitude,
        "phi_r": _format_hours_minutes_ra_right_longitude,
        "theta": _format_deg_minutes_dec_latitude,
    }
}


# ======================================================================================================================#
# PLOTTING FUNCTIONS                                                                                                    #
# ======================================================================================================================#


def scatter_table_2d(table, *args, **kwargs):
    from astropy import units
    from astropy.coordinates import SkyCoord

    # -------------------------------------------------- #
    # Set up the plot
    # ---------------------------------------------------#
    fig, ax, projection = (
        kwargs.pop("fig", None),
        kwargs.pop("ax", None),
        kwargs.pop("projection", "aitoff"),
    )

    if fig is None:
        fig = plt.figure(figsize=kwargs.pop("figsize", (10, 10)))

    if ax is None:
        ax = fig.add_subplot(111, projection=projection)
    else:
        assert (
            ax.projection == projection
        ), f"The provided axis doesn't have projection {projection}."

    # ---------------------------------------------------- #
    # Manage the coordinates
    # ---------------------------------------------------- #
    phi_column, theta_column, coordinate_system, plot_coord_system = (
        kwargs.pop("phi_col", "RA"),
        kwargs.pop("theta_col", "DEC"),
        kwargs.pop("coordinate_system", ICRS),
        kwargs.pop("plot_coord_system", ICRS),
    )

    positions = SkyCoord(
        table[phi_column], table[theta_column], frame=coordinate_system
    )

    # convert the positions to the desired coordinate system
    positions = positions.transform_to(plot_coord_system)
    phi, theta = (
        positions.frame.spherical.lon.wrap_at(180 * units.deg).rad,
        positions.frame.spherical.lat.rad,
    )

    artist = ax.scatter(phi, theta, *args, **kwargs)

    return fig, ax, artist


def scatter_table_3d(table, *args, **kwargs):
    from astropy import units
    from astropy.coordinates import SkyCoord

    # -------------------------------------------------- #
    # Set up the plot
    # ---------------------------------------------------#
    fig, ax = kwargs.pop("fig", None), kwargs.pop("ax", None)

    if fig is None:
        fig = plt.figure(figsize=kwargs.pop("figsize", (10, 10)))

    if ax is None:
        ax = fig.add_subplot(111, projection="3d")
    else:
        assert ax.projection in [
            "3d",
            "3D",
        ], "The provided axis doesn't have projection 3d."

    # ---------------------------------------------------- #
    # Manage the coordinates
    # ---------------------------------------------------- #
    phi_column, theta_column, coordinate_system, plot_coord_system = (
        kwargs.pop("phi_col", "RA"),
        kwargs.pop("theta_col", "DEC"),
        kwargs.pop("coordinate_system", ICRS),
        kwargs.pop("plot_coord_system", ICRS),
    )

    positions = SkyCoord(
        table[phi_column], table[theta_column], frame=coordinate_system
    )

    # convert the positions to the desired coordinate system
    positions = positions.transform_to(plot_coord_system)
    phi, theta = (
        positions.frame.spherical.lon.wrap_at(360 * units.deg).rad,
        positions.frame.spherical.lat.rad,
    )

    X = np.cos(phi) * np.cos(theta)
    Y = np.sin(phi) * np.cos(theta)
    Z = np.sin(theta)

    artist = ax.scatter3D(X, Y, Z, *args, **kwargs)

    return fig, ax, artist
