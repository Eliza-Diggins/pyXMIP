"""
Plotting utilities for use in the pyXMIP backend.
"""
import functools

import healpy as hp
import healpy.projaxes as PA
import numpy as np
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


# ===================================================================================================================== #
# Minor Process Functions                                                                                               #
# ===================================================================================================================== #
_projection_classes = {"mollweide": PA.HpxMollweideAxes}


def _check_healpix(map):
    hp.pixelfunc.check_nside(hp.pixelfunc.get_nside(map))

    return hp.pixelfunc.get_nside(map)


# ======================================================================================================================#
# Plotting Systems                                                                                                      #
# ======================================================================================================================#


def plot_healpix(healpix_map, fig=None, ax=None, projection=None, **kwargs):
    # -------------------------------------------- #
    # Map checking and setup
    # -------------------------------------------- #
    # managing the projection type.
    if projection is None:
        projection = "mollweide"
    projection_class = _projection_classes[projection]

    # Setting up the figure and axes
    if fig is None:
        fig = plt.figure(figsize=kwargs.pop("figsize", (8, 4)))
    if ax is None:
        ax = projection_class(
            fig,
            [0.025, 0.025, 0.95, 0.95],
            coord=kwargs.pop("coord", None),
            rot=kwargs.pop("rot", None),
            format=kwargs.pop("format", "%g"),
            flipconv=kwargs.pop("flip", "geo"),
        )
    elif not isinstance(ax, projection_class):
        # ax is a class, but not the right one.
        _old_axes_position = ax._position
        fig.remove(ax)
        del ax

        # create the new axes.
        ax = projection_class(
            fig,
            _old_axes_position,
            coord=kwargs.pop("coord", None),
            rot=kwargs.pop("rot", None),
            format=kwargs.pop("format", "%g"),
            flipconv=kwargs.pop("flip", "geo"),
        )
    else:
        pass

    fig.add_axes(ax)
    # -------------------------------------------- #
    # Adding the map to the axes.
    # -------------------------------------------- #
    img = ax.projmap(
        healpix_map,
        vmin=kwargs.pop("vmin", np.amin(healpix_map)),
        vmax=kwargs.pop("vmax", np.amax(healpix_map)),
        cmap=kwargs.pop("cmap", "viridis"),
        badcolor=kwargs.pop("fillna", "black"),
        bgcolor=kwargs.pop("facecolor", "w"),
        norm=kwargs.pop("norm", None),
        alpha=kwargs.pop("alpha", None),
    )

    return fig, ax, img


if __name__ == "__main__":
    from pyXMIP.structures.databases import NED

    q = NED.get_poisson_atlas()

    p = q.get_map("IrS")
    h = p.data

    plot_healpix(h, cmap="gnuplot")
    plt.show()
