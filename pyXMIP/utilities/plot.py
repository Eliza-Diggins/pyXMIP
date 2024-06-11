"""
Plotting utilities for use in the pyXMIP backend.
"""
import functools

import healpy as hp
import healpy.projaxes as PA
import numpy as np
from matplotlib import pyplot as plt

from pyXMIP.utilities.core import pxconfig


# ======================================================================================================================#
# STYLE FUNCTIONS                                                                                                      #
# ======================================================================================================================#
def _enforce_style(func):
    """Enforces the mpl style."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        _rcp_copy = plt.rcParams.copy()

        for _k, _v in pxconfig["plotting"]["defaults"].items():
            plt.rcParams[_k] = _v

        out = func(*args, **kwargs)

        plt.rcParams = _rcp_copy
        del _rcp_copy

        return out

    return wrapper


def set_style():
    """
    Uses the ``pyXMIP`` settings in the configuration file and enforces those formatting choices on matplotlib.
    """
    for _k, _v in pxconfig["plotting"]["defaults"].items():
        plt.rcParams[_k] = _v


# ===================================================================================================================== #
# Minor Process Functions                                                                                               #
# ===================================================================================================================== #
_projection_classes = {"mollweide": PA.HpxMollweideAxes}


def _check_healpix(map):
    hp.pixelfunc.check_nside(hp.pixelfunc.get_nside(map))

    return hp.pixelfunc.get_nside(map)


def get_hips_image(
    center, FOV, dims, hips_path="CDS/P/DSS2/red", projection="AIT", **kwargs
):
    """
    Fetch a HIPs image from arbitrary servers and return.

    Parameters
    ----------
    center: :py:class:`astropy.coordinates.SkyCoord`
        The coordinates corresponding to the center of the image.
    FOV: :py:class:`astropy.coordinates.Angle` or str or float
        The field of view for the image. This will be automatically converted the :py:class:`astropy.coordinates.Angle` regardless
        of the input type. If conversion fails to occur, an error is raised.
    dims: tuple of int
        The dimensions of the returned image.
    hips_path: str, optional
        The server path to the HIPs image desired.
    projection: str, optional
        The default return projection for the output image.

    """
    from astropy.coordinates import Angle, Latitude, Longitude
    from astroquery.hips2fits import hips2fits

    # ============================ #
    # Setup / build WCS
    # ============================ #
    parameters = {
        "hips": hips_path,
        "ra": Longitude(center.transform_to("icrs").ra),
        "dec": Latitude(center.transform_to("icrs").dec),
        "width": dims[0],
        "height": dims[1],
        "fov": Angle(FOV),
        "format": "fits",
        "projection": projection,
        **kwargs,
    }

    return hips2fits.query(**parameters)


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


def plot_hips(
    center,
    FOV,
    fig=None,
    projection=None,
    resolution=1000,
    hips_path="CDS/P/DSS2/red",
    hips_kwargs=None,
    **kwargs,
):
    # -- Managing WCS axes -- #
    from astropy.wcs import WCS

    if projection is None:
        projection = "AIT"

    # ========================================== #
    # Pull Image                                 #
    # ========================================== #
    if hips_kwargs is None:
        hips_kwargs = {}

    image = get_hips_image(
        center,
        FOV,
        (resolution, resolution),
        hips_path=hips_path,
        projection=projection,
        **hips_kwargs,
    )

    # -- Throw the WCS coordinates -- #
    wcs = WCS(image[0].header)

    # ========================================== #
    # Setup the figure and axes as needed        #
    # ========================================== #
    if fig is None:
        fig = plt.figure(figsize=kwargs.pop("figsize", (6, 6)))

    ax = fig.add_subplot(111, projection=wcs)
    ax.imshow(image[0].data, **kwargs)

    return fig, ax, wcs
