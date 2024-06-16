"""
Plotting utilities for use in the pyXMIP backend.
"""
import functools
from numbers import Number
from typing import Sequence

import healpy.projaxes as PA
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.io.fits import HDUList, Header
from matplotlib import pyplot as plt

from pyXMIP.utilities.core import pxconfig
from pyXMIP.utilities.logging import devlog


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
# Image Processing and Convenience Functions                                                                            #
# ===================================================================================================================== #
_projection_classes = {
    "mollweide": PA.HpxMollweideAxes,
    "gnomonic": PA.HpxGnomonicAxes,
    "orthographic": PA.HpxOrthographicAxes,
}


def replace_axes(ax: plt.Axes, fig: plt.Figure = None, *args, **kwargs):
    # -- Determine the most fitting figure -- #
    if fig is None:
        fig = plt.gcf()

    # Get the position of the axes object #
    position = ax.get_position(original=True)
    ax.remove()
    del ax

    return fig.add_axes(position, *args, **kwargs)


def image_histogram_equalization(
    image: np.ndarray,
    bins: int | Sequence[Number] = 256,
    vmin: float = 0.0,
    vmax: float = 1.0,
) -> np.ndarray:
    """
    Equalize an image's color so that the PDF is flattened.

    In effect, this means that no particular pixel brightness with be preferred and there will be as many light pixels as
    dark pixels. Typically, this improves visibility but leads to non-tractable data values in the output array.

    Parameters
    ----------
    image: array
        The array representing the image. Must be a 2-D array of floats.
    bins: int or list of float, optional
        The bins to use for the CDF.
    vmin: float, optional
        Value from 0 to 1 for the minimum cutoff of the CDF.

        This causes the ``vmin`` dimmest pixels to be fixed at the value brightness of the dimmest bin above. This can be used
        to remove noise from images, force dim components of the image to become less pronounced, etc. The default is ``0.0``.
    vmax: float, optional
        The maximum value from 0 to 1 for the cutoff of the CDF.

        This causes the ``vmax`` brightest pixels to be fixed at the brightness of the next dimmest bin below. This can be used
        to remove sharp spikes from the image and to improve visibility of dimmer parts of the image. default is ``1.0``.

    Returns
    -------
    array
        The output image array.
    """
    # -- Construct the histogram and CDF -- #
    image_histogram, bins = np.histogram(image.flatten(), bins=bins, density=True)
    cdf = image_histogram.cumsum()  # cumulative distribution function
    cdf = cdf / cdf[-1]  # normalize

    # Apply cutoffs #
    bin_centers = (np.array(bins[1:]) + np.array(bins[:-1])) / 2
    vmin_bin, vmax_bin = bin_centers[cdf >= vmin][0], bin_centers[cdf <= vmax][-1]

    image[image <= vmin_bin] = vmin_bin
    image[image >= vmax_bin] = vmax_bin

    # -- Interpolate -- #
    # use linear interpolation of cdf to find new pixel values
    image_equalized = np.interp(image.flatten(), bin_centers, cdf)

    return image_equalized.reshape(image.shape)


def get_hips_data(
    center: SkyCoord,
    FOV: Angle | str,
    dims: Sequence[Number] = (100, 100),
    hips_path: str = pxconfig.config.plotting.hips_defaults.hips_map,
    projection: str = pxconfig.config.plotting.hips_defaults.projection,
    **kwargs,
) -> HDUList:
    """
    Use ``astroquery`` to seek a HIPs image and return it with the relevant parameters.

    Parameters
    ----------
    center: :py:class:`astropy.coordinates.SkyCoord`
        The coordinates corresponding to the center of the image.
        These must be valid astronomical coordinates.
    FOV: :py:class:`astropy.coordinates.Angle` or str
        The field of view for the image. This will be automatically converted to the :py:class:`astropy.coordinates.Angle` class,
        regardless of the input type. If conversion fails to occur, an error is raised.
    dims: tuple of int
        The dimensions of the returned image.
    hips_path: str, optional
        The server path to the HIPs image desired.
    projection: str, optional
        The default return projection for the output image.

    """
    from astropy.coordinates import Angle, Latitude, Longitude
    from astroquery.hips2fits import hips2fits

    devlog.debug(f"Querying for HIPs at {hips_path}.")
    # ---------------------------------------------------------- #
    # Constructing the parameters                                #
    # ---------------------------------------------------------- #
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


def get_hips_image(
    center: SkyCoord,
    FOV: Angle | str,
    dims: Sequence[Number] = (500, 500),
    eq_kwargs: dict = None,
    equalize: bool = True,
    **kwargs,
) -> tuple[np.ndarray, Header]:
    # ---------------------------------------------------------- #
    # Constructing the parameters                                #
    # ---------------------------------------------------------- #
    if eq_kwargs is None:
        eq_kwargs = {}

    # ---------------------------------------------------------- #
    # Query the image data                                       #
    # ---------------------------------------------------------- #
    image_data = get_hips_data(center, FOV, dims, **kwargs)

    # The image is always returned as a fits HDUList, but may be different shapes.
    header = image_data[0].header
    data = np.array(image_data[0].data, dtype=float)

    if data.ndim == 3:
        # This is a complete image as is.
        return data.T, header

    elif data.ndim == 2:
        # This is a value array but not an image yet.
        if equalize:
            data = image_histogram_equalization(
                data,
                bins=eq_kwargs.pop("eq_bins", 256),
                vmin=eq_kwargs.pop("eq_vmin", 0.0),
                vmax=eq_kwargs.pop("eq_vmax", 1.0),
            )

        cmap = eq_kwargs.pop(
            "cmap", getattr(plt.cm, pxconfig.config.plotting.hips_defaults.cmap)
        )

        image = cmap(data)

        return image, header

    else:
        raise ValueError(
            f"Returned image array had shape {data.shape}, which had a non-standard number of dimensions ({data.ndim})."
        )


# ======================================================================================================================#
# Plotting Systems                                                                                                      #
# ======================================================================================================================#


def plot_healpix(
    healpix_map: np.ndarray,
    fig: plt.Figure = None,
    ax: plt.Axes = None,
    projection: str = None,
    **kwargs,
):
    """
    Plot a HEALPix map.

    Parameters
    ----------
    healpix_map: array
        The HEALPix map to plot.
    fig: :py:class:`matplotlib.Figure`
        The figure to add the axes to.
    ax: :py:class:`matplotlib.Axes`
        The axes to add the image to.
    projection: str
        The projection to select.
    kwargs
        Additional kwargs.

    Returns
    -------
    fig
        The output figure
    ax
        The output axes
    img
        The output image.
    """
    # -------------------------------------------- #
    # Map checking and setup                       #
    # -------------------------------------------- #
    # managing the projection type.
    if projection is None:
        projection = pxconfig.plotting.healpix_defaults.projection
    projection_class = _projection_classes[projection]

    # Setting up the figure and axes
    if fig is None:
        fig = plt.figure()
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
    center: SkyCoord,
    FOV: Angle | str,
    fig: plt.Figure = None,
    ax: plt.Axes = None,
    imshow_kwargs: dict = None,
    **kwargs,
):
    """
    Create a HIPs map figure.

    Parameters
    ----------
    center: :py:class:`astropy.coordinates.SkyCoord`
        The center of the image.
    FOV: :py:class:`astropy.coordinates.Angle`
        The image field of view.
    fig: figure, optional
        The matplotlib figure to create the artist inside of.
    ax: axes: optional
        The matplotlib artist to display the image.

        .. warning::

            Because the axes for the HIP map must be WCS, if an existing
            ``ax`` is passed, it will be removed and replaced with a new one with correct
            projection and in the same position.

            Thus, if you are adding multiple things to the plot, you are advised to start with the HIPs layer
            and then add things to the returned axes.
    imshow_kwargs: dict, optional
        Additional kwargs to pass to ``plt.imshow``.
    kwargs:
        Additional kwargs passed through the :py:func:`get_hips_image` function.

    Returns
    -------
    fig: figure
        The figure with the resulting image.
    ax: axes
        The axes object returned.
    """
    # -- Managing WCS axes -- #
    from astropy.wcs import WCS

    # ---------------------------------------- #
    # Pull the HIPs image                      #
    # ---------------------------------------- #
    image, header = get_hips_image(center, FOV, **kwargs)

    # -- Throw the WCS coordinates -- #
    wcs = WCS(header)

    # ========================================== #
    # Setup the figure and axes as needed        #
    # ========================================== #
    if fig is None:
        fig = plt.figure()

    if ax is None:
        ax = fig.add_subplot(111, projection=wcs)
    else:
        ax = replace_axes(ax, fig=fig, projection=wcs)

    if imshow_kwargs is None:
        imshow_kwargs = {}

    ax.imshow(image, origin="lower", **imshow_kwargs)

    return fig, ax
