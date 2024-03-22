"""
Plotting utilities for use in the pyXMIP backend.
"""
import functools

from matplotlib import pyplot as plt
import numpy as np
from pyXMIP.utilities.core import xsparams
from astropy.coordinates import ICRS, Angle
from astropy.units import rad,deg
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
    val = Angle(-x*rad).wrap_at(360*deg)
    return val.to_string(unit='hour',sep=":")

def _format_hours_minutes_ra_right_longitude(x, pos):
    """Format as H:M:S but with longitude going left."""
    val = Angle(x*rad).wrap_at(360*deg)
    return val.to_string(unit='hour',sep=":")

def _format_deg_minutes_dec_latitude(x, pos):
    """Format as H:M:S but with longitude going left."""
    val = Angle(x*rad).wrap_at(90*deg)
    return val.to_string(unit='deg', sep=('d', ':'))

tickers = {
    "icrs": {
        "phi_l": _format_hours_minutes_ra_left_longitude,
        "phi_r": _format_hours_minutes_ra_right_longitude,
        "theta": _format_deg_minutes_dec_latitude
    }
}


# ======================================================================================================================#
# PLOTTING FUNCTIONS                                                                                                   #
# ======================================================================================================================#
def imshow_healpix(image_data,
                   pixel_positions,
                   fig=None,
                   cmap=None,
                   vmin=None,
                   vmax=None,
                   norm=None,
                   colorbar=True,
                   coordinate_system=None):
    r"""
    Produce a 2D image on a HEALPix grid

    .. note::

      This method uses the ``healpy`` package to manage projection onto HEALPix maps.

    Parameters
    ----------

    image_data: array-like
        This is the 1-D array of values to generate the image from. It should have the same length as the relevant HEALPix grid.
    pixel_positions: SkyCoord
        The SkyCoord array containing the pixel positions in the preferred coordinate system.
    fig: :py:class:`matplotlib.pyplot.Figure`, optional
        The figure in which to embed the image.
    cmap: :py:class:`matplotlib.colors.Colormap`, optional
        The colormap to use.
    vmin: float, optional
        The minimum color value.
    vmax: float, optional
        The maximum color value.
    norm: :py:class:`matplotlib.colors.ColorNorm`, optional
        The normalization for the colorbar.
    colorbar: bool, optional
        If ``True``, a colorbar is added to the axes. Otherwise the :py:class:`matplotlib.pyplot.cm.ScalarMappable` corresponding
        to the colormap is returned.
    coordinate_system: :py:class:`astropy.coordinates.Frame`
        The coordinate system in which to plot the outputs.

    Returns
    -------
    fig:
        The figure.
    ax:
        The axes.
    cbar:
        The colorbar mappable.
    """
    # -- check that healpy is installed -- #
    try:
        import healpy as hp
    except ImportError:
        raise ImportError(
            "Cannot plot a 2d projection plot without healpy installed."
        )

    # -- Managing positions and building healpix specific data -- #
    _n_pixels = len(image_data) # The number of pixels in the HEALPix grid.

    if coordinate_system is not None:
        pixel_positions = pixel_positions.transform_to(coordinate_system)

    _im_phi,_im_theta = pixel_positions.frame.spherical.lon.wrap_at(2*np.pi*rad).rad, pixel_positions.frame.spherical.lon.wrap_at(np.pi*rad).rad
    _im_theta = np.pi/2 - _im_theta # Fix the awkward covention for HEALPix

    # -- Building the figure -- #
    if not fig:
        fig = plt.figure(figsize=(5, 4))
        fig.subplots_adjust(left=0, right=0.95, top=0.92, bottom=0.03)

    hp.mollview(image_data, fig=fig, cbar=False, notext=True, title=None)
    hp.graticule()

    ax = fig.gca()

    if cmap is None:
        cmap = plt.colormaps[plt.rcParams["image.cmap"]]

    if vmin is None:
        vmin = np.amin(image_data)

    if vmax is None:
        vmax = np.amax(image_data)

    if norm is None:
        from matplotlib.colors import Normalize

        norm = Normalize(vmin=vmin, vmax=vmax)

    cbar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    if colorbar:
        cbar = plt.colorbar(cbar, ax=ax)

    return fig, ax, cbar


def _generate_fig_axes(fig=None,ax=None,*args,**kwargs):
    """Generate the figure and axes for a plot."""
    projection = kwargs.pop('projection',None)
    figuresize = kwargs.pop('figsize',(10,10))

    if not fig:
        fig = plt.figure(figsize=figuresize)

    if not ax:
        if projection is None:
            ax = fig.add_subplot(111,*args,**kwargs)
        else:
            ax = fig.add_subplot(111,projection=projection,*args,**kwargs)
    else:
        assert ax.projection == projection

    return fig, ax



if __name__ == '__main__':
    a = [1, 2, 3]
