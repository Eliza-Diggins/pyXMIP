"""
HEALPix based skymapping module for storing and working with functions on the sky.
"""
import healpy as hp
import numpy as np
from numbers import Number
import astropy.units as units
from pyXMIP.utilities.core import mainlog, enforce_units, _enforce_style
from pyXMIP.utilities.terminal import Spinner
from astropy.coordinates import SkyCoord
import pathlib as pt
from time import asctime
from astropy.io import fits

def get_healpix_grid(resolution):
    """
    Generate a grid of HEALPix coordinates with the specified resolution.

    Parameters
    ----------
    resolution: :py:class:`astropy.units.Quantity` or Number
        The resolution of the HEALPix grid. If a value is passed with units, the units are assumed to be in ``rad``.

    Returns
    -------
    :py:class:`numpy.array.ndarray`
        Array containing the HEALPix grid coordinates. ``size = (2,n_pixels)``. The first row corresponds to ``theta`` and the
        second to ``phi``.
    """
    if not isinstance(resolution, units.Quantity):
        resolution = resolution * units.rad

    n_sides = int(np.ceil(1/(resolution.to_value(units.rad)*np.sqrt(3))))
    n_pixels = hp.nside2npix(n_sides)
    hpx = np.zeros((2, n_pixels))

    for u in range(n_pixels):
        hpx[:, u] = hp.pix2ang(n_sides, u)

    return hpx

class SkyAtlas:
    r"""
    Class collection of :py:class:`SkyMap` instances.

    .. warning::

        Be cautious about the coordinate conventions for HEALPix grids vs. astronomical coordinate systems. HEALPix grids have coordinates
        from :math:`0<\phi<2\pi` and :math:`0<\theta<\pi`. In all cases, we convert these to values between :math:`-\pi<\phi<\pi` and
        :math:`-\pi/2<\theta<\pi/2`. It is then by this convention that the coordinates are read into any native coordinate system.

        For example, the point :math:`(\pi,0)` on the HEALPix grid (**which the user never needs to interact with**) corresponds to
        :math:`(0,\pi/2)` in the LON/LAT of whichever coordinate system is configured for this particular instance.

    Notes
    -----

    .. admonition:: formatting

        :py:class:`SkyAtlas` objects are wrappers for ``.fits`` files with very specific formats. The Primary HDU for the fits
        file must have the following headers:

        - ``DATE``: creation date.
        - ``MAPS``: Either ``all`` or ``specified``.
          - If ``all``, then every fits image in the file is interpreted as a map.
          - If ``specified``, then every image must have a header ``"ISMAP" = True`` to be read as a map.
        - ``NSIDE``: Integer specifying the number of sides for the HEALPix grid underlying the atlas.
        - ``NPIX``: The number of pixels in the HEALPix grid.
        - ``RES``: The resolution (in radians) of the image.
        - ``CSYS``: The native coordinate system of the skymaps. These should be string representations.

        In addition to the headers, each of the ``MAP`` images in the atlas must have dimensions ``(NPIX,)``.
    """
    _implemented_atlas_formats = ['fits','hdf5']
    _expected_hudl_headers = ['DATE','MAPS','NSIDE',"NPIX","RES",'CSYS']

    def __init__(self,path,coordinate_system='icrs'):
        """
        Initializes the :py:class:`SkyAtlas` instance.

        Parameters
        ----------
        path: str
            The path to the ``.fits`` file from which to read the atlas.
        coordinate_system: str or :py:class:`astropy.coordinates.GenericFrame`, optional
            The coordinate system to use as the native coordinate frame. by default, this is ``"icrs"``.
        """
        mainlog.debug(f"Loading SkyAtlas object from {path}.")

        # -- loading basic attributes -- #
        self.path = pt.Path(path)
        self.coordinate_system = coordinate_system

    def __len__(self):
        return len(self.available_maps())

    def __repr__(self):
        return f"<SkyAtlas [{self.path}]>"

    def __str__(self):
        return f"<SkyAtlas [N={len(self)}]>"

    def __getitem__(self, item):
        return self.fetch_map(item.upper())

    def __delitem__(self, key):
        with fits.open(self.path, 'update') as hudl:
            del hudl[key.upper()]
            hudl.flush()

    def __contains__(self, item):
        return item in self.available_maps()

    def vstack(self,other):
        """
        Stack / conjoin this :py:class:`SkyAtlas` with another one of the same geometry.

        .. warning::

            This operation occurs inplace.

        Parameters
        ----------
        other: :py:class:`SkyAtlas`
            The other atlas to join.

        Returns
        -------
        None
        """
        assert other.header['NPIX'] == self.header['NPIX'], f"Cannot concatenate SkyAtlas instances with incommensurate sizes [{other.header['NPIX']},{self.header['NPIX']}]."

        with fits.open(self.path,'update') as hudl:
            with fits.open(other.path,'update') as hudl_other:

                for entry in hudl_other:
                    if entry in hudl:
                        # This already exists
                        hudl[entry].data += hudl_other[entry].data
                    else:
                        hudl.append(hudl_other[entry])

            hudl.flush()

    @property
    def header(self):
        """
        The header of the underlying fits file. Only the required headers for formatting are returned.

        Returns
        -------
        dict

        """
        with fits.open(self.path) as hudl:
            return {k:hudl[0].header[k] for k in self.__class__._expected_hudl_headers}

    def set_header(self,key,value):
        with fits.open(self.path,'update') as hudl:
            hudl[0].header.set(key.upper(),value)
            hudl.flush()


    def available_maps(self):
        """
        List the names of the available maps.

        Returns
        -------
        list of str
        """
        if self.header['MAPS'] == 'all':
            with fits.open(self.path) as f:
                return [_f.name for _f in f if isinstance(_f,fits.ImageHDU)]
        else:
            return self.header['MAPS']

    def get_coordinates(self):
        """
        Generate the pixel coordinates in the native coordinate system.

        .. note::

            In HEALPix convention, coordinates go from ``[0,2*pi]`` and ``[0,pi]``. The coordinates returned here are **CONVERTED** to
            the standard astronomical convention. NEVER interact directly with the ``healpy`` base; the coordinates will be inconsistent.

        Returns
        -------
        :py:class:`astropy.coordinates.SkyCoord`
        """
        _number_of_pixels = hp.nside2npix(self.header['NSIDE'])
        _coord_array = np.zeros((2, _number_of_pixels))

        for u in range(_number_of_pixels):
            _coord_array[:, u] = hp.pix2ang(self.header['NSIDE'], u)

        return SkyCoord(_coord_array[1,:]-np.pi,np.pi/2 - _coord_array[0,:],frame=self.header['CSYS'],unit='rad')

    @classmethod
    def generate(cls,path,resolution,overwrite=False):
        """
        Create an empty :py:class:`SkyAtlas` of a given resolution.

        Parameters
        ----------
        path: str
            The path to the ``.fits`` file from which this object is to be loaded / saved.
        resolution: :py:class:`astropy.units.Quantity` or Number
            The resolution of the HEALPix grid. If a value is passed with units, the units are assumed to be in ``rad``.
        overwrite: bool
            Allow file overwrite.

        Returns
        -------
        :py:class:`SkyAtlas`
        """
        from astropy.io import fits
        resolution = enforce_units(resolution,units.rad)

        mainlog.info(f"Generating blank SkyAtlas with resolution {resolution} at {path}.")

        # -- resolving the grid -- #
        n_sides = int(np.ceil(1 / (resolution.to_value(units.rad) * np.sqrt(3))))
        n_pixels = hp.nside2npix(n_sides)

        # -- generating the meta data -- #
        header = fits.Header()

        header['DATE'] = asctime()
        header['MAPS'] = 'all'
        header['NSIDE'] = n_sides
        header['NPIX'] = n_pixels
        header['RES'] = resolution.to_value(units.rad)
        header['CSYS'] = 'icrs'

        # -- creating the fits file -- #
        empty_primary = fits.PrimaryHDU(header=header)
        hudl = fits.HDUList([empty_primary])
        hudl.writeto(path,overwrite=overwrite)

        return cls(path)

    def add_map_from_function(self,mappable,name,overwrite=False,*args,**kwargs):
        """
        Add a new :py:class:`SkyMap` to the :py:class:`SkyAtlas` instance using a function.

        Parameters
        ----------
        mappable: callable
            Function of the signature ``func(phi,theta,*args,**kwargs)``. ``phi`` and ``theta`` should be the native coordinates
            of the atlas (as defined by the user).
        name: str
            The name of the map once it has been added to the atlas.

            .. note::

                ``name`` must be all caps to be compliant with the ``fits`` file convention. If this is not done, it will be done
                automatically.

        overwrite: bool, optional
            If ``True``, the map may be overwritten.
        args
            Additional arguments to pass to the map.
        kwargs
            Additional keyword arguments to pass to the map.

        Returns
        -------
        None
        """
        name = name.upper()
        if not overwrite and name in self.available_maps():
            raise ValueError(f"The map {name} already exists and overwrite=False.")

        with Spinner(text=f"Generating map {name}..."):
            coordinates = self.get_coordinates()
            values = mappable(coordinates.frame.spherical.lon.rad,coordinates.frame.spherical.lat.rad,*args,**kwargs)

            if isinstance(values,units.Quantity):
                unit,values = values.unit,values.value
            else:
                unit,values = "",values

            with fits.open(self.path,'update') as hudl:

                image_hdu = fits.ImageHDU(values)
                image_hdu.header['U'] = str(unit)
                image_hdu.header["ISMAP"] = True

                image_hdu.name = name

                if overwrite and name in self.available_maps():
                    del hudl[name]

                hudl.append(image_hdu)
                hudl.flush()

    def fetch_map(self,name):
        """
        Return the :py:class:`SkyMap` corresponding to the provided name.

        Parameters
        ----------
        name: str
            the name of the image to load.

        Returns
        -------
        :py:class:`SkyMap`
        """
        return SkyMap(self.path,name,coordinate_system=self.coordinate_system)

    def to_file(self,path,inplace=False):
        """
        Write this :py:class:`SkyAtlas` to a new file.

        Parameters
        ----------
        path: str
            The path to the new save location.
        inplace: True
            If ``True``, then the object is reloaded from the new source and not returned. Otherwise a new instance is returned.

        Returns
        -------
        :py:class:`SkyAtlas`

        """
        with fits.open(self.path) as hudl:

            hudl.writeto(path)

        if inplace:
            self.path = pt.Path(path)
        else:
            return self.__class__(path)

    def plot2d(self,map,*args,**kwargs):
        """
        Proxy for :py:meth:`SkyMap.plot2d`.
        """
        return self.fetch_map(map).plot2d(*args,**kwargs)

    def plot3d(self,map,*args,**kwargs):
        """
        Proxy for :py:meth:`SkyMap.plot3d`.
        """
        return self.fetch_map(map).plot3d(*args,**kwargs)


class SkyMap:
    """
    Representation of a 1-D HEALPix grid map corresponding to a given atlas file.
    """

    def __init__(self,path,name,coordinate_system='icrs'):
        mainlog.debug(f"Loading SkyMap object from {path} [Name={name}].")

        # -- loading basic attributes -- #
        self.path = pt.Path(path)
        self.name = name.upper()
        self.coordinate_system = coordinate_system

        # -- reading the values -- #
        with fits.open(self.path) as hudl:
            assert name in [u.name for u in hudl], f"There is no SkyMap {name} in {path}."

            self.data = hudl[name].data

        self._npix = len(self.data)
        self._nside = hp.npix2nside(self._npix)

    def __call__(self, position):
        """
        Evaluate the skymap.

        Parameters
        ----------
        position: :py:class:`astropy.coordinates.SkyCoord`

        """
        lon,lat = position.transform_to(self.coordinate_system).frame.spherical.lon.rad,position.transform_to(self.coordinate_system).frame.spherical.lat.rad
        theta,phi = np.pi/2 - lat, lon + np.pi

        pixels = hp.ang2pix(self._nside,theta,phi)

        return self.data[pixels]

    def __str__(self):
        return f"<SkyMap [name={self.name}]>"

    def __repr__(self):
        return f"<SkyMap [name={self.name}] @ {self.path}>"

    @_enforce_style
    def plot3d(
        self,
        ax=None,
        fig=None,
        gridsize=100,
        cmap=None,
        vmin=None,
        vmax=None,
        norm=None,
        colorbar=True,
        flat=True,
        **kwargs,
    ):
        r"""
        Produce a 3D image of the density model with specific arguments or with the set arguments.

        Parameters
        ----------
        arguments: array-like, optional
            Arguments to pass to the density model. If :py:attr:`DensityModel.args` exists, this option is
            overridden by the set arguments.
        ax: :py:class:`matplotlib.pyplot.Axes`, optional
            The axes onto which to plot the data.
        fig: :py:class:`matplotlib.pyplot.Figure`, optional
            The figure in which to embed the image.
        gridsize: int, optional
            The size of the meshgrid to generate the plot from.
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
        flat: bool
            If ``True``, then the map will be flat and colored only on the surface. If ``False``, then the radius of the surface will
            vary with the value.
        **kwargs
            Additional kwargs to pass.

        Returns
        -------
        fig:
            The figure.
        ax:
            The axes.
        cbar:
            The colorbar mappable.
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D


        # ----------------------------------------------#
        # Computing the grid
        phi, theta = np.mgrid[
            0 : 2 * np.pi : gridsize * 1j, -np.pi / 2 : np.pi / 2 : gridsize * 1j
        ]


        pixels = hp.ang2pix(self._nside,theta + np.pi/2,phi)
        image_array = self.data[pixels]

        if flat:
            X, Y, Z = (
                np.cos(theta) * np.cos(phi),
                np.cos(theta) * np.sin(phi),
                np.sin(theta),
            )
        else:
            X, Y, Z = (
                image_array*np.cos(theta) * np.cos(phi),
                image_array*np.cos(theta) * np.sin(phi),
                image_array*np.sin(theta),
            )


        # ------------------------------------------------#
        # Plotting the figure
        if not fig:
            fig = plt.figure(figsize=(5, 4))
            fig.subplots_adjust(left=0, right=0.95, top=0.92, bottom=0.03)
        if not ax:
            ax = fig.add_subplot(111, projection="3d")
        else:
            assert isinstance(ax, Axes3D)

        ax.axis("off")

        if cmap is None:
            cmap = plt.colormaps[plt.rcParams["image.cmap"]]

        if vmin is None:
            vmin = np.amin(image_array)

        if vmax is None:
            vmax = np.amax(image_array)

        if norm is None:
            from matplotlib.colors import Normalize

            norm = Normalize(vmin=vmin, vmax=vmax)

        ax.plot_surface(X, Y, Z, facecolors=cmap(norm(image_array)), **kwargs)
        ax.set_aspect("equal")
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])
        ax.set_title(f"Source Density ({self.__class__.__name__})")
        cbar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        if colorbar:
            cbar = plt.colorbar(cbar, ax=ax)

        return fig, ax, cbar

    @_enforce_style
    def plot2d(
        self,
        fig=None,
        cmap=None,
        vmin=None,
        vmax=None,
        norm=None,
        colorbar=True,
    ):
        r"""
        Produce a 2D image of the density model with specific arguments or with the set arguments.

        .. note::

          This method uses the ``healpy`` package to manage projection onto HEALPix maps.

        Parameters
        ----------

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

        # -- Building the figure -- #
        if not fig:
            fig = plt.figure(figsize=(5, 4))
            fig.subplots_adjust(left=0, right=0.95, top=0.92, bottom=0.03)

        hp.mollview(self.data, fig=fig, cbar=False, notext=True, title=None)
        hp.graticule()

        ax = fig.gca()

        if cmap is None:
            cmap = plt.colormaps[plt.rcParams["image.cmap"]]

        if vmin is None:
            vmin = np.amin(self.data)

        if vmax is None:
            vmax = np.amax(self.data)

        if norm is None:
            from matplotlib.colors import Normalize

            norm = Normalize(vmin=vmin, vmax=vmax)

        cbar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        if colorbar:
            cbar = plt.colorbar(cbar, ax=ax)

        return fig, ax, cbar


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    u = SkyAtlas.generate('test.fits',5*units.deg,overwrite=True)
    u.add_map_from_function(lambda phi,theta: np.cos(phi),"test")
    h = u['test']
    print(h)
