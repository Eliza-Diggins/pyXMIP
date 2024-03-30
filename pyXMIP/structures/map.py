"""
Skymapping handlers and classes for pyXMIP.
"""
import pathlib as pt
from time import asctime

import healpy as hp
import numpy as np
from astropy import coordinates as astro_coords
from astropy import units
from astropy.io import fits
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from pyXMIP.utilities.core import enforce_units, mainlog


class MapAtlas:
    r"""
    The :py:class:`MapAtlas` class is a generic wrapper for ``.fits`` files which is designed to store sky maps and associated
    data tables.

    Notes
    -----

    On it's surface, :py:class:`MapAtlas` is a standard fits file; except for a few special details regarding the headers and the
    data formats. In the primary HDU, the fits file has additional parameters which determine the behavior of the underlying sky map geometry.

    - ``NPIX``: The number of sky-pixels in the underlying HEALPix grid.
    - ``NSIDE``: Equivalent to ``NPIX`` - specifies the HEALPix grid.
    - ``CSYS``: The native coordinate system of the map.
    - ``CDATE``: The date the :py:class:`MapAtlas` object was created.
    - ``EDATE``: The last date the :py:class:`MapAtlas` was edited.
    - ``RES``: The resolution of the map.

    .. note::

        HEALPix maps use a different coordinate convention than is standard for spherical coordinate systems. Generically, HEALPix uses
        coordinates :math:`0\le \phi < 2\pi` and :math:`0 \le \theta < \pi`, with :math:`\theta` measured from the north-pole.

        The user **NEVER** interacts with this coordinate system when interacting with package objects. The coordinates are automatically converted
        to the standard convention for lon / lat: :math:`[0,2\pi] \times [-\pi/2,\pi/2]`. These are the so-called ``base_coordinates`` of the
        HEALPix grid and are not yet tied to any physical context. The ``CSYS`` header entry specifies the :py:mod:`astropy.coordinates` frame
        on which those coordinates are to be interpreted.

    In addition to the Primary HDU, which has the aforementioned special headers, :py:class:`MapAtlas` instances can contain
    very special Image HDU's called Map HDU's. These are ``1xNPIX`` arrays representing some function on the HEALPix grid.

    The headers of Map HDU's need to be identified in their header with ``ISMAP = True``.
    """

    def __init__(self, filepath):
        """
        Initialize the :py:class:`MapAtlas` from a specified ``.fits`` file.

        Parameters
        ----------
        filepath: str
            The path to the ``.fits`` file to read from.

        """
        self.path = filepath
        # ==========================================#
        # Loading the fits information from header
        # ==========================================#
        with fits.open(self.path, "update") as hudl:
            self._npix = hudl[0].header["NPIX"]
            self._nside = hudl[0].header["NSIDE"]
            self._csys = hudl[0].header["CSYS"]
            self._cd, self._ed = hudl[0].header["CDATE"], hudl[0].header["EDATE"]

    def _set_header(self, key, value):
        assert key not in [
            "NPIX",
            "NSIDE",
        ], "NPIX and NSIDE are invariant structures and cannot be edited."

        with fits.open(self.path, "update") as hudl:
            hudl[0].header[key] = value

            hudl.flush()

    def _get_header(self, key):
        try:
            return getattr(self, f"_{key.lower()}")
        except AttributeError:
            with fits.open(self.path, "update") as hudl:
                setattr(self, f"_{key.lower()}", hudl[0].header[key.upper()])
            return getattr(self, f"_{key.lower()}")

    @property
    def NPIX(self):
        """The number of pixels in the corresponding HEALPix grid."""
        return self._get_header("NPIX")

    @property
    def NSIDE(self):
        """The NSIDES parameter for the HEALPix grid. Equivalent to NPIX."""
        return self._get_header("NSIDE")

    @property
    def CSYS(self):
        """String representation of the standard coordinate system."""
        return self._get_header("CSYS")

    @property
    def CDATE(self):
        """Creation date for the Atlas."""
        return self._get_header("CDATE")

    @property
    def RES(self):
        """The resolution of the map"""
        return self._get_header("RES")

    @property
    def EDATE(self):
        """Last edit date for the Atlas."""
        return self._get_header("EDATE")

    @NPIX.setter
    def NPIX(self, value):
        return self._set_header("NPIX", value)

    @NSIDE.setter
    def NSIDE(self, value):
        self._set_header("NSIDE", value)

    @CSYS.setter
    def CSYS(self, value):
        self._set_header("CSYS", value)

    @CDATE.setter
    def CDATE(self, value):
        self._set_header("CDATE", value)

    @EDATE.setter
    def EDATE(self, value):
        self._set_header("EDATE", value)

    @RES.setter
    def RES(self, value):
        self._set_header("RES", value)

    @property
    def map_names(self):
        """
        Returns a list of the available map HDU names.

        Returns
        -------
        list
        """
        with fits.open(self.path, "update") as hudl:
            return [
                q.name
                for q in hudl
                if isinstance(q, fits.ImageHDU) and q.header["ISMAP"]
            ]

    @property
    def hdus(self):
        """
        Return a list of all of the available HDUs.

        Returns
        -------
        list
        """
        with fits.open(self.path, "update") as hudl:
            return [u.name for u in hudl]

    @property
    def coordinate_frame(self):
        """The :py:class:`astropy.coordinates.GenericFrame` for the Atlas's coordinate system"""
        return getattr(astro_coords, self.CSYS)

    @property
    def pixel_positions(self):
        """The SkyCoord positions of the healpix pixels."""
        _ph, _th = hp.pix2ang(self.NSIDE, np.arange(self.NPIX))
        _ph, _th = _healpix_coordinates_to_spherical(_ph, _th)
        return astro_coords.SkyCoord(_ph, _th, frame=self.coordinate_frame, unit="rad")

    @classmethod
    def generate(cls, path, resolution, overwrite=False):
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

        resolution = enforce_units(resolution, units.rad)

        mainlog.info(
            f"Generating blank SkyAtlas with resolution {resolution} at {path}."
        )

        # -- resolving the grid -- #
        n_sides = int(np.ceil(1 / (resolution.to_value(units.rad) * np.sqrt(3))))
        n_pixels = hp.nside2npix(n_sides)

        # -- generating the meta data -- #
        header = fits.Header()

        header["CDATE"] = asctime()
        header["NSIDE"] = n_sides
        header["NPIX"] = n_pixels
        header["RES"] = resolution.to_value(units.rad)
        header["CSYS"] = "ICRS"

        # -- creating the fits file -- #
        empty_primary = fits.PrimaryHDU(header=header)
        hudl = fits.HDUList([empty_primary])
        hudl.writeto(path, overwrite=overwrite)

        return cls(path)

    def reshape_healpix(self, resolution, force=False):
        """
        Reshape the underlying HEALPix grid for the data.

        .. warning::

            This will result in the deletion of the maps in the Atlas.

        Parameters
        ----------
        resolution: units.Quantity
            The resolution radius of the new HEALPix grid.
        force: bool, optional
            If ``True``, force the deletion of maps instead of failing when they are encountered.

        Returns
        -------
        None
        """
        mainlog.info(
            f"Reshaping HEALPix grid to resolution {resolution} [{self.path}]."
        )

        # ----------------------------------------------------------#
        # Managing existing HEALPix grids that need to be replaced
        # ----------------------------------------------------------#
        if not force and self.has_maps:
            raise ValueError(
                "Maps already exist in this atlas. They will be removed if you proceed. To proceed use force=True."
            )
        elif self.has_maps:
            mainlog.warning(
                f"Deleted {len(self.map_names)} maps from {self.path} to change HEALPix size."
            )
            with fits.open(self.path, "update") as hudl:
                for name in self.map_names:
                    del hudl[name]

                hudl.flush()
        else:
            pass
        # ----------------------------------------------------------#
        # Changing the HEALPix geometry.
        # ----------------------------------------------------------#
        n_sides = int(np.ceil(1 / (resolution.to_value(units.rad) * np.sqrt(3))))
        n_pixels = hp.nside2npix(n_sides)

        with fits.open(self.path, "update") as hudl:
            hudl[0].header["NSIDE"] = n_sides
            hudl[0].header["NPIX"] = n_pixels
            hudl[0].header["RES"] = resolution.to_value(units.rad)
            hudl[0].header["EDATE"] = asctime()

            hudl.flush()


class StatAtlas(MapAtlas):
    """
    The :py:class:`StatAtlas` object is a subclass of :py:class:`MapAtlas` which wraps a ``.fits`` file containing (in addition to the
    standard maps), a table of random samples from a specified database.

    Notes
    -----

    On it's surface, :py:class:`MapAtlas` is a standard fits file; except for a few special details regarding the headers and the
    data formats. In the primary HDU, the fits file has additional parameters which determine the behavior of the underlying sky map geometry.

    - ``NPIX``: The number of sky-pixels in the underlying HEALPix grid.
    - ``NSIDE``: Equivalent to ``NPIX`` - specifies the HEALPix grid.
    - ``CSYS``: The native coordinate system of the map.
    - ``CDATE``: The date the :py:class:`MapAtlas` object was created.
    - ``EDATE``: The last date the :py:class:`MapAtlas` was edited.
    - ``RES``: The resolution of the map.
    - ``DBNAME``: The name of the database.
    """

    def _set_header(self, key, value):
        assert key not in [
            "NPIX",
            "NSIDE",
        ], "NPIX and NSIDE are invariant structures and cannot be edited."

        with fits.open(self.path, "update") as hudl:
            hudl[0].header[key] = value

            hudl.flush()

    def _build_blank_healpix_maps(self):
        """build out blank HEALPix maps for each of the object types available."""
        pass

    def __init__(self, filepath):
        super().__init__(filepath)

        with fits.open(self.path, "update") as hudl:
            self._database = hudl[0].header["DBNAME"]

    @property
    def database(self):
        """The name of the database class."""
        return self._get_header("DBNAME")

    @database.setter
    def database(self, value):
        self._set_header("DBNAME", value)

    @classmethod
    def generate(cls, path, resolution, overwrite=False):
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

        resolution = enforce_units(resolution, units.rad)

        mainlog.info(
            f"Generating blank SkyAtlas with resolution {resolution} at {path}."
        )

        # -- resolving the grid -- #
        n_sides = int(np.ceil(1 / (resolution.to_value(units.rad) * np.sqrt(3))))
        n_pixels = hp.nside2npix(n_sides)

        # -- generating the meta data -- #
        header = fits.Header()

        header["CDATE"] = asctime()
        header["NSIDE"] = n_sides
        header["NPIX"] = n_pixels
        header["RES"] = resolution.to_value(units.rad)
        header["CSYS"] = "ICRS"
        header["DBNAME"] = "NONE"
        header["EDATE"] = asctime()

        # -- creating the fits file -- #
        empty_primary = fits.PrimaryHDU(header=header)
        hudl = fits.HDUList([empty_primary])
        hudl.writeto(path, overwrite=overwrite)

    @property
    def has_maps(self):
        if len(self.map_names):
            return True
        else:
            return False

    def get_points(self):
        import warnings

        from astropy.table import Table

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with fits.open(self.path) as hudl:
                _out = Table(hudl["COUNTS"].data)

        # determining the HEALPix grid
        # !We ALWAYS write counts in RA/DEC for simplicity.
        count_positions = astro_coords.SkyCoord(
            ra=_out["RA"], dec=_out["DEC"], unit="deg"
        )

        _p, _t = (
            count_positions.ra.wrap_at(360 * units.deg).rad,
            count_positions.dec.rad,
        )
        _p, _t = _spherical_coordinates_to_healpix(_p, _t)

        _out["PIX_ID"] = hp.ang2pix(self.NSIDE, _t, _p)

        return _out

    def build_poisson_maps(
        self, method="MLE", prior=None, multiprocess=False, multiprocess_kw=None
    ):
        mainlog.info(
            f"Building poisson density maps for {self.NPIX} HEALPix points from Atlas at {self.path}."
        )

        # ====================================================#
        # Preparing
        # ====================================================#
        if multiprocess:
            mainlog.debug("Constructing poisson maps with multiprocessing enabled.")
        else:
            mainlog.debug("Constructing poisson maps with multiprocessing disabled.")

            with logging_redirect_tqdm(loggers=[mainlog]):
                for _ in tqdm(np.arange(self.NPIX)):
                    pass

    def build_poisson_map(self, object_type, *args, **kwargs):
        mainlog.info(f"Building poisson map in {self.path} for {object_type}.")

        # --------------------------------------------------#
        # Setup and argument management
        # --------------------------------------------------#
        # manage arguments
        mode = kwargs.pop("mode", "LOCAL_UNIFORM")

        # pull points
        point_table = self.get_points()[[object_type, "PIX_ID", "RAD"]]
        point_table["AREA"] = np.pi * point_table["RAD"] ** 2

        # --------------------------------------------------#
        # Building counts information
        # --------------------------------------------------#
        if mode in ["GLOBAL_UNIFORM", "GU", "gu", "global_uniform"]:
            # use a global uniform method
            map = self._bpm_gu(point_table, object_type, *args, **kwargs)
        elif mode in ["LOCAL_UNIFORM", "LU", "lu", "local_uniform"]:
            map = self._bpm_lu(point_table, object_type, *args, **kwargs)
        elif mode in ["kernel"]:
            map = None
        else:
            raise ValueError(f"The method {mode} is not recognized.")

        # ------------------------------------------------#
        # Fix Broken
        # ------------------------------------------------#
        fixna = kwargs.pop("fixna", None)

        if fixna == "average":
            map[np.isnan(map)] = np.sum(point_table[object_type]) / np.sum(
                point_table["AREA"]
            )
        elif fixna == "zero":
            map[np.isnan(map)] = 0
        else:
            pass

        return map

    def _bpm_gu(self, point_table, object_type, *args, **kwargs):
        mainlog.debug("Using GLOBAL_UNIFORM to determine poisson map.")
        # -- Manage arguments and kwargs -- #
        prior = kwargs.pop("prior", None)

        # -- computing globals -- #
        count, area = np.sum(point_table[object_type]), np.sum(point_table["AREA"])

        # -- break method on type -- #
        if prior is None:
            mainlog.debug("No prior detected: MLE estimator (analytic)")
            # no prior was specified; we can simply proceed via MLE.
            _map_value = count / area
        elif callable(prior):
            mainlog.debug("Prior detected: Using a MAP scheme.")
            # the prior is callable
            from scipy.optimize import minimize
            from scipy.stats import poisson

            _min_func = lambda x: -poisson.pmf(count, area * x) * prior(x)
            _map_value = minimize(_min_func, count / area, *args, **kwargs).x
        else:
            raise ValueError(
                f"{prior} is not a valid prior for global uniform poisson mapping."
            )

        return _map_value * np.ones(self.NPIX)

    def _bpm_lu(self, point_table, object_type, *args, **kwargs):
        mainlog.debug("Using LOCAL_UNIFORM to determine poisson map.")
        # -- Manage arguments and kwargs -- #
        prior = kwargs.pop("prior", None)
        map = np.zeros(self.NPIX)
        index_array = [i for i in np.arange(self.NPIX) if i in point_table["PIX_ID"]]

        # ------------------------------------------------------------------#
        # No PRIOR
        # ------------------------------------------------------------------#
        if prior is None:
            mainlog.debug("No prior detected: MLE estimator (analytic)")
            with logging_redirect_tqdm(loggers=[mainlog]):
                for ind in tqdm(index_array):
                    _pt = point_table[point_table["PIX_ID"] == ind]
                    count, area = np.sum(_pt[object_type]), np.sum(_pt["AREA"])
                    map[ind] = count / area
        elif callable(prior):
            from scipy.optimize import minimize
            from scipy.stats import poisson

            mainlog.debug("Prior detected: Using a MAP scheme.")
            with logging_redirect_tqdm(loggers=[mainlog]):
                for ind in tqdm(index_array):
                    # the prior is callable
                    _pt = point_table[point_table["PIX_ID"] == ind]
                    count, area = np.sum(_pt[object_type]), np.sum(_pt["AREA"])
                    _min_func = lambda x, c=count, a=area: -poisson.pmf(
                        c, a * x
                    ) * prior(x)
                    map[ind] = minimize(_min_func, count / area, *args, **kwargs).x
        else:
            raise ValueError(
                f"{prior} is not a valid prior for global uniform poisson mapping."
            )

        return map

    def _bpm_kern(self, point_table, *args, **kwargs):
        pass

    def append_to_fits(self, table, hudl):
        _self_hudl = fits.table_to_hdu(table)

        with fits.open(self.path, "update") as hudl_list:
            if hudl in hudl_list:
                _hudl, _len_hudl = hudl_list[hudl], len(hudl_list[hudl].data)
                new_hudl = fits.BinTableHDU.from_columns(
                    _hudl.columns, nrows=_len_hudl + len(table)
                )
                for colname in hudl_list[hudl].columns.names:
                    new_hudl.data[colname][_len_hudl:] = _self_hudl.data[colname]

                del hudl_list[hudl]
            else:
                new_hudl = fits.table_to_hdu(table)

            new_hudl.name = hudl
            hudl_list.append(new_hudl)
            hudl_list.flush()

    def plot_points_2d(self, object_type=None, *args, **kwargs):
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm, Normalize

        from pyXMIP.utilities.plot import scatter_table_2d

        # -- Pull the necessary points for the plotting procedure -- #
        points = self.get_points
        points["RA"] *= units.deg
        points["DEC"] *= units.deg

        # -- Setting up coloring -- #
        if object_type is not None:
            assert (
                object_type in points.columns
            ), f"The object type {object_type} is not in the corresponding table."

            color_scale = {"linear": Normalize, "log": LogNorm}[
                kwargs.pop("color_scale", "linear")
            ]
            ccmin, ccmax = np.amin(points[object_type]), np.amax(points[object_type])
            color_norm = color_scale(vmin=ccmin, vmax=ccmax)
            kwargs["c"] = color_norm(points[object_type])
            color_mappable = plt.cm.ScalarMappable(
                norm=color_norm, cmap=kwargs.get("cmap", None)
            )
        else:
            kwargs["c"] = "k"
            color_mappable = None

        fig, ax, _ = scatter_table_2d(points, *args, **kwargs)

        return fig, ax, color_mappable

    def plot_points_3d(self, object_type=None, *args, **kwargs):
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm, Normalize

        from pyXMIP.utilities.plot import scatter_table_3d

        # -- Pull the necessary points for the plotting procedure -- #
        points = self.get_points()
        points["RA"] *= units.deg
        points["DEC"] *= units.deg

        # -- Setting up coloring -- #
        if object_type is not None:
            assert (
                object_type in points.columns
            ), f"The object type {object_type} is not in the corresponding table."

            color_scale = {"linear": Normalize, "log": LogNorm}[
                kwargs.pop("color_scale", "linear")
            ]
            ccmin, ccmax = np.amin(points[object_type]), np.amax(points[object_type])
            color_norm = color_scale(vmin=ccmin, vmax=ccmax)
            kwargs["c"] = color_norm(points[object_type])
            color_mappable = plt.cm.ScalarMappable(
                norm=color_norm, cmap=kwargs.get("cmap", None)
            )
        else:
            kwargs["c"] = "k"
            color_mappable = None

        fig, ax, _ = scatter_table_3d(points, *args, **kwargs)

        return fig, ax, color_mappable


class Map:
    """
    Representation of a callable HEALPix map in an Atlas.
    """

    def __init__(self, path, name):
        """
        Initialize the specified Map HDU from the atlas located at ``path``.

        Parameters
        ----------
        path: str
            The path to the :py:class:`MapAtlas` instance.
        name: str
            The name of the HDU to load.
        """
        mainlog.debug(f"Loading Map object from {path} [Name={name}].")

        # -- loading basic attributes -- #
        self.path = pt.Path(path)
        self.name = name.upper()

        # -- reading the values -- #
        with fits.open(self.path) as hudl:
            assert name in [
                u.name for u in hudl
            ], f"There is no SkyMap {name} in {path}."

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
        lon, lat = (
            position.transform_to(self.coordinate_frame).frame.spherical.lon.rad,
            position.transform_to(self.coordinate_frame).frame.spherical.lat.rad,
        )
        phi, theta = _spherical_coordinates_to_healpix(lon, lat)

        pixels = hp.ang2pix(self._nside, theta, phi)

        return self.data[pixels]

    @property
    def coordinate_frame(self):
        """The :py:class:`astropy.coordinates.GenericFrame` for the Atlas's coordinate system"""
        with fits.open(self.path) as hudl:
            return getattr(astro_coords, hudl[0].header["CSYS"])


def _healpix_coordinates_to_spherical(phi, theta):
    return phi, np.pi / 2 - theta


def _spherical_coordinates_to_healpix(phi, theta):
    return _healpix_coordinates_to_spherical(
        phi, theta
    )  # --> this is a self-inverse transformation.
