r"""
Module for representing functions on the sky.

Notes
-----

The :py:mod:`structures.map` module provides an interface for managing data on the sky. The :py:class:`structures.map.MapAtlas` class
represents a collection of :py:class:`structures.map.Map` objects. Each of these maps acts like a function of sky position.
"""
import pathlib as pt
import warnings
from time import asctime
from types import SimpleNamespace

import healpy as hp
import numpy as np
from _collections_abc import Collection
from astropy import coordinates as astro_coords
from astropy import units
from astropy.io import fits
from tqdm.auto import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from pyXMIP.utilities.core import enforce_units
from pyXMIP.utilities.geo import convert_coordinates, convert_skycoord
from pyXMIP.utilities.logging import mainlog
from pyXMIP.utilities.plot import _enforce_style, plot_healpix


class _AtlasHeaderParam:
    def __set_name__(self, owner, name):
        self._name = name

    def __get__(self, instance, owner):
        try:
            val = getattr(instance, f"_{self._name}")
            if val is not None:
                return val
            else:
                raise AttributeError
        except AttributeError:
            with fits.open(instance.path) as hudl:
                setattr(instance, f"_{self._name}", hudl[0].header[self._name])
                return hudl[0].header[self._name]

    def __set__(self, instance, value):
        with fits.open(instance.path, "update") as hudl:
            hudl[0].header[self._name] = value
            hudl.flush()


class _MapHeaderParam:
    def __set_name__(self, owner, name):
        self._name = name

    def __get__(self, instance, owner):
        try:
            val = getattr(instance, f"_{self._name}")
            if val is not None:
                return val
            else:
                raise AttributeError
        except AttributeError:
            with fits.open(instance.path) as hudl:
                setattr(
                    instance, f"_{self._name}", hudl[instance.name].header[self._name]
                )
                return hudl[instance.name].header[self._name]

    def __set__(self, instance, value):
        with fits.open(instance.path, "update") as hudl:
            hudl[instance.name].header[self._name] = value
            hudl.flush()


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
    NPIX = _AtlasHeaderParam()
    NSIDE = _AtlasHeaderParam()
    CSYS = _AtlasHeaderParam()
    CD = _AtlasHeaderParam()
    ED = _AtlasHeaderParam()

    def __init__(self, filepath):
        """
        Initialize the :py:class:`MapAtlas` from a specified ``.fits`` file.

        Parameters
        ----------
        filepath: str
            The path to the ``.fits`` file to read from.

        """
        self.path = filepath

        if not pt.Path(self.path).exists():
            raise FileNotFoundError(f"There is no Atlas file at {self.path}.")

    def get_map(self, name):
        """
        Obtain an instance of the map corresponding the the name specified from this Atlas.

        Parameters
        ----------
        name: str
            The name of the map.

        Returns
        -------
        :py:class:`Map`

        Examples
        --------

        Let's try loading the galaxies poisson map for the SIMBAD database.

        .. code-block:: python

            from pyXMIP.structures.databases import SIMBAD
            poisson_atlas = SIMBAD.get_default_poisson_atlas() # load the default poisson map.
            map = poisson_atlas.get_map('G')

        """
        return Map(self.path, name.upper())

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
    def has_maps(self):
        return True if len(self.map_names) else False

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
        _th, _ph = hp.pix2ang(self.NSIDE, np.arange(self.NPIX))
        _ph, _th = convert_coordinates(
            _ph, _th, from_system="healpix", to_system="latlon"
        )
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
        header["EDATE"] = asctime()
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

    def remove(self):
        mainlog.info(f"Removing {self.path}.")
        pt.Path(self.path).unlink()
        del self

    @classmethod
    def _get_descriptors(cls):
        # fetching descriptors.
        return [m for m, v in cls.__dict__.items() if isinstance(v, _AtlasHeaderParam)]

    def _update_attributes(self):
        # fetch the attached descriptors
        descriptors = self._get_descriptors()

        for descriptor in descriptors:
            setattr(self, f"_{descriptor._name}", None)  # resets everything.


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

    DBNAME = _AtlasHeaderParam()

    def __init__(self, filepath):
        super().__init__(filepath)

    def __len__(self):
        return len(self.COUNTS)

    @property
    def COUNTS(self):
        from astropy.table import Table

        with warnings.catch_warnings(), fits.open(self.path) as hudl:
            warnings.simplefilter("ignore")
            return Table(hudl["COUNTS"].data)

    def get_points(self):
        """
        Fetch the underlying ``COUNTS`` table for the atlas.

        Returns
        -------
        :py:class:`astropy.table.Table`
        """
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
            count_positions.transform_to(self.coordinate_frame).spherical.lon.rad,
            count_positions.transform_to(self.coordinate_frame).spherical.lat.rad,
        )
        _p, _t = convert_coordinates(_p, _t, from_system="latlon", to_system="healpix")

        _out["PIX_ID"] = hp.ang2pix(self.NSIDE, _t, _p)

        return _out

    def sample_from_database(self, npoints, search_radius, *args, **kwargs):
        """
        Randomly sample points from the database specified in the meta-data of this atlas and add the count
        data to the pre-existing COUNTS table.

        Parameters
        ----------
        npoints: int
            The number of random sample points to query for.
        search_radius: :py:class:`astropy.units.Quantity`
            Angular search radius for counting around each point.
        args
            Additional arguments to pass through.
        kwargs
            Additional key-word arguments to pass through.

        Returns
        -------
        None
        """
        from pyXMIP.structures.databases import DEFAULT_DATABASE_REGISTRY

        mainlog.info(f"Adding {npoints} counting queries to {self.path} atlas.")

        # -- fetching the correct database -- #
        _registry = kwargs.pop("registry", DEFAULT_DATABASE_REGISTRY)
        if self.DBNAME == "NONE":
            raise KeyError(
                "This PoissonAtlas doesn't appear to be linked to a database. It may be corrupted or it was never set. Use instance.DBNAME = 'new_name' to set a database name."
            )
        elif self.DBNAME not in _registry:
            raise KeyError(
                "This Database doesn't correspond to a known database in the registry provided."
            )
        else:
            database = _registry[self.DBNAME]
        mainlog.debug(f"LINKED DATABASE: {database}")

        # -- Requesting pull from database -- #
        sample_count_table = database.random_sample_count(
            npoints, search_radius, *args, **kwargs
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.append_to_fits(sample_count_table, "COUNTS")

        mainlog.info(f"Appended {npoints} samples to COUNTS.")

    def reset(self, resolution, **kwargs):
        mainlog.info(f"Resetting Atlas at {self.path}.")
        path = self.path
        kwargs["database"] = kwargs.get("database", self.DBNAME)
        self.remove()
        self.generate(path, resolution, **kwargs)

    def _write_build_output_to_fits(self, object_type, output_object, overwrite=False):
        # =========================================== #
        # Writing the map to an HDU of the correct name
        object_type = object_type.upper()
        with fits.open(self.path, "update") as hdul:
            if object_type in [hdu.name for hdu in hdul] and not overwrite:
                # The map already exists and we cannot overwrite.
                raise ValueError(
                    f"The map {object_type} already exists in {self.path} and overwrite = False."
                )
            elif object_type in [hdu.name for hdu in hdul]:
                # The map exists we need to remove it.
                del hdul[object_type]

            # -- Writing the new map -- #
            image_hdu = fits.ImageHDU(output_object.map)
            image_hdu.name = object_type
            image_hdu.header["ISMAP"] = True
            image_hdu.header["METH"] = output_object.method
            image_hdu.header["DATE"] = asctime()
            image_hdu.header["CSYS"] = self.CSYS

            hdul.append(image_hdu)
            hdul.flush()

    @classmethod
    def generate(cls, path, resolution, overwrite=False, database="NONE"):
        obj = super().generate(path, resolution, overwrite=overwrite)
        obj.DBNAME = database

        return obj

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


class PoissonAtlas(StatAtlas):
    """
    Special case class of :py:class:`StatAtlas` for Poisson maps.
    """

    def build_poisson_map(
        self, object_type, method="KNN", overwrite=True, inplace=True, *args, **kwargs
    ):
        r"""
        Build a Poisson map in the Atlas for a specific object type.

        Parameters
        ----------
        object_type: str
            A particular object type in this PoissonAtlas from which to generate the map.
        method: str
            The method by which to calculate the map.
        overwrite: bool, optional
            [Default ``True``]. Allow the overwriting of pre-existing maps of the same name?
        inplace: bool, optional
            [Default ``True``]. If ``True``, the map is written directly to file. If ``False`` it is returned to the user.
        *args
            Additional arguments.
        **kwargs
            Additional keyword arguments.

        Returns
        -------

        """
        # ======================================================= #
        # SETUP
        # ======================================================= #
        _method = getattr(self, f"build_poisson_map_{method}")
        mainlog.info(f"Generating {object_type} map using {method}.")

        # ======================================================= #
        # RUNNING
        # ======================================================= #
        _output = _method(object_type, *args, **kwargs)

        # ===================================================== #
        # Writing data
        # ===================================================== #
        if inplace:
            self._write_build_output_to_fits(object_type, _output, overwrite=overwrite)
        else:
            return _output

    def build_poisson_maps(
        self,
        object_types="all",
        methods="MAP",
        overwrite=True,
        build_args=None,
        build_kwargs=None,
        **kwargs,
    ):
        r"""
        Construct all or a subset of the available object type Poisson maps for this Atlas.

        Parameters
        ----------
        object_types: list of str or str, optional
            The object types to generate Poisson maps for. If ``all`` (default), then all of the object types in the ``COUNTS``
            table are used.
        methods: list of str or str, optional
            The method to use for the map generation. If this is a list, it must match the length of ``object_types`` and each
            of the Poisson maps will be generated using the specified method. If ``methods`` is a single ``str``, then a single
            method will be used for all of the generation procedures.
        build_args: list of list, optional
            List of additional arguments to pass to each of the generation methods individually. By default, this is ``None``, and will result
            in empty lists being appended to the arguments of each generation method. If specified, this parameter must be a list of length ``object_types.size``
            with each element also being a list (or arbitrary size) containing any additional args to pass.
        build_kwargs: list of dict, optional
            Similar to ``build_args``, ``build_kwargs`` may be specified to pass particular keyword arguments along to the individual generation methods.
            If left unspecified, no kwargs are passed to the sub-generation processes.
        overwrite: bool, optional
            [Default ``True``]. Allow the overwriting of pre-existing maps of the same name?
        **kwargs
            Additional keyword arguments.

        Returns
        -------

        """
        # ======================================================================== #
        # Setup
        # ======================================================================== #
        _count_data = self.get_points()  # fetch the point data to process.

        # -- managing args / kwargs -- #
        if object_types == "all":
            object_types = list(
                _count_data.columns[:-5]
            )  # Removing RA,DEC,RAD,PIX_ID,TIME
        else:
            object_types = [
                object_type
                for object_type in object_types
                if object_type in _count_data.columns
            ]

        mainlog.info(f"Generating Poisson maps for {len(object_types)} object types...")

        # -- Managing computation methods -- #
        if isinstance(methods, str):
            methods = [methods] * len(object_types)
        elif isinstance(methods, Collection):
            assert len(methods) == len(
                object_types
            ), f"Attempted to pass {len(methods)} methods for {len(object_types)} object types."
        else:
            raise TypeError(
                f"Kwarg `methods` has type {type(methods)} when it was expected to be `str` or `Collection`."
            )

        # Managing the build args / kwargs
        if build_args is None:
            build_args = [[] for _ in range(len(object_types))]
        else:
            assert isinstance(
                build_args, Collection
            ), f"Kwarg `build_args` has type {type(build_args)} which is not a subclass of `Collection`."

        if build_kwargs is None:
            build_kwargs = [{} for _ in range(len(object_types))]
        else:
            assert isinstance(
                build_args, Collection
            ), f"Kwarg `build_kwargs` has type {type(build_kwargs)} which is not a subclass of `Collection`."

        # ======================================================================== #
        # Map Generation Processes
        # ======================================================================== #
        with logging_redirect_tqdm(loggers=[mainlog]):
            for object_type, method, bargs, bkwargs in tqdm(
                zip(object_types, methods, build_args, build_kwargs),
                total=len(object_types),
                desc="Constructing Poisson Maps",
                disable=(not kwargs.get("progress_bar", True)),
                leave=True,
            ):
                bkwargs["progress_bar"] = kwargs.get(
                    "progress_bar", bkwargs.get("progress_bar", True)
                )

                self.build_poisson_map(
                    object_type,
                    method=method,
                    overwrite=overwrite,
                    inplace=True,
                    *bargs,
                    **bkwargs,
                )

    def build_poisson_map_regressor(
        self,
        regressor,
        object_type,
        retrain=True,
        cross_validate=True,
        training_kw=None,
        param_kw=None,
        cv_kw=None,
        **kwargs,
    ):
        """
        Build a poisson map with a regressor.
        Parameters
        ----------
        regressor
        retrain
        cross_validate
        training_kw
        param_kw

        Returns
        -------

        """
        # =================================================== #
        # SETUP
        # =================================================== #
        phi, theta = convert_skycoord(
            self.pixel_positions.transform_to("icrs"), "math_convention"
        )
        map_positions = np.vstack(
            [theta, phi]
        ).T  # --> This is done in oposite order for haversine.

        # get the points that are being fitted.
        count_data = self.get_points()

        if cv_kw is None:
            cv_kw = {}
        # ================================================= #
        # Cross Validation
        # ================================================= #
        if cross_validate:
            cross_validation_data = regressor.cross_validate(
                count_data,
                object_type,
                training_kw=training_kw,
                param_kw=param_kw,
                **cv_kw,
            )

            # Parse the data and find the best score
            test_scores = cross_validation_data["mean_test_score"]
            best_score, best_params = (
                np.amax(test_scores),
                cross_validation_data["params"][
                    list(cross_validation_data["rank_test_score"]).index(1)
                ],
            )

            if kwargs.get("cv_score_threshold", None) is not None:
                # the user has specified a cross validation score threshold.
                _threshold = kwargs.pop("cv_score_threshold")
                assert (
                    best_score >= _threshold
                ), f"Score threshold was {_threshold} and best CV score was {best_score}."

            mainlog.debug(f"CROSS-VALIDATED {object_type}. PARAMS:")
            for k, v in best_params.items():
                mainlog.debug(f"{k} = {v}")

            regressor.regressor.set_params(**best_params)
        else:
            if param_kw is not None:
                regressor.regressor.set_params(**{k: v[0] for k, v in param_kw.items()})
        # ================================================= #
        # Training
        # ================================================= #
        if retrain:
            training_kw = {} if training_kw is None else training_kw
            score = regressor.train_model(count_data, object_type, **training_kw)

            if kwargs.get("score_threshold", None) is not None:
                # the user has specified a cross validation score threshold.
                _threshold = kwargs.pop("score_threshold")
                assert (
                    score >= _threshold
                ), f"Score threshold was {_threshold} and best test score was {score}."

        # ================================================= #
        # Output
        # ================================================= #
        _out_map = regressor.regressor.predict(map_positions)

        return SimpleNamespace(map=_out_map, method=regressor.__class__.__name__)

    def build_poisson_map_KNN(self, object_type, *args, **kwargs):
        from pyXMIP.stats.map_regression import KNNeighborMapRegressor

        _ = args

        return self.build_poisson_map_regressor(
            KNNeighborMapRegressor(), object_type, **kwargs
        )

    def build_poisson_map_RNN(self, object_type, *args, **kwargs):
        from pyXMIP.stats.map_regression import RNNeighborMapRegressor

        _ = args

        return self.build_poisson_map_regressor(
            RNNeighborMapRegressor(), object_type, **kwargs
        )

    def build_poisson_map_MAP(self, object_type, *args, **kwargs):
        from pyXMIP.stats.map_regression import BayesianPoissonMapRegressor

        _ = args
        # -- pull counts data -- #
        count_table = self.get_points()

        # -- setup the regressor -- #
        regressor = BayesianPoissonMapRegressor(
            prior=kwargs.pop("prior", None), model=kwargs.pop("model", None)
        )

        estimates = regressor.build_map_MAP(count_table, object_type, **kwargs)

        # -- returning the map -- #
        map = np.zeros(self.NPIX)
        map[: len(estimates)] = estimates

        return SimpleNamespace(map=map, method=regressor.__class__.__name__)


class Map:
    """
    Representation of a callable HEALPix map in an Atlas.
    """

    NPIX = _AtlasHeaderParam()
    NSIDE = _AtlasHeaderParam()
    CSYS = _MapHeaderParam()
    CD = _MapHeaderParam()
    ED = _MapHeaderParam()

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

            self._data = hudl[name].data

    def __call__(self, position):
        """
        Evaluate the skymap.

        Parameters
        ----------
        position: :py:class:`astropy.coordinates.SkyCoord`

        """
        pixels = self.get_healpix_id(position)

        return self.data[pixels]

    @classmethod
    def _get_descriptors(cls):
        # fetching descriptors.
        return [
            m
            for m, v in cls.__dict__.items()
            if isinstance(v, (_AtlasHeaderParam, _MapHeaderParam))
        ]

    def _update_attributes(self):
        # fetch the attached descriptors
        descriptors = self._get_descriptors()

        for descriptor in descriptors:
            setattr(self, f"_{descriptor}", None)  # resets everything.

    @property
    def data(self):
        try:
            val = self._data
            if val is not None:
                return val
            else:
                raise AttributeError
        except AttributeError:
            with fits.open(self.path) as hudl:
                self._data = hudl[self.name].data
                return self._data

    @data.setter
    def data(self, value):
        with fits.open(self.path, "update") as hudl:
            hudl[self.name].data = value
            hudl.flush()

        self.ED = asctime()

    @property
    def coordinate_frame(self):
        """The :py:class:`astropy.coordinates.GenericFrame` for the Atlas's coordinate system"""
        with fits.open(self.path) as hudl:
            return getattr(astro_coords, hudl[self.name].header["CSYS"])

    def transform_map_coordinates(self, frame, inplace=False):
        # -- coercing the frame type -- #
        if isinstance(frame, str):
            frame = getattr(astro_coords, frame)
        else:
            pass

        # -- converting -- #
        _th, _ph = hp.pix2ang(self.NSIDE, np.arange(self.NPIX))
        _ph, _th = convert_coordinates(
            _ph, _th, from_system="healpix", to_system="latlon"
        )
        coords = astro_coords.SkyCoord(_ph, _th, frame=frame, unit="rad")
        output_array = self(coords)

        if inplace:
            with fits.open(self.path, "update") as hudl:
                hudl[self.name].header["CSYS"] = frame.__name__
                hudl[self.name].data = output_array
                hudl.flush()
            self.ED = asctime()
            self._update_attributes()
            self._data = output_array

        return output_array

    @property
    def pixel_positions(self):
        """The SkyCoord positions of the healpix pixels."""
        _th, _ph = hp.pix2ang(self.NSIDE, np.arange(self.NPIX))
        _ph, _th = convert_coordinates(
            _ph, _th, from_system="healpix", to_system="latlon"
        )
        return astro_coords.SkyCoord(_ph, _th, frame=self.coordinate_frame, unit="rad")

    def get_healpix_id(self, positions, frame=None):
        if frame is None:
            frame = self.coordinate_frame

        positions = positions.transform_to(frame)
        _ph, _th = convert_skycoord(positions, "healpix")
        print(_ph, _th)
        ids = hp.ang2pix(self.NSIDE, _th, _ph)

        return ids

    @_enforce_style
    def plot(self, *args, **kwargs):
        return plot_healpix(self.data, *args, **kwargs)


def _parse_default_kwarg_groups(defaults, kwargs):
    # -- setting args and kwargs -- #
    out = {k: kwargs.pop(k, None) for k in defaults}

    for def_k, def_v in defaults.items():
        out[def_k] = {
            k: (v if (out[def_k] is None or k not in out[def_k]) else out[def_k][k])
            for k, v in def_v.items()
        }

    return tuple([v for v in out.values()])
