"""
Skymapping handlers and classes for pyXMIP.
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
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from pyXMIP.utilities.core import enforce_units, mainlog


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

        """
        return Map(self.path, name)

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

    def reset(self, *args, **kwargs):
        mainlog.info(f"Resetting Atlas at {self.path}.")
        path = self.path
        self.remove()
        self.generate(path, *args, **kwargs)

    def _write_build_output_to_fits(self, object_type, output_object, overwrite=False):
        # =========================================== #
        # Writing the map to an HDU of the correct name
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

            hdul.append(image_hdu)
            hdul.flush()

    def build_poisson_map(
        self, object_type, method="GP", overwrite=True, inplace=True, *args, **kwargs
    ):
        """
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

        _count_data = self.get_points()
        _method = getattr(self, f"build_poisson_map_{method}")
        mainlog.info(f"Generating {object_type} map using {method}.")
        # ======================================================= #
        # Managing different methods
        # ======================================================= #
        if method in ["MAP"]:
            # These methods take object type, *args, **kwargs as their arguments.
            #
            #
            _output = _method(object_type, args, **kwargs)
        else:
            # These methods take X,Y, *args, **kwargs
            #
            #
            positions = astro_coords.SkyCoord(
                ra=_count_data["RA"], dec=_count_data["DEC"], unit="deg"
            )
            X = np.vstack(
                [positions.frame.spherical.lon.rad, positions.frame.spherical.lat.rad]
            )
            _r = _count_data["RAD"] * units.arcmin
            Y = _count_data[object_type] / (np.pi * (_r.to_value("rad")) ** 2)

            _output = _method(X, Y, *args, **kwargs)

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
        """
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
        _count_data = self.get_points()

        if object_types == "all":
            # the object types need to be grabbed. Removing RA,DEC,RAD,PIX_ID,TIME
            object_types = list(_count_data.columns[:-5])
        else:
            object_types = [
                object_type
                for object_type in object_types
                if object_type in _count_data.columns
            ]

        mainlog.info(f"Generating Poisson maps for {len(object_types)} object types...")

        # Managing the methods
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
            ), f"Kwarg `build_kwargs` has type {type(build_args)} which is not a subclass of `Collection`."

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

    def build_poisson_map_gp(
        self, X, Y, training_kw=None, parallel_kw=None, *args, **kwargs
    ):
        r"""
        Build the poisson map using a Haversine based Gaussian Process regression methodology.

        Parameters
        ----------
        X: array-like
            An array of size ``N,2`` containing the coordinates :math:`\phi,\theta` (in that order).

            .. warning::

                The coordinates MUST be in valid LAT/LON format: :math:`-\pi<\phi<\pi` and :math:`-\pi/2<\theta<\pi/2`.
                If this is not done, then the Gaussian Process regression will fail to correctly compute distances between
                points on the sphere and the resulting map may be entirely erroneous.

        Y: array-like
            The density values (``N,1``) at each of the provided coordinate positions.

            .. hint::

                Behind the scene, this is converted to the log-density to avoid issues with the GP regressor allowing
                negative densities. It is then rectified after the GP has been created.
        training_kw: dict
            Training key-word arguments specific to this ML algorithm.

            +-----------------------------------+-------------------------------------------------+--------------------+------------------------+
            | Training Key Word Argument        | Description                                     | Expected Type      | Default                |
            +===================================+=================================================+====================+========================+
            | ``training_proportion``           | The proportion of the provided data to use as   | ``float``          | ``0.2``                |
            |                                   | the training set.                               |                    |                        |
            +-----------------------------------+-------------------------------------------------+--------------------+------------------------+
            | ``cv_groups``                     | The number of groups to split the training      | ``int``            | ``5``                  |
            |                                   | set into for cross-validation.                  |                    |                        |
            +-----------------------------------+-------------------------------------------------+---------------------------------------------+
            | ``length_scale``                  | The initial guess for the Haversine Kernel's    | ``float``          | ``1.0``                |
            |                                   | characteristic length scale.                    |                    |                        |
            +-----------------------------------+-------------------------------------------------+---------------------------------------------+
            | ``length_scale_bounds``           | Bounds on the characteristic length scale.      | ``tuple``          | ``(1e-2,1e2)``         |
            |                                   |                                                 |                    |                        |
            +-----------------------------------+-------------------------------------------------+---------------------------------------------+
            | ``alpha_cv``                      | The :math:`\alpha` values for the GP regressor  | ``list``           | Log spaced from        |
            |                                   | to be compared during cross-validation.         | [if not bounds,    | ``1e-8`` to ``1e0``    |
            |                                   | If a 2-``tuple``, then linearly spaced between  | then must match    |                        |
            |                                   | bounds. Else, individual values.                | ``cv_groups`` ]    |                        |
            +-----------------------------------+-------------------------------------------------+---------------------------------------------+

        parallel_kw: dict
            Additional kwargs to pass for parallelism capabilities.

            +-----------------------------------+-------------------------------------------------+--------------------+------------------------+
            | Training Key Word Argument        | Description                                     | Expected Type      | Default                |
            +===================================+=================================================+====================+========================+
            | ``multiprocess``                  | If ``True``, then multi-core processing is      | ``bool``           | ``False``              |
            |                                   | enabled. Otherwise single-core.                 |                    |                        |
            +-----------------------------------+-------------------------------------------------+--------------------+------------------------+
            | ``nproc``                         | The number of processors to utilize.            | ``int``            | ``n_cpu``              |
            |                                   | [Cannot exceed ``cv_groups``]                   |                    |                        |
            +-----------------------------------+-------------------------------------------------+---------------------------------------------+

        Returns
        -------

        """
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.model_selection import train_test_split

        from pyXMIP.stats._gaussian_process import Haversine_RBF

        # -- Defaults -- #
        _default_training_kw = {
            "training_proportion": 0.2,
            "cv_groups": 20,
            "length_scale": 1,
            "length_scale_bounds": (1e-5, 1e5),
            "alpha_cv": np.geomspace(1e-1, 1e3, 20),
        }
        _default_parallel_kw = {"multiprocess": False, "nproc": None}

        mainlog.debug("Generating map from GPR (HVSINE-RBF).")
        with logging_redirect_tqdm(loggers=[mainlog]):
            _process_progress_bar = tqdm(
                total=100,
                desc="GP HVSINE-RBF - Setup",
                bar_format="{desc}: {percentage:3.0f}%|{bar}|",
                disable=(not kwargs.pop("progress_bar", True)),
                leave=False,
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # ================================================== #
                # Arguments / kwargs
                training_kw = {
                    k: (
                        v
                        if (training_kw is None or k not in training_kw)
                        else training_kw[k]
                    )
                    for k, v in _default_training_kw.items()
                }
                parallel_kw = {
                    k: (
                        v
                        if (parallel_kw is None or k not in parallel_kw)
                        else parallel_kw[k]
                    )
                    for k, v in _default_parallel_kw.items()
                }

                # X should be Nx2, Y should be Nx1 in shape.
                X, Y = X.reshape((X.size // 2, 2)), Y.reshape(
                    (Y.size, 1)
                )  # This will yield an error if the arrays are not correct.

                # ================================================== #
                # splitting training data and testing data.
                X_TRAIN, X_TEST, Y_TRAIN, Y_TEST = train_test_split(
                    X,
                    Y,
                    train_size=training_kw["training_proportion"],
                    shuffle=kwargs.pop("shuffle_features", False),
                    random_state=kwargs.pop("random_state", None),
                )

                # ================================================== #
                # Generating the Kernel and Cross-validating.
                kernel = Haversine_RBF(
                    length_scale=training_kw["length_scale"],
                    length_scale_bounds=training_kw["length_scale_bounds"],
                )

                _process_progress_bar.update(10)
                _process_progress_bar.desc = "GP HVSINE-RBF - Cross Validating"
                # If there are no CV groups, we just fit the regressor and proceed.
                if training_kw["cv_groups"] in [0, 1]:
                    # cross-validation doesn't occur. We simply proceed to fit.
                    _gaussian_process_regressor = GaussianProcessRegressor(
                        kernel=kernel, *args, **kwargs
                    )
                    _cv_return_data = None
                else:
                    # cross-validation does occur
                    from sklearn.model_selection import GridSearchCV

                    _ = kwargs.pop(
                        "alpha", None
                    )  # Remove alpha if it is present to allow for X-val.
                    _gaussian_process_regressor = GaussianProcessRegressor(
                        kernel=kernel, *args, **kwargs
                    )

                    if (
                        len(training_kw["alpha_cv"]) == 2
                        and training_kw["cv_groups"] != 2
                    ):
                        # we have a query region.
                        _gscv = GridSearchCV(
                            _gaussian_process_regressor,
                            {
                                "alpha": np.geomspace(
                                    *training_kw["alpha_cv"], training_kw["cv_groups"]
                                )
                            },
                        )
                    else:
                        _gscv = GridSearchCV(
                            _gaussian_process_regressor,
                            {"alpha": training_kw["alpha_cv"]},
                        )

                    # pulling scores
                    _gscv.fit(X_TRAIN, Y_TRAIN)
                    _cv_return_data = _gscv.cv_results_
                    _alpha_fixed = _cv_return_data["param_alpha"].data[
                        np.where(
                            _cv_return_data["mean_test_score"]
                            == np.amax(_cv_return_data["mean_test_score"])
                        )
                    ]

                    _process_progress_bar.update(50)
                    _process_progress_bar.desc = (
                        "GP HVSINE-RBF - Fitting to training set"
                    )
                    mainlog.debug(
                        f"CV determined alpha to be {_alpha_fixed} with mean test score {np.amax(_cv_return_data['mean_test_score'])}."
                    )

                    # creating the regressor
                    kwargs["alpha"] = _alpha_fixed
                    _gaussian_process_regressor = GaussianProcessRegressor(
                        kernel=kernel, *args, **kwargs
                    )

                # ===============================================================
                # Training fully
                _gaussian_process_regressor.fit(X_TRAIN, Y_TRAIN)

                # ===============================================================
                # Getting scores.
                _return_score = _gaussian_process_regressor.score(X_TEST, Y_TEST)
                _process_progress_bar.update(20)
                mainlog.debug(f"Fit to training set yielded score of {_return_score}.")
                _process_progress_bar.desc = "GP HVSINE-RBF - Predicting Poisson Map"
                # =================================================================
                # Generating the actual sky map.
                map_positions = np.vstack(
                    [
                        self.pixel_positions.frame.spherical.lat.rad,
                        self.pixel_positions.frame.spherical.lon.rad,
                    ]
                ).T

                map_values, map_std = _gaussian_process_regressor.predict(
                    map_positions, return_std=True
                )

                _process_progress_bar.update(20)
                _process_progress_bar.close()
        mainlog.debug("GP HVSINE-RBF - [COMPLETE]")
        return SimpleNamespace(
            map=map_values,
            std_map=map_std,
            training_info={"score": _return_score, "cv_info": _cv_return_data},
            method="GP",
            parameters={"alpha": _alpha_fixed},
        )

    def build_poisson_map_NNR(
        self, X, Y, training_kw=None, parallel_kw=None, *args, **kwargs
    ):
        r"""
        Build the poisson map using a Haversine based KNN radius matching approach.

        Parameters
        ----------
        X: array-like
            An array of size ``N,2`` containing the coordinates :math:`\phi,\theta` (in that order).

            .. warning::

                The coordinates MUST be in valid LAT/LON format: :math:`-\pi<\phi<\pi` and :math:`-\pi/2<\theta<\pi/2`.
                If this is not done, then the Gaussian Process regression will fail to correctly compute distances between
                points on the sphere and the resulting map may be entirely erroneous.

        Y: array-like
            The density values (``N,1``) at each of the provided coordinate positions.

            .. hint::

                Behind the scene, this is converted to the log-density to avoid issues with the GP regressor allowing
                negative densities. It is then rectified after the GP has been created.
        training_kw: dict
            Training key-word arguments specific to this ML algorithm.

            +-----------------------------------+-------------------------------------------------+--------------------+------------------------+
            | Training Key Word Argument        | Description                                     | Expected Type      | Default                |
            +===================================+=================================================+====================+========================+
            | ``training_proportion``           | The proportion of the provided data to use as   | ``float``          | ``0.2``                |
            |                                   | the training set.                               |                    |                        |
            +-----------------------------------+-------------------------------------------------+--------------------+------------------------+
            | ``cv_groups``                     | The number of groups to split the training      | ``int``            | ``5``                  |
            |                                   | set into for cross-validation.                  |                    |                        |
            +-----------------------------------+-------------------------------------------------+---------------------------------------------+
            | ``radius``                        | The initial guess for the KNN search radius.    | ``float``          | ``1.0``                |
            |                                   |                                                 |                    |                        |
            +-----------------------------------+-------------------------------------------------+---------------------------------------------+
            | ``radius_bounds``                 | Bounds on the characteristic KNNR scale (deg)   | ``tuple``          | ``(0.5,10)``           |
            |                                   |                                                 |                    |                        |
            +-----------------------------------+-------------------------------------------------+---------------------------------------------+

        parallel_kw: dict
            Additional kwargs to pass for parallelism capabilities.

            +-----------------------------------+-------------------------------------------------+--------------------+------------------------+
            | Training Key Word Argument        | Description                                     | Expected Type      | Default                |
            +===================================+=================================================+====================+========================+
            | ``multiprocess``                  | If ``True``, then multi-core processing is      | ``bool``           | ``False``              |
            |                                   | enabled. Otherwise single-core.                 |                    |                        |
            +-----------------------------------+-------------------------------------------------+--------------------+------------------------+
            | ``nproc``                         | The number of processors to utilize.            | ``int``            | ``n_cpu``              |
            |                                   | [Cannot exceed ``cv_groups``]                   |                    |                        |
            +-----------------------------------+-------------------------------------------------+---------------------------------------------+

        Returns
        -------

        """
        from sklearn.model_selection import train_test_split
        from sklearn.neighbors import RadiusNeighborsRegressor

        # -- Defaults -- #
        _default_training_kw = {
            "training_proportion": 0.2,
            "cv_groups": 20,
            "radius": 1,
            "radius_bounds": tuple((np.array([0.1, 30]) * units.deg).to_value("rad")),
        }
        _default_parallel_kw = {"multiprocess": False, "nproc": None}

        mainlog.debug("Generating map from NNR (HVSINE).")
        with logging_redirect_tqdm(loggers=[mainlog]):
            _process_progress_bar = tqdm(
                total=100,
                desc="NNR HVSINE - Setup",
                bar_format="{desc}: {percentage:3.0f}%|{bar}|",
                disable=(not kwargs.pop("progress_bar", True)),
                leave=False,
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # ================================================== #
                # Arguments / kwargs
                training_kw = {
                    k: (
                        v
                        if (training_kw is None or k not in training_kw)
                        else training_kw[k]
                    )
                    for k, v in _default_training_kw.items()
                }
                parallel_kw = {
                    k: (
                        v
                        if (parallel_kw is None or k not in parallel_kw)
                        else parallel_kw[k]
                    )
                    for k, v in _default_parallel_kw.items()
                }

                # X should be Nx2, Y should be Nx1 in shape.
                X, Y = X.reshape((X.size // 2, 2)), Y.reshape(
                    (Y.size, 1)
                )  # This will yield an error if the arrays are not correct.

                # ================================================== #
                # splitting training data and testing data.
                X_TRAIN, X_TEST, Y_TRAIN, Y_TEST = train_test_split(
                    X,
                    Y,
                    train_size=training_kw["training_proportion"],
                    shuffle=kwargs.pop("shuffle_features", False),
                    random_state=kwargs.pop("random_state", None),
                )

                # ================================================== #
                # Generating the Kernel and Cross-validating.

                _process_progress_bar.update(10)
                _process_progress_bar.desc = "NNR HVSINE - Cross Validating"
                # If there are no CV groups, we just fit the regressor and proceed.
                if training_kw["cv_groups"] in [0, 1]:
                    # cross-validation doesn't occur. We simply proceed to fit.
                    _nnr_regressor = RadiusNeighborsRegressor(
                        radius=training_kw["radius"],
                        metric="haversine",
                        algorithm="brute",
                        weights=kwargs.get("weights", "distance"),
                        n_jobs=(
                            parallel_kw["nproc"]
                            if parallel_kw["multiprocess"]
                            else None
                        ),
                    )
                    _cv_return_data = None
                else:
                    # cross-validation does occur
                    from sklearn.model_selection import GridSearchCV

                    _nnr_regressor = RadiusNeighborsRegressor(
                        metric="haversine",
                        algorithm="brute",
                        weights=kwargs.get("weights", "distance"),
                        n_jobs=(
                            parallel_kw["nproc"]
                            if parallel_kw["multiprocess"]
                            else None
                        ),
                    )

                    # we have a query region.
                    _gscv = GridSearchCV(
                        _nnr_regressor,
                        {
                            "radius": np.geomspace(
                                *training_kw["radius_bounds"], training_kw["cv_groups"]
                            )
                        },
                    )

                    # pulling scores
                    _gscv.fit(X_TRAIN, Y_TRAIN)
                    _cv_return_data = _gscv.cv_results_
                    _cv_return_data["param_radius"] = _cv_return_data["param_radius"][
                        ~np.isnan(_cv_return_data["mean_test_score"])
                    ]
                    _cv_return_data["mean_test_score"] = _cv_return_data[
                        "mean_test_score"
                    ][~np.isnan(_cv_return_data["mean_test_score"])]
                    _radius_fixed = _cv_return_data["param_radius"].data[
                        np.where(
                            _cv_return_data["mean_test_score"]
                            == np.amax(_cv_return_data["mean_test_score"])
                        )
                    ]

                    _process_progress_bar.update(50)
                    _process_progress_bar.desc = "NNR HVSINE - Fitting to training set"
                    mainlog.debug(
                        f"CV determined radius to be {_radius_fixed} with mean test score {np.amax(_cv_return_data['mean_test_score'])}."
                    )

                    if isinstance(_radius_fixed, np.ndarray) and len(_radius_fixed):
                        _radius_fixed = _radius_fixed[0]
                    elif isinstance(_radius_fixed, np.ndarray) and not len(
                        _radius_fixed
                    ):
                        mainlog.error(
                            f"Failed to find a valid radius from CV. Using default {training_kw['radius']}."
                        )
                        _radius_fixed = training_kw["radius"]
                    else:
                        pass

                    # creating the regressor
                    _nnr_regressor = RadiusNeighborsRegressor(
                        radius=_radius_fixed,
                        metric="haversine",
                        algorithm="brute",
                        weights=kwargs.get("weights", "distance"),
                        n_jobs=(
                            parallel_kw["nproc"]
                            if parallel_kw["multiprocess"]
                            else None
                        ),
                    )

                # ===============================================================
                # Training fully
                _nnr_regressor.fit(X_TRAIN, Y_TRAIN)

                # ===============================================================
                # Getting scores.
                _return_score = _nnr_regressor.score(X_TEST, Y_TEST)
                _process_progress_bar.update(20)
                mainlog.debug(f"Fit to training set yielded score of {_return_score}.")
                _process_progress_bar.desc = "NNR HVSINE - Predicting Poisson Map"
                # =================================================================
                # Generating the actual sky map.
                map_positions = np.vstack(
                    [
                        self.pixel_positions.frame.spherical.lat.rad,
                        self.pixel_positions.frame.spherical.lon.rad,
                    ]
                ).T

                map_values = _nnr_regressor.predict(map_positions)
                map_values = map_values.reshape((map_values.size,))
                _process_progress_bar.update(20)
                _process_progress_bar.close()
        mainlog.debug("NNR HVSINE - [COMPLETE]")

        return SimpleNamespace(
            map=map_values,
            std_map=None,
            training_info={"score": _return_score, "cv_info": _cv_return_data},
            method="NNR",
            parameters={"radius": _radius_fixed},
        )

    def build_poisson_map_MAP(self, object_id, *args, **kwargs):
        mainlog.debug("Generating map from posterior (MAP).")
        with logging_redirect_tqdm(loggers=[mainlog]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")

                # ===============================================#
                # Constructing grids
                _pix_ids = np.arange(self.NPIX)
                _map_array = np.zeros(self.NPIX)
                _point_data = self.get_points()

                _non_empty_pix_ids = [
                    pid for pid in _pix_ids if pid in _point_data["PIX_ID"]
                ]

                for pix_id in tqdm(
                    _non_empty_pix_ids,
                    desc="posterior (MAP) - Calculating",
                    disable=(not kwargs.pop("progress_bar", True)),
                    leave=False,
                ):
                    _c = np.sum(
                        _point_data[object_id][
                            np.where(_point_data["PIX_ID"] == pix_id)
                        ]
                    )
                    _r = _point_data["RAD"][np.where(_point_data["PIX_ID"] == pix_id)]

                    # -- convert radii to rad -- #
                    _r = (np.array(_r) * units.arcmin).to_value("rad")
                    _a = np.pi * _r**2
                    _map_array[pix_id] = _c / np.sum(_a)

        return SimpleNamespace(
            map=_map_array,
            std_map=None,
            training_info=None,
            method="MAP",
            parameters=None,
        )

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
