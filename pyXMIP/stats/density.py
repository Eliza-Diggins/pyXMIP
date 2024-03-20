r"""
Statistics methods and classes for point density estimation and density mapping.

Overview
--------

Density mapping plays a critical role in the cross-matching process by allowing the statistical determination of our
confidence in the identification of a source based only on its type and distance. In principle, if a source at position
:math:`\textbf{r}` is matched to a catalog / database source at position :math:`\textbf{r}'` of type :math:`j` then we can
model the probability of the match being spurious if we know the density of such objects around the identification point.

To accomplish this, we construct position density maps for matching catalogs; the methods for which are contained in this
module.
"""
import operator
import random
from collections.abc import Collection
from numbers import Number
from string import ascii_letters
from types import SimpleNamespace

import astropy.units as units
import numpy as np
from astropy.coordinates import ICRS, SkyCoord
from tqdm.contrib.logging import tqdm_logging_redirect

from pyXMIP.utilities.core import _enforce_style, mainlog, xsparams


def _do_op(operation, function1, function2, narg):
    return lambda phi, theta, *args, f1=function1, f2=function2, n1=narg: operation(
        f1(phi, theta, *args[:n1]), f2(phi, theta, *args[n1:])
    )


class PopulationDensity:
    r"""
    Class representing the density of a particular population of sources on the sky.
    """

    def __init__(
        self, coordinates, counts, area=np.pi * units.arcmin**2, population_name=None
    ):
        r"""
        Initializes the :py:class:`stats.density.PopulationDensity` instance.

        Parameters
        ----------
        coordinates: SkyCoord or list of SkyCoord
            Sky coordinates of the sampling positions. Length :math:`N`.
        counts: array-like, list of array-like
            The counts in each area, must be of length :math:`N`.
        area: array-like or float or int, optional
            The area of the search region for each sample. If an array, separate values corresond to each of the sample points,
            if a single value is provided, then only one area is attributed to each point. By default, area is a circular area of
            radius 1 arcmin on the sky.
        population_name: str, optional
            The name of this population.
        """
        #: The coordinate array for the population.
        self.coordinates = coordinates
        #: the counts array for the population.
        self.counts = counts
        #: The sampling area for the population
        self.area = area
        #: The name of the population.
        self.name = (
            population_name
            if population_name is not None
            else "".join(random.choice(ascii_letters) for i in range(10))
        )
        # --------------------------------------------------------------------------------------------------------------#
        # Sanity Check for sizing
        # --------------------------------------------------------------------------------------------------------------#
        if not isinstance(self.coordinates, SkyCoord):
            raise ValueError(
                f"coordinates had type {type(self.coordinates)}, must be a collection of SkyCoords or single SkyCoord."
            )

        if (isinstance(self.area, Collection)) and (
            not isinstance(self.area, units.Quantity)
        ):
            self.area = [
                k if isinstance(k, units.Quantity) else k * units.arcmin**2
                for k in self.area
            ]
        elif not isinstance(self.area, units.Quantity):
            # the area is scalar
            mainlog.warning(
                "A scalar area was provided, units are defaulted to arcmin^2."
            )
            self.area = self.area * units.arcmin**2

        if self.area.isscalar:
            # we have units, the area is scalar.
            self.area = self.area * np.ones(self.coordinates.size)

        assert (
            len(self.counts) == len(self.area) == len(self.coordinates)
        ), f"Input arrays are not of the same length: {len(self.counts)},{len(self.area)},{len(self.coordinates)}"

        # printing some basic logging info.
        mainlog.info(f"Loading population of size {len(self)} with name {self.name}")

        # --------------------------------------------------------------------------------------------------------------#
        # Additional Attributes
        # --------------------------------------------------------------------------------------------------------------#
        self._model = None
        self._map_params = None

    def __len__(self):
        return len(self.counts)

    def __repr__(self):
        return f"<PopulationDensity {self.name} [{len(self)}]>"

    def __str__(self):
        return self.__repr__()

    @property
    def model(self):
        if self._model is None:
            mainlog.info(
                f"Population {self} doesn't have a model yet. Using a UniformDensityModel."
            )
            self._model = UniformDensityModel()

        return self._model

    @property
    def positions(self):
        return (
            self.coordinates.frame.spherical.lon.wrap_at(2 * np.pi * units.rad).rad,
            self.coordinates.frame.spherical.lat.wrap_at(units.rad * np.pi / 2).rad,
        )

    def fit_model(self, *args, **kwargs):
        self.model.fit_to_population(self, *args, **kwargs)

        return self.model

    def get_gaussian_kde(self, PHI, THETA):
        r"""
        Compute the gaussian kde on the domain :math:`0 \le \phi \le 2\pi,\;-\pi/2\le \theta <\pi/2`.

        .. warning::

            The gaussian KDE approach uses a gaussian kernel to smooth the density of points. It should be noted that,
            in most cases, this is not a particularly well constrained approach and may yield results with do not represent
            the heuristic intention of it's use.

        Parameters
        ----------
        PHI: array-like
            The grid of :math:`\phi` coordinates.
        THETA: array-like
            The grid of :math:`\theta` coordinates.

        Returns
        -------
        array-like
            The array of outputs.

        """
        from scipy.stats import gaussian_kde

        from pyXMIP.utilities.terminal import Spinner

        mainlog.info("Computing the Gaussian KDE...")
        dphi, dtheta = self.positions

        dphi = np.concatenate(
            [
                dphi + (2 * np.pi),
                dphi - (2 * np.pi),
                dphi + (2 * np.pi),
                dphi - (2 * np.pi),
                dphi,
                dphi,
                dphi,
                dphi + (2 * np.pi),
                dphi - (2 * np.pi),
            ]
        )
        dtheta = np.concatenate(
            [
                dtheta + (np.pi),
                dtheta - (np.pi),
                dtheta,
                dtheta,
                dtheta,
                dtheta - (np.pi),
                dtheta + (np.pi),
                dtheta - (np.pi),
                dtheta + (np.pi),
            ]
        )

        mainlog.debug(
            f"Constructed position buffers of lengths {len(dphi)},{len(dtheta)}."
        )

        with Spinner(text="Computing Gaussian KDE"):
            kernel = gaussian_kde(np.vstack([dphi, dtheta]), weights=self.counts)
            positions = np.vstack([PHI.ravel(), THETA.ravel()])

            Z = np.sum(self.counts) * np.reshape(kernel(positions).T, PHI.shape)

        return Z


class DensityModel:
    r"""
    Abstract class representation of a sky density model :math:`f(\phi,\theta; \mathbf{\Omega})` with a set of model parameters
    specified by :math:`\mathbf{\Omega}`.

    Attributes
    ----------

    class_optimization_methods: dict
        Classifies the various optimization methods available for the specific density model. Each entry in the dictionary
        should be a key of the form ``"<framework>_<method>"`` followed by a subdictionary containing keys ``"level"`` corresponding
        to the priority of the method and ``requires``, corresponding to the necessary additional information the user must specify to take
        advantage of that method. The :py:attr:`~class_optimization_methods` attribute is used to deduce the correct solution method.

        The formatting is as follows:

        .. code-block:: python

            class_optimization_methods = {
                "mle_analytical": {
                    "level": 0 # --> indicates that this method takes top priority.
                    "requires": ["model_function", "analytically_tractable"] #--> for this to be possible, the model must be analytic and it must have a model function.
                }
            }

    coordinate_system: :py:class:`astropy.coordinates.GenericFrame`
        The coordinate system in which this density model is native. Concretely, it is only in this coordinate system that
        :math:`\phi,\theta` correspond to the :math:`\phi,\theta` of the model.

    class_function: callable
        The characteristic density function for the model. This must be a callable instance with signature ``function(phi,theta,arg1,...,argN)``,
        where each ``arg`` is a specific fit parameter of the model.

    class_parameters: dict
        The parameter map for the density model.

        The ``class_parameters`` attribute (if specified) is a dictionary containing the parameters of the model, a description of
        the parameter, the latex representation of the parameter, and units for the parameter. In general, not all of these need
        to be specified and may be filled with ``"None"`` if unused. The dictionary should take the form

        .. code-block:: python

            class_parameters = {
            "parameter_name": ("desc","latex","units")
            }

    """
    coordinate_system = ICRS
    function_string = None
    class_function = None
    class_parameters = None
    class_optimization_approaches = {
        "mle_analytical": {
            "level": 0,
            "requires": ["model_function", "analytically_tractable"],
        },
        "mle_nonlin": {"level": 1, "requires": ["model_function", "grad"]},
    }
    _class_analytically_tractable = False

    def __init__(self, density_function=None, parameters=None, derivative_function=None, second_derivative_function=None, prior_function=None):
        r"""
        Initializes the :py:class:`stats.density.DensityModel` instance.

        Parameters
        ----------
        density_function: callable, optional
            The representative density function. This must be a callable function with signature ``func(phi,theta,arg1,...,argN)``.
            The coordinates ``phi`` and ``theta`` represent spherical coordinates (using the physicist's convention) corresponding to
            model's coordinate system.

            .. warning::

                If the :py:attr:`~density_function` contains unit-ed parameters, the density function **must** be of a form
                such that when the parameter is converted to a float (using units as provided in ``parameters``), the model is
                self-consistent.

        parameters: dict, optional
            A dictionary containing the names of the functions parameters (excluding ``phi`` and ``theta``) with values containing strings
            of a description of the parameter, the latex representation of the parameter, and units for the parameter. In general, not all of these need
            to be specified and may be filled with ``"None"`` if unused. The dictionary should take the form

            .. code-block:: python

                class_parameters = {
                "parameter_name": ("desc","latex","units")
                }

            .. note::

                If ``parameters`` is left unspecified, object inspection is used to determine argument names and
                parameter information is populated with blank descriptions, latex repr., and units.

        derivative_function: callable, optional
            Optional parameter to specify the first order partial derivatives of the density function. The function must have a signature
            ``func(phi,theta,arg1,...,argN)``. It must return an array of shape ``(N,)`` containing the partial derivatives of the density function for
            each of the arguments.

            .. math::

                g(\phi,\theta,a_1,\cdots,a_N)_{i} = \left.\frac{\partial f}{\partial a_i}\right|_{\phi,\theta,a_1,\cdots,a_N}

        second_derivative_function: callable, optional
            Optional parameter to specify the second order partial derivatives of the density function. The function must have a signature
            ``func(phi,theta,arg1,...,argN)``. It must return an array of shape ``(N,N)`` containing the partial derivatives of the density function for
            each of the arguments. Thus,

            .. math::

                g(\phi,\theta,a_1,\cdots,a_N)_{ij} = \left.\frac{\partial^2 f}{\partial a_i \partial a_j}\right|_{\phi,\theta,a_1,\cdots,a_N}
        prior_function: callable, optional
            Optional parameter to specify the Bayesian priors for each of the free parameters in the model. This should be a function
            ``func(arg1,arg2,...,argN)`` which returns a ``float`` indicating the *a priori* likelihood of the parameter set.

            .. warning::

                If you choose to use the ``prior_function`` argument, you **must** remember that the function returns a float, not an array. This means
                that the priors are combined:

                .. math::

                    P(\mathbf{\Omega}|\mathcal{D}) = P(\mathcal{D}|\mathbf{\Omega})\frac{P(\mathbf{\Omega})}{P(\mathcal{D})} \neq P(\mathcal{D}|\mathbf{\Omega})P(\mathcal{D})^{-1} \prod_{i = 1}^{N} P(\mathbf{\Omega}_i)
        """
        self.function = (
            density_function
            if density_function is not None
            else self.__class__.class_function
        )
        self.params = (
            parameters if parameters is not None else self.__class__.class_parameters
        )
        self._fit_info = None

        if parameters is None:
            import inspect

            self.params = {
                k: ("No Description", f"${k}$", "None")
                for k in inspect.getfullargspec(self.function)[0]
                if k not in ["phi", "theta"]
            }

    def __str__(self):
        return f"<{self.__class__.__name__} {self.args}>"

    def __repr__(self):
        return self.__str__()

    def __call__(self, phi, theta):
        return self.function(phi, theta, *self.args)

    def __add__(self, other):
        return self._do_op(self, other, operator.add)

    def __sub__(self, other):
        return self._do_op(self, other, operator.sub)

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        return self.__sub__(other)

    def __mul__(self, other):
        return self._do_op(self, other, operator.mul)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __getitem__(self, item):
        return self.args[item]

    def __setitem__(self, key, value):
        self.args[key] = value

    def __len__(self):
        return len(self.args)

    def __bool__(self):
        return len(self.args) != 0

    @classmethod
    def is_analytic(cls):
        """
        Indicates if the density model has analytical solutions. If ``True``, then there is an analytically tractable solution
        to the parameter determination problem. If ``False``, then numerical methods are required.

        Returns
        -------
        bool

        """
        return cls._class_analytically_tractable

    @classmethod
    def latex(cls):
        """
        Returns the class's native latex representation of the density function.

        .. note::

            Not to be confused with :py:meth:`get_latex`, which will fill parameter values from best fit into the formula.

        Returns
        -------
        str

        """
        return cls.function_string

    def get_latex(self, fill_args=False, formatter=None, **kwargs):
        r"""
        Return the latex representation of the model function with the option of filled arguments.

        Parameters
        ----------
        fill_args: bool, optional
            If ``True``, the formula will display the argument values. These are formatted using the ``formatter`` kwarg.
        formatter: callable, optional
            A function ``f(argument_value,**kwargs)`` which formats each of the argument values into a string to insert into the equation.
        **kwargs
            Additional keyword arguments to pass to the formatter.
        """
        if not formatter:
            formatter = lambda x, **kw: x

        if fill_args:
            output_string = self.function_string % {
                k: f"{v[1]} = {formatter(self.args[k],**kwargs)}"
                for k, v in self.params.items()
            }
        else:
            output_string = self.function_string % {
                k: v[1] for k, v in self.params.items()
            }

        return output_string

    def keys(self):
        """Return the available arguments."""
        return self.args.keys()

    def values(self):
        """Return the available argument values."""
        return self.args.values()

    def items(self):
        """Return iterable of ``(key,value)`` from ``self.args``."""
        return self.args.items()

    @classmethod
    def _do_op(cls, o1, o2, oper):
        if isinstance(o2, Number):
            return cls(
                _do_op(oper, o1.function, lambda x, o=o2: x * 0 + o, len(o1.params)),
                parameters=o1.params,
            )
        elif isinstance(o2, cls):
            return cls(
                _do_op(oper, o1.function, o2.function, len(o1.params)),
                parameters={**o1.params, **o2.params},
            )

    @property
    def args(self):
        """
        The arguments of the model.

        Returns
        -------
        dict
            The arguments of the model. If the model is fit, then ``{param:value}``. Else, ``{param:nan}``.
        """
        if self._fit_info is not None:
            return self._fit_info.args
        else:
            return {k: np.nan for k in self.params}

    @property
    def has_args(self):
        """
        ``True`` if arguments have been fit. ``False`` otherwise.

        Returns
        -------
        bool

        """
        return all([not np.isnan(v) for v in self.args.values()])

    @property
    def fit(self):
        """
        Return the model fit if it has been fit. Otherwise ``None``.

        Returns
        -------
        SimpleNamespace

        """
        if self._fit_info is not None:
            return self._fit_info
        else:
            return None

    @property
    def has_fit(self):
        """
        ``True`` if the model has a fit. ``False`` otherwise.

        Returns
        -------
        bool

        """
        if self._fit_info is not None:
            return True

    @fit.setter
    def fit(self, value):
        self._fit_info = value

    def _get_opt_meths(self):
        """Determines the optimal methods."""
        available_methods = self.class_optimization_approaches

        # -- establish what is provided -- #
        _known_information = []

        if self.function is not None:
            _known_information.append(
                "model_function"
            )  # the function is a known information piece of information.
        if self.__class__._class_analytically_tractable:
            _known_information.append(
                "analytically_tractable"
            )  # The function can be solved analytically.

        # -- search the available methods for a good fit -- #
        available_methods = {
            m: v
            for m, v in available_methods.items()
            if all(r in _known_information for r in v["requires"])
        }
        assert (
            len(available_methods) != 0
        ), f"No available methods were found with known information {_known_information}."

        levels = np.arange(
            np.amax([v["level"] for v in available_methods.values()]) + 1
        )
        methods = []

        for l in levels:
            methods += [m for m, v in available_methods.items() if v["level"] == l]

        return methods

    def print_fit_summary(self, use_rich=xsparams["system"]["add_ons"]["rich"]):
        """
        Print a summary of the model fit.

        .. note::

            For ``rich`` output, you must have the ``rich`` package installed.

        Parameters
        ----------
        use_rich: bool
            If ``True``, will use ``rich`` to print the output.

        Returns
        -------
        None

        """
        if use_rich:
            self._rich_print_fit_summary()
        else:
            self._print_fit_summary()

    def _rich_print_fit_summary(self):
        try:
            from rich.console import Console
            from rich.table import Table
        except ImportError:
            self._print_fit_summary()
            return None

        rich_console = Console()
        rich_console.print(
            "[cyan]|----------------------------------------------------------|",
            justify="center",
        )
        rich_console.print(
            f"[magenta]{self.__class__.__name__}[reset] Fit Summary", justify="center"
        )
        rich_console.print(
            "[cyan]|----------------------------------------------------------|",
            justify="center",
        )

        if not self.has_fit:
            rich_console.print("[bold green]HAS FIT[reset]: [red] FALSE[reset]")
            return None
        else:
            rich_console.print("[bold green]HAS FIT[reset]: [green] TRUE[reset]")

        rich_console.print(
            f"[bold green]Computational Method[reset]: {self.fit.method}"
        )
        rich_console.print(
            f"[bold green]Number of Parameters[reset]: {len(self.params)}"
        )

        # -- creating the table -- #
        param_table = Table(
            title=f"{self.__class__.__name__} -- Fit Parameters",
            highlight=True,
            expand=True,
            style="bold",
            header_style="bold italic red",
        )

        param_table.add_column("Parameter Name", style="bold underline green")
        param_table.add_column("Value", style="")
        param_table.add_column("Description", style="")
        param_table.add_column("Unit", style="")
        param_table.add_column("Variance", style="")

        for param in self.params:
            param_table.add_row(
                param,
                np.format_float_scientific(self.fit.args[param], precision=5),
                self.params[param][0],
                self.params[param][2],
                np.format_float_scientific(self.fit.vargs[param], precision=5),
            )

        rich_console.print(param_table)

    def _print_fit_summary(self):
        print("|----------------------------------------------------------|")
        print(f"\t{self.__class__.__name__} Fit Summary")
        print("|----------------------------------------------------------|")

        if not self.has_fit:
            print("HAS FIT: FALSE")
            return None
        else:
            print("HAS FIT: TRUE")

        print(f"Computational Method: {self.fit.method}")
        print(f"Number of Parameters: {len(self.params)}")

        # -- creating the table -- #
        print(f"{self.__class__.__name__} -- Fit Parameters")

        for param in self.params:
            print(
                f"{param} -- {self.params[param][0]}: {np.format_float_scientific(self.fit.args[param],precision=5)},{np.format_float_scientific(self.fit.vargs[param],precision=5)} "
            )

    def fit_to_population(
        self,
        population,
        *args,
        optimization_method=None,
        inplace=False,
        **kwargs,
    ):
        r"""
        Fit the :py:class:`DensityModel` instance to a provided population.

        Parameters
        ----------
        population: :py:class:`PopulationDensity`
            The population to which this model should be fit.
        *args: optional
            Additional arguments to pass to the optimization method.
        optimization_method: str, optional
            The name of the optimization method to use. If unspecified, it will be deduced.
        inplace: bool, optional
            If ``True``, the fit will be immediately associated with the model instead of simply returned.
        **kwargs: optional
            Additional kwargs to pass to the optimization method.

        Returns
        -------
        SimpleNamespace

        """
        mainlog.info(f"Fitting {self.__class__.__name__} to {len(population)} events.")

        # =========================================================================#
        # DETERMINING THE OPTIMAL FIT METHOD
        # =========================================================================#
        with tqdm_logging_redirect(
            total=100,
            desc="Determining optimal fit method",
            dynamic_ncols=True,
            leave=False,
            bar_format="|{bar}| {percentage:3.0f}% - {desc}",
            loggers=[mainlog],
        ) as pbar:
            # ------------------------------------------------#
            # determining the best fit method
            if optimization_method is not None:
                _opt_meth = [getattr(self, f"_{optimization_method}")]
                mainlog.debug(f"Optimization Method: [FIXED] {optimization_method}.")
            else:
                # determine from the available information what optimization method to use.
                _opt_meths = [
                    getattr(self, f"_{optimization_method}")
                    for optimization_method in self._get_opt_meths()
                ]
                mainlog.debug(
                    "Optimization Methods:\n"
                    + "".join(
                        [
                            f"\t\t {k+1} - {meth.__name__}"
                            for k, meth in enumerate(_opt_meths)
                        ]
                    )
                )

            pbar.update(n=10)
            mainlog.debug("Determining fit method: [COMPLETE]")
            pbar.desc = "Coercing population data..."
            # ------------------------------------------------#
            # Managing the population data.
            phi, theta = (
                population.coordinates.transform_to(self.coordinate_system)
                .frame.spherical.lon.wrap_at(2 * np.pi * units.rad)
                .rad,
                population.coordinates.transform_to(self.coordinate_system)
                .frame.spherical.lat.wrap_at(units.rad * np.pi / 2)
                .rad,
            )
            counts = population.counts
            areas = population.area.to_value(units.sr)
            pbar.update(n=10)
            mainlog.debug("Coercing population data: [COMPLETE]")
            pbar.desc = "Fitting to population..."
            # ------------------------------------------------#
            # Cycling through optimizers
            with tqdm_logging_redirect(
                total=len(_opt_meths),
                leave=False,
                desc="",
                loggers=[mainlog],
            ) as opt_pbar:
                for _opt_meth in _opt_meths:
                    opt_pbar.desc = f"Attempting to fit using {_opt_meth.__name__}..."

                    status, data = _opt_meth(
                        self, phi, theta, counts, areas, self.function, *args, **kwargs
                    )

                    if status:
                        data["method"] = _opt_meth.__name__
                        continue

                    opt_pbar.update()
                    mainlog.warning(
                        f"Failed to fit using {_opt_meth.__name__}. (REASON: {data})"
                    )
            assert (
                status
            ), "The optimization failed after exhausting available optimization schemes."

            pbar.desc = f"[COMPLETE] (Method = {data['method']})"
            pbar.update(n=80)

        # --------------------------------------------------#
        # Setting optimization parameters.
        if inplace:
            self.fit = SimpleNamespace(**data)

        mainlog.info(f"Completed fitting {self.__class__.__name__} to {population}.")
        return self.fit

    def _mle_nonlinear(
        self, phi, theta, counts, areas, function, guess_functions=None, ls_kwargs=None
    ):
        r"""
        Determines best fit parameters of the model for the provided details from a maximum likelihood estimation predicated
        on a poisson distribution.

        Parameters
        ----------
        phi: array-like
            Array of length ``N`` containing the latitude points of the data in the relevant coordinate system.
        theta: array-like
            Array of length ``N`` containing the longitudinal points of the data in the relevant coordinate system.
        counts: array-like
            Array of length ``N`` containing the counts for each position.
        areas: array-like
            Array of length ``N`` of solid angles from which each count sample was pulled. If not a unit carrying array,
            the units are assumed to be steradians.
        function: callable
            The density function of the Poisson distribution:

            .. math::

                \lambda_i = A_if(\phi_i,\theta_i,\textbf{P})

            The function must have signature ``func(phi,theta,param1,param2,...,paramN)``.
        guess_functions: Collection of callable
            The guess functions (``f_i(phi,theta)``) for each of the estimators. These should each return a float, which is
            used as the initial guess for the MLE minimization.
        ls_kwargs: dict
            Additional kwargs to pass to ``scipy.optimize.least_squares``.
        """
        from scipy.optimize import least_squares

        if not ls_kwargs:
            ls_kwargs = {}
        # -- building initial guess -- #
        assert guess_functions is None or (
            len(guess_functions) == len(self.params)
        ), f"There are too many guess functions {len(guess_functions)} for the model (N parameters = {len(self.params)})."

        try:
            assert guess_functions is not None
            guess = [f(phi, theta) for f in guess_functions]
        except AssertionError:
            mainlog.warning(
                "Guess functions were not supplied for the NL-LS algorithm. Assuming best guess is 0."
            )
            guess = [0 for _ in self.params]

        # -- setting up the minimizations -- #
        try:
            # defining the cost function.
            cost_function = lambda x, f=function: f(phi, theta, *x) - (counts / areas)

            # running the minimization.
            fit_info = least_squares(cost_function, guess, **ls_kwargs)

            mainlog.info(
                f"Fit {self.__class__.__name__} to {len(phi)} points via Poisson MLE least squares [Non-Linear]."
                + self._parameter_summary_table(fit_info.x)
            )

        except Exception as exception:
            mainlog.info(f"P-MLE NLLS Failed. [{exception.__repr__()}]")

            return False, None

        return True, {
            "args": {k: v for k, v in zip(self.params, fit_info.x)},
            "full": fit_info,
            "vargs": {k: np.nan for k in self.params},
            "cov": None,
        }

    @_enforce_style
    def plot3d(
        self,
        arguments=None,
        ax=None,
        fig=None,
        gridsize=100,
        cmap=None,
        vmin=None,
        vmax=None,
        norm=None,
        colorbar=True,
        coordinate_system=None,
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
        coordinate_system: :py:class:`astropy.coordinates.GenericFrame`, optional
            The coordinate system to plot in.
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
        # Argument Management

        if self.has_args:
            arguments = list(self.args.values())
        else:
            arguments = arguments

        if arguments is None:
            arguments = []

        # ----------------------------------------------#
        # Computing the grid
        phi, theta = np.mgrid[
            0 : 2 * np.pi : gridsize * 1j, -np.pi / 2 : np.pi / 2 : gridsize * 1j
        ]
        coordinates = SkyCoord(phi, theta, unit="rad", frame=self.coordinate_system)

        # converting to new coordinates.
        if coordinate_system is not None:
            phi = (
                getattr(coordinates, coordinate_system.name)
                .frame.spherical.lon.wrap_at(2 * np.pi * units.rad)
                .rad
            )
            theta = (
                getattr(coordinates, coordinate_system.name)
                .frame.spherical.lat.wrap_at(units.rad * np.pi / 2)
                .rad
            )

        X, Y, Z = (
            np.cos(theta) * np.cos(phi),
            np.cos(theta) * np.sin(phi),
            np.sin(theta),
        )

        # Image generation
        image_array = self.function(phi, theta, *arguments)
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
            cbar = plt.colorbar(cbar, ax=ax, label=self.get_latex())

        return fig, ax, cbar

    @_enforce_style
    def plot2d(
        self,
        arguments=None,
        fig=None,
        coordinate_system=None,
        gridsize=64,
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
        arguments: array-like, optional
            Arguments to pass to the density model. If :py:attr:`DensityModel.args` exists, this option is
            overridden by the set arguments.
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
        coordinate_system: :py:class:`astropy.coordinates.GenericFrame`, optional
            The coordinate system to plot in.

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

        if coordinate_system is None:
            coordinate_system = self.coordinate_system

        # -- managing arguments -- #
        if self.has_args:
            arguments = list(self.args.values())
        else:
            arguments = arguments

        if arguments is None:
            arguments = []

        # -- constructing grid and plotting -- #
        _number_of_pix = hp.nside2npix(gridsize)
        _t, _p = hp.pix2ang(gridsize, np.arange(_number_of_pix))

        # transform _t and _p to mathematical coordinates.
        _t -= np.pi / 2

        # construct the theta and phi
        coordinates = SkyCoord(_p, _t, unit="rad", frame=self.coordinate_system)

        # convert
        _p = (
            getattr(coordinates, coordinate_system.name)
            .frame.spherical.lon.wrap_at(2 * np.pi * units.rad)
            .rad
        )
        _t = (
            getattr(coordinates, coordinate_system.name)
            .frame.spherical.lat.wrap_at(units.rad * np.pi / 2)
            .rad
        )
        image_values = self.function(_p, _t, *arguments)
        # -- Building the figure -- #
        if not fig:
            fig = plt.figure(figsize=(5, 4))
            fig.subplots_adjust(left=0, right=0.95, top=0.92, bottom=0.03)

        hp.mollview(image_values, fig=fig, cbar=False, notext=True, title=None)
        hp.graticule()

        ax = fig.gca()

        if cmap is None:
            cmap = plt.colormaps[plt.rcParams["image.cmap"]]

        if vmin is None:
            vmin = np.amin(image_values)

        if vmax is None:
            vmax = np.amax(image_values)

        if norm is None:
            from matplotlib.colors import Normalize

            norm = Normalize(vmin=vmin, vmax=vmax)

        cbar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        if colorbar:
            cbar = plt.colorbar(cbar, ax=ax, label=self.get_latex())

        return fig, ax, cbar


class UniformDensityModel(DensityModel):
    r"""
    Density model with number density following

    .. math::

        \hat{n}(\phi,\theta) = \alpha

    """

    coordinate_system = ICRS
    function_string = (
        r"$\hat{N}(\phi,\theta;\alpha) = \alpha\; / \;\left[\mathrm{Srd}^{-1}\right]$"
    )
    class_function = lambda phi, theta, alpha: alpha * np.ones(np.array(phi).shape)
    class_parameters = {
        "alpha": ("Density parameter", r"$\alpha$", str(units.rad**-2))
    }
    class_optimization_approaches = {
        "poisson_mle_analytical": {
            "level": 0,
            "requires": ["model_function", "analytically_tractable"],
        },
        "poisson_mle_nonlin": {"level": 1, "requires": ["model_function", "grad"]},
    }
    _class_analytically_tractable = True

    def __init__(self):
        super().__init__(
            self.__class__.class_function, parameters=self.__class__.class_parameters
        )

    def _mle_analytical(
        self, phi, theta, counts, areas, function, *args, **kwargs
    ):
        _, _, _, _, _ = (
            phi,
            theta,
            function,
            args,
            kwargs,
        )  # --> prevent issues with flake8
        return True, {
            "args": {"alpha": np.sum(counts) / np.sum(areas)},
            "full": None,
            "vargs": {"alpha": np.sum(counts) / (np.sum(areas) ** 2)},
            "cov": np.array([[np.sum(counts) / (np.sum(areas) ** 2)]]),
        }


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    ra, dec = np.random.uniform(-180, 180, 1000), np.random.uniform(-90, 90, 1000)
    coords = SkyCoord(ra=ra, dec=dec, unit="deg")
    counts = np.random.randint(1, 100, 1000)
    h = PopulationDensity(coords, counts)
    h.fit_model(inplace=True)
    h.model.print_fit_summary(use_rich=True)
    func = lambda phi, theta, a, b, c: a * np.exp(-b * (np.abs(theta) - c))
    phi, theta = np.random.uniform(0, 2 * np.pi, 100), np.random.uniform(
        -np.pi / 2, np.pi / 2, 100
    )
    counts = np.random.randint(0, 10, 100)
