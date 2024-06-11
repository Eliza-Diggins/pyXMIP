"""
Regression tools for pyXMIP
"""
from types import SimpleNamespace

import numpy as np
from astropy.coordinates import SkyCoord
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsRegressor, RadiusNeighborsRegressor
from tqdm.contrib.concurrent import process_map

from pyXMIP.utilities.geo import convert_skycoord
from pyXMIP.utilities.logging import mainlog


class PoissonMapRegressor:
    """
    Generic wrapper class for :py:mod:`sklearn` regression classes for generating :py:class:`pyXMIP.structures.map.Map` instances.

    Attributes
    ----------

    regressor_parameters: dict
        The complete set of parameters for this regressor class. This includes hyperparameters and also
        parameters which are not permitted to vary.
    default_training_fraction: float
        The standard training fraction for the regressor.
    hyperparameters: dict
        Subset of ``regressor_parameters`` containing only the (variation permissible) hyperparameters of the regressor.
    """

    regressor_parameters = {}
    default_training_fraction = 0.5
    hyperparameters = {}

    def __init__(self, regressor_class, **kwargs):
        """
        Initialize the :py:class:`PoissonMapRegressor` instance.

        Parameters
        ----------
        regressor_class:
            The :py:mod:`sklearn` specific regressor type to build this class around.
        **kwargs
            regression parameters.
        """
        self.regressor = regressor_class(
            **{
                k: (v if k not in kwargs else kwargs[k])
                for k, v in self.__class__.regressor_parameters.items()
            }
        )

    def train_model(self, count_table, object_type, training_fraction=None, **kwargs):
        """
        Train the regression model on a ``COUNTS`` table.

        Parameters
        ----------
        count_table: :py:class:`astropy.table.Table`
            The table of counts to use for the training.
        object_type: str
            The specific object type to model.
        training_fraction: float
            The training fraction.
        **kwargs
            Additional kwargs.
        """
        # -- Setup and sanity check -- #
        training_fraction = (
            training_fraction
            if training_fraction is not None
            else self.__class__.default_training_fraction
        )
        assert (
            object_type in count_table.columns
        ), f"Failed to find {object_type} in provided table."

        # -- Fetching data -- #
        X, Y = self._sanitize_points(count_table, object_type)
        mainlog.debug(
            f"Training {self.regressor} on {int(training_fraction * Y.size)} data-points"
        )

        positions_train, positions_test, values_train, values_test = train_test_split(
            X,
            Y,
            train_size=training_fraction,
            random_state=kwargs.pop("random_state", None),
        )
        self.regressor.fit(positions_train, values_train)

        if not kwargs.pop("return_all", False):
            return SimpleNamespace(
                score=self.regressor.score(positions_test, values_test),
                x_train=None,
                x_test=None,
                y_train=None,
                y_test=None,
            )
        else:
            return SimpleNamespace(
                score=self.regressor.score(positions_test, values_test),
                x_train=positions_train,
                x_test=positions_test,
                y_train=values_train,
                y_test=values_test,
            )

    def cross_validate(
        self, count_table, object_type, training_kw=None, param_kw=None, **kwargs
    ):
        """
        Cross-validate the regressor against the available count data.

        Parameters
        ----------
        count_table: :py:class:`astropy.table.Table`
            The table of data from ``COUNTS``.
        object_type: str
            The object type to pull from the table.
        training_kw: dict
            Training kwargs.
        param_kw: dict
            Parameter kwargs
        kwargs:
            Additional kwargs to pass.
        """
        from sklearn.model_selection import GridSearchCV

        mainlog.debug(f"Cross-validating {self}.")
        # -- SETUP -- #
        # grab data
        X, Y = self._sanitize_points(count_table, object_type)

        # fix parameters
        if training_kw is None:
            training_kw = {
                "training_fraction": self.__class__.default_training_fraction,
                "cv": 5,
            }

        if param_kw is None:
            param_kw = self.__class__.hyperparameters

        # split data
        positions_train, positions_test, values_train, values_test = train_test_split(
            X, Y, train_size=training_kw.get("training_fraction")
        )

        # -- RUN GS -- #
        cv_estimator = GridSearchCV(self.regressor, param_kw, **kwargs)
        cv_estimator.fit(positions_train, values_train)

        return cv_estimator.cv_results_

    def _sanitize_points(self, count_table, object_type):
        """Pull count data and fix issues."""
        RA, DEC, N, R = (
            np.array(count_table["RA"]),
            np.array(count_table["DEC"]),
            np.array(count_table[object_type]),
            np.array(count_table["RAD"]),
        )

        # RA and DEC need to be converted to standardized lat-lon. USE MATH CONVENTION
        positions = SkyCoord(ra=RA, dec=DEC, unit="deg")
        phi, theta = convert_skycoord(positions, "math_convention")
        # phi,theta = np.deg2rad(RA),np.deg2rad(DEC)
        # Fixing R and getting areas.
        # NOTE: we always use N/arcmin^2. The RAD is always already in ARCMIN.
        A = np.pi * (R**2)
        density = N / A

        X = np.vstack([theta, phi]).T  # --> Enforced haversine convention.
        Y = density

        return X, Y


class BayesianPoissonMapRegressor:
    r"""
    The :py:class:`BaysianPoissonMapRegressor` class provides a variety of posterior estimation methods for estimating
    sky-map densities.

    Notes
    -----

    The :py:class:`BaysianPoissonMapRegressor` takes 2 kwargs: ``prior`` and ``model``. Generically, the posterior distribution is

    .. math::

        P(\beta|f,\textbf{r},\textbf{N},\textbf{A}) = \frac{p(\beta)}{p(\textbf{N}|\textbf{A},\textbf{r},\beta)} \prod_{i=1}^{i = |N|} \frac{(A_if(\textbf{r_i},\beta))^{N_i} \exp\left(-A_if(\textbf{r_i},\beta)\right)}{N_i!}

    Thus, the log posterior is

    .. math::

        \log P(\beta|f,\textbf{r},\textbf{N},\textbf{A}) = p(\beta) + \sum_{i=1}^{i=|N|} N_i \log(A_if(\textbf{r_i},\beta)) + A_if(\textbf{r_i},\beta) - N_i! + C.

    If no model is specified (``model=None``), then the posterior instead takes the form

    .. math::

        \log P(\lambda_k|f,\textbf{r}_k,\textbf{N}_k,\textbf{A}_k) = p(\lambda_k) + \sum_{i=1}^{i=|N_k|} N_i \log(A_if(\textbf{r_i},\beta)) + A_if(\textbf{r_i},\beta) - N_i! + C.

    Where :math:`\lambda_k` is the rate in each of the HEALPix cells. Thus, when no model is specified, the problem is discretized between cells.
    """

    def __init__(self, prior=None, model=None):
        self.model = model
        self.prior = prior

    def build_map_MAP(self, count_table, object_type, **kwargs):
        # ======================================= #
        # Setup
        # ======================================= #
        mainlog.debug(f"Building MAP estimate for {self}.")

        # ======================================= #
        # organize the pipeline, select method
        # ======================================= #
        if (self.model is None) and (self.prior is None):
            # We can get away with an analytical discretization
            mainlog.debug("MAP estimator -- ANALYTIC [no p, no m]")
            return self._build_map_MAP_analytic(count_table, object_type, **kwargs)
        else:
            raise NotImplementedError(
                "Cannot currently generate MAP estimates when priors / model are used."
            )

    def _build_map_MAP_analytic(self, count_table, object_type, **kwargs):
        _ = kwargs
        count_table = count_table
        # ==================================== #
        # Discretizing data observations
        # ==================================== #
        PIX_ID, N, RAD = (
            count_table["PIX_ID"],
            count_table[object_type],
            count_table["RAD"],
        )
        A = np.array(np.pi * (RAD**2), dtype="float64")
        PIX_ID = np.array(PIX_ID, dtype="uint32")
        N = np.array(N, dtype="uint32")

        # ==================================== #
        # Getting Estimates
        # ==================================== #
        from itertools import repeat

        from pyXMIP.utilities.optimize import created_shared_memory_equivalent, split

        estimates = np.zeros(
            np.amax(PIX_ID) + 1, dtype="float64"
        )  # generate PIX_ID count. ! Not necessarily all HP, might cut off end.
        mp_kwargs = kwargs.pop(
            "mp_kwargs", {"multiprocess": True, "nproc": 8, "chunksize": 100}
        )
        if not mp_kwargs.get("multiprocess", True):
            mp_kwargs = {"multiprocess": False, "nproc": 1, "chunksize": len(estimates)}
        # Produce the shared memory regions #
        SMI, SMN, SMA, SME = (
            created_shared_memory_equivalent(PIX_ID),
            created_shared_memory_equivalent(N),
            created_shared_memory_equivalent(A),
            created_shared_memory_equivalent(estimates),
        )
        # Chunk generation #
        PIX_GROUPS = split(PIX_ID, len(PIX_ID) // mp_kwargs.get("chunksize", 100))
        # run #
        process_map(
            self._get_map_analytic_estimates_mp,
            PIX_GROUPS,
            repeat(SMI),
            repeat(SMN),
            repeat(SMA),
            repeat(SME),
            max_workers=mp_kwargs.get("nproc", 8),
            desc="Computing MAP estimate (ANYTC)",
        )

        # recover shared memory regions #
        ee = np.frombuffer(SME.buf)[:]
        estimates = ee.copy()
        del ee
        for u in [SMI, SMN, SMA, SME]:
            u.close()
            u.unlink()
        return estimates

    @staticmethod
    def _get_map_analytic_estimates_mp(pxg, PIX_ID, N, A, E):
        n, a, ee, i = (
            np.frombuffer(N.buf, dtype="uint32"),
            np.frombuffer(A.buf),
            np.frombuffer(E.buf),
            np.frombuffer(PIX_ID.buf, dtype="uint32"),
        )

        for pi in pxg:
            if pi in i:
                ee[pi] = np.sum(n[np.where(i == pi)]) / np.sum(a[np.where(i == pi)])


class KNNeighborMapRegressor(PoissonMapRegressor):
    regression_class = KNeighborsRegressor
    regressor_parameters = {
        "n_neighbors": 10,
        "weights": "distance",
        "algorithm": "auto",
        "metric": "haversine",
    }
    hyperparameters = {"n_neighbors": np.array([5, 10, 30, 60, 120, 300])}

    def __init__(self, **kwargs):
        super().__init__(self.__class__.regression_class, **kwargs)


class RNNeighborMapRegressor(PoissonMapRegressor):
    regression_class = RadiusNeighborsRegressor
    regressor_parameters = {
        "radius": 0.01,
        "weights": "distance",
        "algorithm": "auto",
        "metric": "haversine",
    }
    hyperparameters = {"radius": np.array([0.001, 0.005, 0.01, 0.05, 0.1])}

    def __init__(self, **kwargs):
        super().__init__(self.__class__.regression_class, **kwargs)


def train_map_regressor(
    regressor, map_phi, map_theta, map_target_values, training_fraction=0.5, **kwargs
):
    """
    Train the map regressor ``regressor`` on the provided positional data.

    Parameters
    ----------
    regressor
    map_phi
    map_theta
    map_target_values
    kwargs

    Returns
    -------

    """
    mainlog.debug(
        f"Training {regressor} on {int(training_fraction*len(map_target_values))} data-points"
    )
    positions = np.vstack(map_phi, map_theta).T  # ready for sklearn now.
    positions_train, positions_test, values_train, values_test = train_test_split(
        positions,
        map_target_values,
        train_size=training_fraction,
        random_state=kwargs.pop("random_state", None),
    )

    regressor.fit(positions_train, values_train)
    return regressor.score(positions_test, values_test)
