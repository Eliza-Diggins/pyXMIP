"""
Reduction processes module for ``pyXMIP``.
"""
from abc import ABC


class ReductionProcessParameter:
    """
    Descriptor representation of a standard reduction-process setting.
    """

    def __init__(
        self,
        default=None,
        required=False,
        enabled=True,
        expected_types=None,
        allowed_values=None,
    ):
        self.default = default
        self._is_required = required
        self._is_enabled = enabled
        self._expected_types = expected_types
        self._allowed_values = allowed_values

    def __set_name__(self, owner, name):
        self._name = name

    def __get__(self, instance, owner):
        if self._name in instance._parameter_dictionary:
            return instance._parameter_dictionary[self._name]
        else:
            return self.default

    def __set__(self, instance, value):
        instance._parameter_dictionary[self._name] = value

        self.validate()

    def validate(self):
        pass


class ReductionProcess(ABC):
    """
    Abstract class implementation of a generic reduction process.
    """

    def __init__(self, process_name, **kwargs):
        self.name = process_name

        self._parameter_dictionary = kwargs

    def __call__(self, cross_match_database):
        """
        Run the reduction process contained by this class instance on a given cross-matching database.

        Parameters
        ----------
        cross_match_database: str or :py:class:`cross_reference.CrossMatchDatabase`
            The cross-matching database to run the reduction on. If `str` is possed, then it is assumed to be the corresponding path
            to the SQL database containing the cross-matching database. Otherwise, it is assumed that this is a :py:class:`cross_reference.CrossMatchDatabase` instance
            in its own right.

        Returns
        -------
        None
        """

    @classmethod
    def _get_parameter_descriptor_instances(cls):
        """Searches and returns all descriptors attached to this class or inherited."""
        _descriptor_instances = {
            k: v
            for k, v in cls.__dict__.items()
            if isinstance(v, ReductionProcessParameter)
        }

        # iterate through all parent classes.
        _bases = cls.__bases__  # list of all parents.
        while object not in _bases:
            _nbases = []
            for _base in _bases:
                _descriptor_instances = {
                    **_descriptor_instances,
                    **{
                        k: v
                        for k, v in _base.__dict__.items()
                        if isinstance(v, ReductionProcessParameter)
                    },
                }
                _nbases += _base.__bases__

            _bases = _nbases[:]

        return _descriptor_instances


# class _REDUCTION_DESCRIPTOR:
#    """
#    Descriptor class for reduction process parameters.
#    """
#
#    def __init__(
#        self, dict_location, default=None, required=False, dtype=None, allow_set=True
#    ):
#        self.default = default
#        self.dict_location = dict_location
#        self.required = required
#        self.dtype = dtype
#        self.allow_set = allow_set
#
#    def __set_name__(self, owner, name):
#        self._name = name
#        self.dict_location += [self._name]
#
#    def __get__(self, instance, owner):
#        try:
#            return getFromDict(instance._kwargs, self.dict_location)
#        except KeyError:
#            return self.default
#
#    def __set__(self, instance, value):
#        setInDict(instance._kwargs, self.dict_location, value)
#
#
# class ReductionProcess:
#    """
#    Base class for all ``pyXMIP`` :py:class:`ReductionProcess` types.
#
#    .. warning::
#
#        This class is abstract and should not be used except for development purposes.
#
#    Attributes
#    ----------
#    DB_REGISTRY: :py:class:`structures.databases.DBRegistry`
#        The database registry to use for looking up database classes by name. By default, ``DB_REGISTRY`` is just
#        ``DEFAULT_DATABASE_REGISTRY``.
#
#        .. hint::
#
#            If you're trying to run reductions on custom databases which are not part of the main codebase,
#            this registry should include your custom database.
#    CHUNKSIZE: int
#        The chunksize of SQL operations. Lower chunksizes increase memory efficiency while decreasing speed; higher
#        chunksize decreases memory efficiency but increases speed. By default, ``CHUNKSIZE=10000``.
#
#
#    Notes
#    -----
#
#    """
#
#    _process_meta_tag = "ABSTRACT"
#    _output_column_name = "ABSTRACT_RES"
#
#    DB_REGISTRY = DEFAULT_DATABASE_REGISTRY
#    CHUNKSIZE = _REDUCTION_DESCRIPTOR([], default=10000, dtype=int, allow_set=True)
#
#    def __init__(self, **kwargs):
#        """
#        Initialize a reduction process.
#
#        Parameters
#        ----------
#        **kwargs:
#            keyword arguments to specify for the reduction process.
#        """
#        mainlog.debug(f"[{self.__class__.__name__}] Initializing instance.")
#
#        # -- Manage kwarg loading -- #
#        self._kwargs = {}  # create empty kwarg dict to fill.
#        self._manage_kwarg_loading(kwargs)  # fill karg dict.
#
#    def __call__(self, cross_match_database, table, overwrite=False):
#        """
#        Run the reduction process contained in this class on a given cross matching table.
#
#        Parameters
#        ----------
#        cross_match_database: CrossMatchDatabase
#            The :py:class:`CrossMatchDatabase` object to reduce.
#        table: str
#            The table to run the reduction on. Should be ``<Database_Name>_MATCH``.
#        overwrite: bool
#            If ``True``, then the process is run and will overwrite the existing data. Otherwise,
#            it will fail.
#
#        Returns
#        -------
#        None
#        """
#        # =================================== #
#        # Setup                               #
#        # =================================== #
#        mainlog.info(
#            f"[{self.__class__.__name__}] Running on {cross_match_database} table {table}."
#        )
#
#        # -- Manage the meta status -- #
#        meta_status = cross_match_database.check_meta(self._process_meta_tag, table)
#        if meta_status and overwrite:
#            mainlog.debug(
#                f"[{self.__class__.__name__}] overwriting existing data in {table}."
#            )
#
#            with cross_match_database.connect() as conn:
#                conn.execute(
#                    sql.text(f"ALTER TABLE {table} DROP {self._output_column_name}")
#                )
#            mainlog.info(
#                f"[{self.__class__.__name__}] DROPPED {self._output_column_name} FROM {table}."
#            )
#        elif meta_status and not overwrite:
#            raise ValueError(
#                f"[{self.__class__.__name__}] {cross_match_database} META indicates this reduction already run, and overwrite=False."
#            )
#        else:
#            pass
#
#        # =================================== #
#        # Run reduction process               #
#        # =================================== #
#        # Setup
#        vals = self.setup_run_reduction(cross_match_database, table, **self._kwargs)
#
#        # run
#        self.run_reduction(cross_match_database, table, *vals)
#        # =================================== #
#        # Manage META                         #
#        # =================================== #
#        cross_match_database.meta_add(self._process_meta_tag, table)
#
#    def __str__(self):
#        return f"<{self.__class__.__name__} Instance>"
#
#    def __repr__(self):
#        return f"<{self.__class__.__name__} Instance, CHUNKSIZE = {self.CHUNKSIZE}>"
#
#    def _manage_kwarg_loading(self, kwargs):
#        for p in self.get_parameter_names():
#            if p in self.__class__.__dict__:
#                p_cls = self.__class__.__dict__[p]
#            else:
#                # necessary for loading from superclass. Leads to inheritance depth limit.
#                p_cls = self.__class__.__base__.__dict__[p]
#
#            try:
#                setattr(self, p, kwargs[p])
#            except (KeyError, AttributeError):
#                # that attribute doesn't have a value.
#                if p_cls.required:
#                    raise ValueError(f"The kwarg {p} is required.")
#
#        for k, v in kwargs.items():
#            if k not in self.get_parameter_names() and v is not None:
#                mainlog.warning(
#                    f"Detected kwarg {k}, which is not recognized. Check that it is spelled correctly."
#                )
#
#    @property
#    def parameters(self):
#        """Dictionary of the parameters and their values for this object."""
#        return {m: getattr(self, m) for m in self._get_descriptors()}
#
#    @classmethod
#    def _get_descriptors(cls):
#        # fetching descriptors.
#        self_instances = [
#            m for m, v in cls.__dict__.items() if isinstance(v, _REDUCTION_DESCRIPTOR)
#        ]
#        other_instances = [
#            m
#            for m, v in cls.__base__.__dict__.items()
#            if isinstance(v, _REDUCTION_DESCRIPTOR)
#        ]
#
#        return self_instances + other_instances
#
#    @classmethod
#    def get_parameter_names(cls):
#        """Get the names of the available parameters."""
#        return cls._get_descriptors()
#
#    def _enforce_parameters(self):
#        for k, v in self.parameters.items():
#            if v is None and self.__class__.__dict__[k].required:
#                raise ValueError(f"Parameter {k} is not set.")
#
#    @classmethod
#    def setup_run_reduction(cls, cross_match_database, table, **kwargs):
#        """
#        Setup a reduction run. Provides a ``tuple`` of the arguments needed to run the reduction.
#
#        Parameters
#        ----------
#        cross_match_database: :py:class:`cross_reference.CrossMatchDatabase`
#            The cross matching database to run on.
#        table: str
#            The table to run on.
#        **kwargs:
#            Keyword arguments for the necessary parameters in the reduction.
#
#        Returns
#        -------
#        tuple:
#            Ordered tuple of results to pass into :py:meth:`ReductionProcess.run_reduction`.
#
#        """
#        _, _ = cross_match_database, table
#        return tuple(kwargs.values())
#
#    @classmethod
#    def run_reduction(cls, cross_match_database, table, *args):
#        """
#        Run a reduction process.
#
#        Parameters
#        ----------
#        cross_match_database: :py:class:`cross_reference.CrossMatchDatabase`
#            The cross matching database to run on.
#        table: str
#            The table to run on.
#        *args:
#            The arguments for the reduction run. These can be obtained from :py:meth:`ReductionProcess.setup_run_reduction`.
#
#        Returns
#        -------
#        None
#
#        """
#        with logging_redirect_tqdm(loggers=[mainlog]):
#            # -- Setting up the temporary database names -- #
#            _temp_database_table_name = (
#                f"{table}_TMP"  # Where we are putting the full output.
#            )
#
#            # ================================================== #
#            # Setup chunking
#            # ================================================== #
#            _chunksize, args = args[0], args[1:]
#
#            with cross_match_database.connect() as conn:
#                _table_size = conn.execute(
#                    sql.text(f"SELECT COUNT(*) FROM {table}")
#                ).scalar()
#
#            offsets = np.arange(0, _table_size, _chunksize, dtype="uint32")
#            _chunk_pbar = tqdm(
#                total=_table_size, desc="Chunking Execution", leave=False
#            )
#            # ================================================== #
#            # Iterative run.
#            # ================================================== #
#            for offset in offsets:
#                cls._run_reduction(
#                    cross_match_database,
#                    table,
#                    _temp_database_table_name,
#                    offset,
#                    _chunksize,
#                    _chunk_pbar,
#                    *args,
#                )
#
#            # ================================================== #
#            # Drop old tables.
#            # ================================================== #
#            with cross_match_database.connect() as conn:
#                conn.execute(sql.text(f"DROP TABLE {table}"))
#                mainlog.debug(f"DROPPED TABLE {table}.")
#                conn.execute(
#                    sql.text(
#                        f"ALTER TABLE {_temp_database_table_name} RENAME TO {table}"
#                    )
#                )
#                mainlog.debug(
#                    f"ALTERED TABLE {_temp_database_table_name}, RENAMED AS {table}."
#                )
#
#    @classmethod
#    def _run_reduction(
#        cls, cross_match_database, table, temp_table, offset, chunksize, pbar, *args
#    ):
#        pass
#
#
# class PSFReductionProcess(ReductionProcess):
#    """
#    Reduction process to compute the probability of a match at a distance :math:`d` given the PSF of the instrument.
#
#    Parameters
#    ----------
#    PSF_FUNCTION: callable, optional
#        The functional representation of the PSF in the form ``f(x,scale)``.
#        :math:`x` is the distance from the source, and ``scale`` is the characteristic length scale. By default, this is a Gaussian.
#
#        .. warning::
#
#            At this point in development, the PSF cannot take any additional parameters aside from :math:`x`.
#
#    PSF_SCALE: Quantity, optional
#        The scale length for the PSF. Only applies if ``PSF_FUNCTION`` is not specified, in which case, the default Gaussian PDF is used and
#        the ``PSF_SCALE`` parameter is the standard deviation of the distribution.
#
#    EXCLUSION_CRITERIA: tuple, optional
#        An optional parameter to exclude some objects from this reduction (giving them ``score=0``). The tuple should be a valid
#        column in ``CATALOG`` and a callable returning a boolean value for exclusion.
#
#        .. hint::
#
#            Typically, this is most useful for excluding extended sources. In that case, you might have a column with a likelihood of
#            being extended (e.g. ``EXT_LIKE``). You could then set a threshold beyond which a source is deemed extended and automatically
#            given a perfect score for this reduction by setting ``EXCLUSION_CRITERIA = ("EXT_LIKE",lambda x: x>3)``.
#
#    Notes
#    -----
#
#    The basic principle behind the PSF reduction is that for a point source at true sky location :math:`\mathbf{r}`, the instrumental
#    PSF will spread the source out such that the source photons appear to be distributed according to
#
#    .. math::
#
#    P(|\textbf{r}'-\textbf{r}|) \sim f(|\textbf{r}'-\textbf{r}|)
#
#    This :math:`f(x)` is effectively the PSF. As such, we might observe a "detection" anywhere within the PSF with a given
#    probability. That is the functionality of this class.
#
#    """
#
#    _process_meta_tag = "PSF_REDUCTION"
#    _output_column_name = "PSF_REDUCTION_RES"
#
#    _default_psf = lambda x, scale: np.exp(-0.5 * (x / scale) ** 2)
#
#    # -- Parameters -- #
#    PSF_FUNCTION = _REDUCTION_DESCRIPTOR(
#        [], required=False, allow_set=True, default=_default_psf
#    )
#    PSF_SCALE = _REDUCTION_DESCRIPTOR([], required=True, allow_set=True)
#    EXCLUSION_CRITERIA = _REDUCTION_DESCRIPTOR(
#        [], required=False, allow_set=True, default=None
#    )
#
#    def __init__(self, **kwargs):
#        super().__init__(**kwargs)
#
#    @classmethod
#    def setup_run_reduction(cls, cross_match_database, table, **kwargs):
#        # ================================================== #
#        # Setup the run settings                             #
#        # ================================================== #
#        _, _ = cross_match_database, table  # --> Trick IDE.
#
#        EXCLUSION_CRITERIA = kwargs.pop("EXCLUSION_CRITERIA", None)
#        PSF_FUNCTION = kwargs.pop("PSF_FUNCTION", None)
#        PSF_SCALE = enforce_units(kwargs["PSF_SCALE"], "arcmin").to_value("arcmin")
#
#        if PSF_FUNCTION is None:
#            PSF_FUNCTION = cls._default_psf
#
#        # Make sure that the psf function is a single variable function for ease.
#        _psf_function = lambda x, q=PSF_SCALE: PSF_FUNCTION(x, q)
#
#        if EXCLUSION_CRITERIA is not None:
#            _exclusion_column, _exclusion_callable = EXCLUSION_CRITERIA
#        else:
#            # set a trivial exclusion system.
#            _exclusion_column, _exclusion_callable = "CATALOG_OBJECT", lambda x: 0
#
#        return (
#            kwargs.get("CHUNKSIZE", 1000),
#            _psf_function,
#            _exclusion_column,
#            _exclusion_callable,
#        )
#
#    @classmethod
#    def _run_reduction(
#        cls, cross_match_database, table, temp_table, offset, chunksize, pbar, *args
#    ):
#        # =============================================== #
#        # Setup                                           #
#        # =============================================== #
#        psf, exc_c, exc_call = args
#        sql_call = "SELECT %(table)s.*,CATALOG.%(ext_col)s FROM %(table)s INNER JOIN CATALOG USING(CATALOG_OBJECT) LIMIT %(chunksize)s OFFSET %(offset)s"
#
#        # =============================================== #
#        # Runtime                                         #
#        # =============================================== #
#        with cross_match_database.connect() as conn:
#            _memory_reference_table = pd.read_sql(
#                sql.text(
#                    sql_call
#                    % dict(
#                        table=table, chunksize=chunksize, offset=offset, ext_col=exc_c
#                    )
#                ),
#                conn,
#            )
#
#            # -- calculate separations -- #
#            positions_reference = SkyCoord(
#                ra=_memory_reference_table["RA"],
#                dec=_memory_reference_table["DEC"],
#                unit="deg",
#            )
#            positions_catalog = SkyCoord(
#                ra=_memory_reference_table["CATALOG_RA"],
#                dec=_memory_reference_table["CATALOG_DEC"],
#                unit="deg",
#            )
#            displacement = positions_reference.separation(positions_catalog).to_value(
#                "arcmin"
#            )
#
#            # -- compute the PSF value -- #
#            scores = 1 - psf(displacement)
#
#            # -- compute the exceptions -- #
#            exceptions = exc_call(_memory_reference_table[exc_c])
#
#            scores[exceptions == 1] = 0
#
#            output_reference_table = _memory_reference_table.loc[:, :]
#            output_reference_table[cls._output_column_name] = scores
#
#        output_reference_table.to_sql(
#            temp_table,
#            cross_match_database._sql_engine,
#            index=False,
#            if_exists="append",
#        )
#        pbar.update(len(_memory_reference_table))
#
#
# class PoissonReductionProcess(ReductionProcess):
#    pass
#
#
# class ObjectTypeReductionProcess(ReductionProcess):
#    pass
#
