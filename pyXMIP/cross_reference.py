"""
Cross-Referencing toolkit
"""
import operator
import os
import pathlib as pt
from abc import ABC, abstractmethod

import astropy.units as units
import numpy as np
import pandas as pd
import sqlalchemy as sql
from _collections_abc import Collection
from astropy.coordinates import SkyCoord
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from pyXMIP.schema import ReductionSchema
from pyXMIP.structures.databases import DEFAULT_DATABASE_REGISTRY
from pyXMIP.structures.table import SourceTable
from pyXMIP.utilities.core import getFromDict, mainlog, setInDict

# ============================================================================================ #
# X-Matching Processes                                                                         #
# ============================================================================================ #

_operators = {
    "<": operator.lt,
    "<=": operator.le,
    "=": operator.eq,
    ">=": operator.ge,
    ">": operator.gt,
}


def cross_match(
    input_path, output_path, databases=None, overwrite=False, *args, **kwargs
):
    mainlog.info(f"X-Matching {input_path} into {output_path}.")

    # ========================================================#
    # Managing arguments and kwargs. Enforcing types.
    input_path, output_path = pt.Path(input_path), pt.Path(output_path)
    source_table = SourceTable.read(input_path)

    if databases is None:
        databases = DEFAULT_DATABASE_REGISTRY.as_list()
    elif not isinstance(databases, Collection):
        databases = [databases]
    else:
        pass

    return cross_match_table(
        source_table,
        output_path,
        databases=databases,
        overwrite=overwrite,
        *args,
        **kwargs,
    )


def cross_match_table(
    table, output_path, databases=None, overwrite=False, *args, **kwargs
):
    import sqlalchemy as sql

    # =====================================================#
    # Managing args and kwargs.
    output_path = pt.Path(output_path)

    if databases is None:
        databases = DEFAULT_DATABASE_REGISTRY.as_list()
    mainlog.info(
        f"Cross matching with {len(databases)} databases: {[k.__name__ for k in databases]}."
    )

    # -- Dealing with the sql data -- #
    if output_path.exists():
        # check if the specific table is there
        _tmp_engine = sql.create_engine(f"sqlite:///{output_path}")
        _insp = sql.inspect(_tmp_engine)
        table_names = _insp.get_table_names()

        for table in table_names:
            if table in [f"{db.__name__}_MATCH" for db in databases]:
                if not overwrite:
                    raise ValueError(
                        f"Table {table} exists in {output_path} and overwrite=False."
                    )
                else:
                    mainlog.warning(
                        f"Table {table} exists in {output_path}. Overwrite = True -> deleting."
                    )

                    with _tmp_engine.connect() as conn:
                        _exec = sql.text(f"DROP TABLE '{table}'")
                        conn.execute(_exec)

    # ---------------------------------------- #
    # Cross Matching
    # ---------------------------------------- #
    with logging_redirect_tqdm(loggers=[mainlog]):
        for database in tqdm(
            databases, desc=f"Matching from {len(databases)} databases"
        ):
            database.source_match(output_path, table, *args, **kwargs)


# ================================================================================================ #
# REDUCTION MAIN FUNCTION                                                                          #
# ================================================================================================ #


class _REDUCTION_DESCRIPTOR:
    def __init__(
        self, dict_location, default=None, required=False, dtype=None, allow_set=True
    ):
        self.default = default
        self.dict_location = dict_location
        self.required = required
        self.dtype = dtype
        self.allow_set = allow_set

    def __set_name__(self, owner, name):
        self._name = name
        self.dict_location += [self._name]

    def __get__(self, instance, owner):
        try:
            return getFromDict(instance._kwargs, self.dict_location)
        except KeyError:
            return self.default

    def __set__(self, instance, value):
        setInDict(instance._kwargs, self.dict_location, value)


class ReductionProcess(ABC):
    """
    Abstract class representation for a standard cross-match reduction process.
    """

    _base_sql_string = (
        "SELECT * FROM %(table_name)s LIMIT %(chunksize)s OFFSET %(offset)s"
    )

    # === PARAMETERS === #
    CHUNKSIZE = _REDUCTION_DESCRIPTOR([], default=10000, dtype=int, allow_set=True)

    def __init__(self, database_path, databases="all", db_registry=None, **kwargs):
        self._kwargs = {}
        self._sql_string = (
            self.__class__._base_sql_string
        )  # Allow for custom sql string generation inside of process setup.
        self.included_databases = databases
        self.path = database_path

        for p in self.get_parameter_names():
            p_cls = self.__class__.__dict__[p]

            try:
                setattr(self, p, getFromDict(kwargs, p_cls.dict_location))
            except (KeyError, AttributeError):
                # that attribute doesn't have a value.
                if p_cls.required:
                    raise ValueError(f"The kwarg {p} is required.")

        for k, v in kwargs.items():
            if k not in self.get_parameter_names() and v is not None:
                mainlog.warning(
                    f"Detected kwarg {k}, which is not recognized. Check that it is spelled correctly."
                )

        # =================================== #
        # Managing args
        if db_registry is None:
            db_registry = DEFAULT_DATABASE_REGISTRY

        self._engine = sql.create_engine(f"sqlite:///{database_path}")
        # we need to lookup the available tables.
        with self._engine.connect() as conn:
            self._sql_tables = sql.inspect(conn).get_table_names()

        if self.included_databases == "all":
            self.included_databases = [
                db_registry[k.split("_MATCH")[0]]
                for k in self._sql_tables
                if "_MATCH" in k
            ]
        assert all(
            f"{k.__name__}_MATCH" in self._sql_tables for k in self.included_databases
        ), f"Databases {[k.__name__ for k in self.included_databases if f'{k}_MATCH' not in self._sql_tables]} were not found in the match database."

    def __repr__(self):
        return f"<ReductionProcess {self.__class__.__name__}>"

    def __str__(self):
        return self.__repr__()

    def __call__(self, *args, **kwargs):
        """
        Provides the core of all reduction processes.


        """
        mainlog.info(f"Performing {self} on {self.path}.")

        # --------------------------------------------------------------------------------- #
        # PROCESS SETUP
        # --------------------------------------------------------------------------------- #
        self._setup_reduction(*args, **kwargs)
        _progress_bar = tqdm(total=len(self.included_databases), desc="", leave=False)

        # --------------------------------------------------------------------------------- #
        # PROCESS MAIN
        # --------------------------------------------------------------------------------- #
        with logging_redirect_tqdm(loggers=[mainlog]):
            for (
                database
            ) in self.included_databases:  # iterate through each of the databases.
                # -- Setting up the reduction process -- #
                _progress_bar.set_description(
                    f"{self.__class__.__name__} - {database.__name__}"
                )
                _database_table_name = (
                    f"{database.__name__}_MATCH"  # Holds the actual match data.
                )
                _temp_database_table_name = f"{_database_table_name}_TMP"  # Where we are putting the full output.

                # getting the table size
                with self._engine.connect() as conn:
                    _table_size = conn.execute(
                        sql.text(f"SELECT COUNT(*) FROM {_database_table_name}")
                    ).scalar()

                offsets = np.arange(0, _table_size, self.CHUNKSIZE, dtype="uint32")
                _chunk_pbar = tqdm(
                    total=_table_size, desc="Chunking Execution", leave=False
                )

                # --------------------------------------------------------------------- #
                # REDUCTION
                # --------------------------------------------------------------------- #
                for offset in offsets:
                    _mem_table = pd.read_sql(
                        self.get_sql_query(
                            self.CHUNKSIZE, offset, _database_table_name
                        ),
                        self._engine,
                    )
                    self._process(_mem_table).to_sql(
                        _temp_database_table_name,
                        self._engine,
                        index=False,
                        if_exists="append",
                    )
                    _chunk_pbar.update(len(_mem_table))

                with self._engine.connect() as conn:
                    conn.execute(sql.text(f"DROP TABLE {_database_table_name}"))
                    conn.execute(
                        sql.text(
                            f"ALTER TABLE {_temp_database_table_name} RENAME TO {_database_table_name}"
                        )
                    )

                _progress_bar.update()
                mainlog.info(
                    f"{self.__class__.__name__} - {database.__name__} - [COMPLETE]"
                )

    def _setup_reduction(self, *args, **kwargs):
        pass

    @abstractmethod
    def _process(self, table):
        return table

    def get_sql_query(self, cs, off, db):
        return sql.text(f"SELECT * FROM {db} LIMIT {cs} OFFSET {off}")

    @property
    def parameters(self):
        return {m: getattr(self, m) for m in self._get_descriptors()}

    @classmethod
    def _get_descriptors(cls):
        # fetching descriptors.
        return [
            m for m, v in cls.__dict__.items() if isinstance(v, _REDUCTION_DESCRIPTOR)
        ]

    @classmethod
    def get_parameter_names(cls):
        return cls._get_descriptors()

    def _enforce_parameters(self):
        for k, v in self.parameters.items():
            if v is None and self.__class__.__dict__[k].required:
                raise ValueError(f"Parameter {k} is not set.")

    def _pull_table(self, table):
        with self._engine.connect() as conn:
            return pd.read_sql_table(table, conn)


class IdentityReductionProcess(ReductionProcess):
    def _process(self, table):
        return table


class SeparationReductionProcess(ReductionProcess):
    def _process(self, table):
        catalog_coords = SkyCoord(
            ra=table["CATALOG_RA"], dec=table["CATALOG_DEC"], unit="deg"
        )
        match_coords = SkyCoord(ra=table["RA"], dec=table["DEC"], unit="deg")

        table["SEPARATION"] = catalog_coords.separation(match_coords).to_value("arcmin")

        return table


class InstrumentReductionProcess(ReductionProcess):
    column_name = "InstrumentRedScore"

    PSF = _REDUCTION_DESCRIPTOR([], required=True, allow_set=True, dtype=units.Quantity)

    EXT_FLAG_COL = _REDUCTION_DESCRIPTOR([], required=False, allow_set=True, dtype=str)
    EXT_FLAG_THRESH = _REDUCTION_DESCRIPTOR(
        [], required=False, allow_set=True, dtype=(float, int, str)
    )
    EXT_FLAG_COMP_OP = _REDUCTION_DESCRIPTOR([], required=False, default="=")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # -- adding additional objects -- #
        self._ext_func = None
        self._psf_func = None

    def _setup_reduction(self, *args, **kwargs):
        # ==================================================== #
        # Check Columns are present
        with self._engine.connect() as conn:
            for table in [f"{db.__name__}_MATCH" for db in self.included_databases]:
                assert any(
                    [
                        k["name"] == "SEPARATION"
                        for k in sql.inspect(conn).get_columns(table)
                    ]
                ), f"Failed to find SEPARATION column in {table}. You should run a SeparationReductionProcess..."

            # Check that a catalog reduction is already there
            tables = sql.inspect(conn).get_table_names()

            if self.EXT_FLAG_COL is not None:
                assert (
                    "CATALOG" in tables
                ), f"Failed to find required table [because EXT_FLAG_COL is set] CATALOG in {self.path}."
                assert any(
                    _cname["name"] == self.EXT_FLAG_COL
                    for _cname in sql.inspect(conn).get_columns("CATALOG")
                ), f"Failed to find {self.EXT_FLAG_COL} in CATALOG."

        # ==================================================== #
        # Managing custom SQL pull information.
        if self.EXT_FLAG_COL is not None:
            # --> the sql call needs to be altered to allow for mergers.
            self._sql_string = f" %(table_name)s.*,CATALOG.{self.EXT_FLAG_COL} FROM %(table_name)s LIMIT %(chunksize)s OFFSET %(offset)s INNER JOIN "

        if self.EXT_FLAG_COL is not None:
            # set the special SQL call to include the EXT_FLAG_COL

            self._ext_func = lambda table: _operators[self.EXT_FLAG_COMP_OP](
                table[self.EXT_FLAG_COL], self.EXT_FLAG_THRESH
            )
        else:
            self._ext_func = lambda table: np.array([False] * len(table))

        # -- build the PSF function -- #
        self._psf_func = lambda table, sigma=self.PSF.to_value("arcmin"): np.exp(
            -0.5 * (table["SEPARATION"] / sigma) ** 2
        )

    def _process(self, table):
        table["IS_EXT"] = self._ext_func(table)
        table[self.__class__.column_name] = np.ones(len(table))
        table.loc[
            table["IS_EXT"] is False, self.__class__.column_name
        ] = self._psf_func(table.loc[table["IS_EXT"] is False, :])

        del table["IS_EXT"]
        if self.EXT_FLAG_COL is not None:
            del table[self.EXT_FLAG_COL]

        return table

    def get_sql_query(self, cs, off, db):
        if self.EXT_FLAG_COL is not None:
            # we need to generate a special one.
            return sql.text(
                f"WITH TMP AS (SELECT * FROM {db} LIMIT {cs} OFFSET {off}) SELECT TMP.*,CATALOG.{self.EXT_FLAG_COL} FROM TMP INNER JOIN CATALOG USING(CATALOG_OBJECT)"
            )
        else:
            return super().get_sql_query(cs, off, db)


def reduce_match_database(schema, overwrite=False):
    """
    Reduce the directed match database as characterized by the schema.

    Parameters
    ----------
    schema: str or :py:class:`pyXMIP.schema.ReductionSchema`
        The schema for the reduction process.
    overwrite: bool, optional
        If ``True``, the function will allow for an overwrite of the ``REDUCED`` table in the sql database.

    Returns
    -------
    None
    """
    import sqlalchemy as sql

    # ---------------------------------------------- #
    # Schema Processing and loading
    if isinstance(schema, str):
        # We need to load the schema from file.
        schema = ReductionSchema.from_file(schema)
    else:
        # the schema is already of type schema.
        pass

    # ---------------------------------------------- #
    # Setup
    mainlog.info(f"Reducing X-match database {schema.DBPATH}.")
    with logging_redirect_tqdm(loggers=[mainlog]):
        _process_progress_bar = tqdm(
            total=100,
            desc="X-Match Reduction - Setup",
            bar_format="{desc}: {percentage:3.0f}%|{bar}|",
        )

        # -- Check the database -- #
        _sql_engine = sql.create_engine(f"sqlite:///{schema.DBPATH}")
        assert os.path.exists(
            schema.DBPATH
        ), f"The DBPATH {schema.DBPATH} doesn't exist."

        _sql_tables = sql.inspect(_sql_engine).get_table_names()

        if "REDUCED_MATCH" in _sql_tables and not overwrite:
            raise ValueError(
                f"Database {schema.DBPATH} contains REDUCED_MATCH table already and overwrite = False."
            )
        elif "REDUCED_MATCH" in _sql_tables and overwrite:
            mainlog.warning(
                f"Found REDUCED_MATCH in {schema.DBPATH}. Overwrite = True --> Overwriting."
            )
            with _sql_engine.connect() as conn:
                conn.execute(sql.text("DROP TABLE REDUCED"))
        else:
            pass

        _database_sql_tables = [k for k in _sql_tables if "_MATCH" in k]
        _database_sql_names = [k.split("_MATCH")[0] for k in _database_sql_tables]

        mainlog.info(
            f"Located {len(_database_sql_names)} database tables: {_database_sql_tables}."
        )

        # -- check included databases -- #
        if schema.REFDBS == "all":
            schema.REFDBS = _database_sql_names
        else:
            schema.REFDBS = [db for db in schema.REFDBS if db in _database_sql_names]

        mainlog.debug(
            f"Reduction processing on {len(schema.REFDBS)} databases: {schema.REFDBS}."
        )

        _process_progress_bar.desc = "X-Match Reduction - "
        _process_progress_bar.update(1)


def add_catalog(catalog_path, database_path, identifier_column, overwrite=False):
    from astropy.table import Table

    engine = sql.create_engine(f"sqlite:///{database_path}")

    tbl = Table.read(catalog_path)

    with engine.connect() as conn:
        table_names = sql.inspect(conn).get_table_names()

    if "CATALOG" in table_names:
        if overwrite:
            mainlog.info(f"Overwriting existing CATALOG table in {database_path}.")
            with engine.connect() as conn:
                conn.execute(sql.text("DROP TABLE CATALOG"))

        else:
            raise ValueError(
                f"Found existing CATALOG table in {database_path} and overwrite=False."
            )
    else:
        pass

    tbl = tbl.to_pandas()
    tbl.rename(columns={identifier_column: "CATALOG_OBJECT"}, inplace=True)

    # -- fixing datatypes -- #
    # This must be done because strings might come back as bytestrings -> np.objects -> BLOB in sql.
    for col, dtype in tbl.dtypes.items():
        if dtype == object:
            tbl[col] = tbl.loc[:, col].astype("string")

    tbl.to_sql("CATALOG", engine, index=False)
