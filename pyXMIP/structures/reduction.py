"""
Reduction processes and related structures.

Notes
-----

For guidance on working with reduction processes, we recommend looking at the guide: :ref:`cross_Referencing`.

"""
from abc import ABC, abstractmethod
from typing import Callable

import numpy as np
import pandas as pd
from pydantic import BaseModel, model_validator

from pyXMIP.cross_reference import PydanticCMD
from pyXMIP.utilities.logging import mainlog
from pyXMIP.utilities.sql import chunk_sql_query_operation
from pyXMIP.utilities.types import ICRSCoordinateStdErrorSpecifier, TableColumn

# -- Importing SELF -- #
try:
    from typing import Self  # noqa
except ImportError:
    from typing_extensions import Self as Self  # noqa


class ReductionProcess(BaseModel, ABC):
    """
    Abstract class structure for reduction processes.
    """

    # ------------------------------ #
    # Schema Specifiers              #
    # ------------------------------ #
    settings_header: str = None
    """str: The settings header under which to find the parameters for this reduction process.

    When initializing a reduction process or set of reduction processes from a schema, the schema's underlying data dictionary
    will contain all of the relevant parameters under a header (i.e. ``ABSTRACT_PROCESS_PARAMS``). This attribute allows ``pyXMIP``
    to identify the :py:class:`ReductionProcess` associated with the settings and set everything up correctly.
    """
    settings_flag: str = None
    """str: The schema flag in ``RUNTIME_PARAMS`` which marks this process as enabled.
    """
    # ------------------------------ #
    # Process Specifiers             #
    # ------------------------------ #
    process_name: str
    """str: The name of this reduction process.
    """
    cross_match_database: PydanticCMD
    """:py:class:`cross_reference.CrossMatchDatabase`: The CMD this reduction process will run on.
    """
    table: str
    """str: The table name this reduction process will run on."""
    fill_unknown: bool = False
    """bool: If ``True``, then process will attempt to generate values for the process."""

    def __init__(self, **kwargs):
        """
        Initialize the reduction process.

        Parameters
        ----------
        process_name: str
            The name of the final :py:class:`ReductionProcess` instance.
        **kwargs:
            parameter and assigned value pairs corresponding to the parameters of this reduction process.

        Notes
        -----
        The initialization process simply assigns the ``process_name`` to ``self.name`` and adds the kwargs to the
        ``self._parameter_dictionary``, which is then accessed by the descriptors. Parameters are validated during instantiation.

        """
        super().__init__(**kwargs)

    def __str__(self):
        return f"<ReductionProcess: {self.process_name}>"

    def __repr__(self):
        return self.__str__()

    @property
    def debug_flag(self) -> str:
        """
        The string annotation used in debugging statements for this process.
        """
        return f"[{self.process_name}|{self.table}]"

    @property
    def proc_tmp_table(self) -> str:
        """str: The temp table name where process results are tabulated."""
        return f"{self.process_name}_TMP"

    @property
    def op_tmp_table(self) -> str:
        """str: The temp table where the newly joined outputs are stored while replace occurs."""
        return f"{self.table}_TMP"

    @property
    def obj_col(self) -> str:
        """str: The column of the ``table`` containing the object identifier."""
        return self.cross_match_database.schema.table_schema[self.table].NAME

    @property
    def score_col(self) -> str:
        """str: The column name storing the scores for the process."""
        return f"{self.process_name}_SCORE"

    @model_validator(mode="after")
    def validate(self) -> Self:
        """
        Custom validator for this reduction process.

        The ``validate`` method is a standard ``pydantic`` validator (runs after ``__init__`` and returns ``self``) which
        parses the provided attributes and utilizes the :py:attr:`ReductionProcess.cross_match_database` and :py:attr:`ReductionProcess.table`
        values to fill in other missing values. This method should only succeed if there is enough information to generate the operator and
        SQL query.
        """
        # ----------------------------------------------- #
        # Constructing Missing Settings                   #
        # ----------------------------------------------- #
        if self.fill_unknown:
            # The unknowns are going to be filled.
            mainlog.info(
                f"{self.debug_flag} Attempting to generate any missing attributes."
            )
            self._fill_unknown_parameters()

        # ----------------------------------------------- #
        # Validating the Model                            #
        # ----------------------------------------------- #
        mainlog.info(f"{self.debug_flag} Validating model parameters.")

        # -- default checks -- #
        assert (
            self.table in self.cross_match_database.match_tables
        ), f"{self.debug_flag} TABLE {self.table} not in {self.cross_match_database}."

        # Non-standard checks
        self._validate()

        return self

    @abstractmethod
    def _fill_unknown_parameters(self):
        pass

    @abstractmethod
    def _validate(self):
        return self

    @abstractmethod
    def generate_sql_query(self) -> str:
        """
        Generate the correct ``sql`` query for this reduction process.

        Returns
        -------
        str
            The output SQL query.
        """
        pass

    @abstractmethod
    def generate_operator(self) -> Callable[[pd.DataFrame], pd.DataFrame | pd.Series]:
        """
        Generate the correct ``operator`` function for this reduction process.

        Returns
        -------
        Callable[[pd.DataFrame],pd.DataFrame|pd.Series]
            The function which operates on the table.

        Notes
        -----
        The output operator is a function with signature ``func(table: pd.DataFrame) -> pd.DataFrame``. The output dataframe
        **must** have 3 columns (``"<REDUCTION_PROCESS_NAME>_SCORE"``) with ``float`` values corresponding to the respective scores
        of each entry in the input table as well as the ``CATOBJ`` column specifying the catalog object and the column identifying
        the potential match.
        """
        pass

    def __call__(self, *args, overwrite=False, **kwargs):
        import sqlalchemy as sql

        # ------------------------------------- #
        # Setup and validation                  #
        # ------------------------------------- #
        tqdm_kwargs = {"desc": self.debug_flag}
        _sql_engine = self.cross_match_database._sql_engine

        # -- check that we can actually perform this reduction -- #
        _check_status = self.cross_match_database.check_meta(
            self.process_name, self.table
        )

        if _check_status and overwrite:
            # The process has already been run, but overwrite = True.
            with _sql_engine.connect() as conn:
                mainlog.warning(
                    f"{self.debug_flag} Process has already been performed. Overwritting."
                )
                conn.execute(
                    sql.text(
                        f"ALTER TABLE `{self.table}` DROP COLUMN `{self.process_name}_SCORE`"
                    )
                )
        elif _check_status and not overwrite:
            raise ValueError(
                f"{self.debug_flag} The process has already been performed and overwrite=False."
            )
        else:
            pass

        # -- construct the query and the operator -- #
        # These call to the two abstract generator methods for the setup procedure.
        query, _operator = self.generate_sql_query(), self.generate_operator()

        # ------------------------------------- #
        # Create the callable operator          #
        # ------------------------------------- #
        # The chunked SQL operator is called to produce a version of OPERATOR which runs efficiently on
        # chunks.
        #
        # The new signature will be f(engine,query,output_table,*args,**kwargs)
        operator = chunk_sql_query_operation(**dict(tqdm_kwargs=tqdm_kwargs))(_operator)

        # -------------------------------------- #
        # Running the operator                   #
        # -------------------------------------- #
        operator(_sql_engine, query, f"{self.process_name}_TMP", *args, **kwargs)

        # -------------------------------------- #
        # Post-Processing                        #
        # -------------------------------------- #
        self._post_process()

    def _post_process(self):
        # ---------------------------------------------- #
        # Constructing queries                           #
        # ---------------------------------------------- #
        import sqlalchemy as sql

        JOIN_QUERY = sql.text(
            f"CREATE TABLE {self.op_tmp_table} AS "
            f"SELECT {self.table}.* , {self.proc_tmp_table}.`{self.score_col}` "
            f"FROM {self.table} LEFT JOIN {self.proc_tmp_table} "
            f"ON (({self.table}.CATOBJ = {self.proc_tmp_table}.CATOBJ) AND "
            f"({self.table}.`{self.obj_col}` = {self.proc_tmp_table}.`{self.obj_col}`))"
        )

        DROP_QUERY = sql.text(f"DROP TABLE {self.table}")
        RENAME_QUERY = sql.text(
            f"ALTER TABLE {self.op_tmp_table} RENAME TO {self.table}"
        )

        with self.cross_match_database._sql_engine.connect() as conn:
            mainlog.info(
                f"{self.debug_flag} Cleaning up SQL components... (Joining tables)"
            )
            conn.execute(JOIN_QUERY)
            mainlog.info(
                f"{self.debug_flag} Cleaning up SQL components... (Dropping tables)"
            )
            conn.execute(DROP_QUERY)
            mainlog.info(
                f"{self.debug_flag} Cleaning up SQL components... (Renaming tables)"
            )
            conn.execute(RENAME_QUERY)

        # -- Managing the META data -- #
        self.cross_match_database.meta_add(
            self.process_name, self.table, is_reduction=True
        )


class AstrometricReductionProcess(ReductionProcess):
    """
    Reduction process for refining match confidence on the basis of astrometric precision and separation
    between sources.
    """

    settings_header: str = "ASTROMETRY_PARAMS"
    settings_flag: str = "enable_astrometry"
    process_name: str = "ASTROMETRIC_REDUCTION"

    # ----------------------------------------------------------- #
    # Settings                                                    #
    # ----------------------------------------------------------- #
    CATALOG_ERR: ICRSCoordinateStdErrorSpecifier | None = None
    """:py:class:`utilities.types.CoordinateStdErrorSpecifier`: The error specifier for the catalog.
    """
    DATABASE_ERR: ICRSCoordinateStdErrorSpecifier | None = None
    """:py:class:`utilities.types.CoordinateStdErrorSpecifier`: The error specifier for the database.
    """
    prior: Callable[[], float] = lambda x: 1 / np.array(x["NMATCH"])

    def _get_displacement_op(
        self,
    ) -> Callable[[pd.DataFrame, ...], tuple[np.ndarray, np.ndarray]]:
        """
        Constructs the displacement operator for this reduction process
        """
        from astropy.coordinates import SkyCoord
        from astropy.units import Quantity

        # ---------------------------------------------------- #
        # Construct the names for important columns            #
        # ---------------------------------------------------- #
        # Necessary columns are the positions and units for each of the
        # coordinate positions.
        _cra, _cdec, _dra, _ddec = "CATRA", "CATDEC", "DBRA", "DBDEC"
        _cra_unit, _cdec_unit = "deg", "deg"
        _dra_unit, _ddec_unit = (
            self.cross_match_database.schema.table_schema[
                self.table
            ].column_map.RA.unit,
            self.cross_match_database.schema.table_schema[
                self.table
            ].column_map.DEC.unit,
        )

        # ---------------------------------------------------- #
        # Construct the operator                               #
        # ---------------------------------------------------- #
        def _coord_func(
            table_chunk: pd.DataFrame,
            cra=_cra,
            cdec=_cdec,
            dra=_dra,
            ddec=_ddec,
            cra_unit=_cra_unit,
            cdec_unit=_cdec_unit,
            dra_unit=_dra_unit,
            ddec_unit=_ddec_unit,
        ) -> tuple[np.ndarray, np.ndarray]:
            # -- generate the coordinate positions -- #
            cp = SkyCoord(
                ra=Quantity(table_chunk[cra], unit=cra_unit),
                dec=Quantity(table_chunk[cdec], unit=cdec_unit),
            )
            dp = SkyCoord(
                ra=Quantity(table_chunk[dra], unit=dra_unit),
                dec=Quantity(table_chunk[ddec], unit=ddec_unit),
            )

            # -- return the separations -- #
            return (cp.ra - dp.ra).to_value("arcsec"), (cp.dec - dp.dec).to_value(
                "arcsec"
            )

        return _coord_func

    def generate_operator(self) -> Callable[[pd.DataFrame], pd.DataFrame | pd.Series]:
        # ----------------------------------- #
        # Get the column separation operator  #
        # ----------------------------------- #
        # the displacement_function maps the table to the separation between RA and DEC coords.
        displacement_function = self._get_displacement_op()

        # ------------------------------------ #
        # Construct main operator              #
        # ------------------------------------ #
        # There are 3 separate formulae depending on the nature of the specified errors.
        if self.DATABASE_ERR is None or (self.DATABASE_ERR.mode is None):
            # The database is assumed to be ideal. Thus, we model error domain as a single gaussian.
            crau, cdecu = (
                self.CATALOG_ERR.lon_error.unit,
                self.CATALOG_ERR.lat_error.unit,
            )

            # -- Pull additional meta-data -- #
            def _func(
                table_chunk: pd.DataFrame,
                df: Callable[
                    [pd.DataFrame, ...], tuple[np.ndarray, np.ndarray]
                ] = displacement_function,
            ):
                # -- compute displacement vector -- #
                ra_dis, dec_dis = df(table_chunk)
                # -- compute error values -- #
                sra, sdec = (
                    (np.array(table_chunk["CATRA_ERR"]) * crau).to_value("arcsec"),
                    (np.array(table_chunk["CATDEC_ERR"]) * cdecu).to_value("arcsec"),
                )

                coef = ((2 * np.pi) ** (-1)) * (sra * sdec) ** (-1)
                table_chunk["bf"] = coef * np.exp(
                    (-0.5) * ((ra_dis / sra) ** 2 + (dec_dis / sdec) ** 2)
                )

                return table_chunk

        elif self.CATALOG_ERR is None or (self.CATALOG_ERR.mode is None):
            # The catalog is assumed to be ideal. Thus, we model error domain as a single gaussian.
            crau, cdecu = (
                self.DATABASE_ERR.lon_error.unit,
                self.DATABASE_ERR.lat_error.unit,
            )

            # -- Pull additional meta-data -- #
            def _func(table_chunk, df=displacement_function):
                # -- compute displacement vector -- #
                ra_dis, dec_dis = df(table_chunk)
                # -- compute error values -- #
                sra, sdec = (
                    (np.array(table_chunk["DBRA_ERR"]) * crau).to_value("arcsec"),
                    (np.array(table_chunk["DBDEC_ERR"]) * cdecu).to_value("arcsec"),
                )

                coef = ((2 * np.pi) ** (-1)) * (sra * sdec) ** (-1)
                table_chunk["bf"] = coef * np.exp(
                    (-0.5) * ((ra_dis / sra) ** 2 + (dec_dis / sdec) ** 2)
                )
                return table_chunk

        else:
            # -- Complex formalism -- #
            crau, cdecu = (
                self.CATALOG_ERR.lon_error.unit,
                self.CATALOG_ERR.lat_error.unit,
            )
            drau, ddecu = (
                self.DATABASE_ERR.lon_error.unit,
                self.DATABASE_ERR.lat_error.unit,
            )

            def _func(table_chunk, df=displacement_function):
                # -- compute displacement vector -- #
                ra_dis, dec_dis = df(table_chunk)

                # -- compute error values -- #
                cra, cdec, dra, ddec = (
                    (np.array(table_chunk["CATRA_ERR"]) * crau).to_value("arcsec"),
                    (np.array(table_chunk["CATDEC_ERR"]) * cdecu).to_value("arcsec"),
                    (np.array(table_chunk["DBRA_ERR"]) * drau).to_value("arcsec"),
                    (np.array(table_chunk["DBDEC_ERR"]) * ddecu).to_value("arcsec"),
                )
                iscra, iscdec, isdra, isddec = (
                    cra**-2,
                    cdec**-2,
                    dra**-2,
                    ddec**-2,
                )

                coef = 1 / (
                    (2 * np.pi)
                    * (
                        (cra * cdec * dra * ddec)
                        / (np.sqrt((cra**-2 + dra**-2) * (cdec**-2 + ddec**-2)))
                    )
                )

                exp = (-0.5) * (
                    ((iscra * isdra) / (iscra + isdra)) * ra_dis**2
                    + ((iscdec * isddec) / (iscdec + isddec)) * dec_dis**2
                )

                table_chunk["bf"] = coef * np.exp(exp)
                return table_chunk

        # ------------------------------------------------ #
        # Implement the prior                              #
        # ------------------------------------------------ #
        def _pfunc(table_chunk, f=_func, df=displacement_function, p=self.prior):
            table_chunk = f(table_chunk, df=df)

            table_chunk[self.score_col] = (
                1 + (1 - p(table_chunk)) / (p(table_chunk) * table_chunk["bf"])
            ) ** (-1)

            return table_chunk

        return _pfunc

    def generate_sql_query(self) -> str:
        # ------------------------------------------------ #
        # Adding fixed columns to the column list          #
        # ------------------------------------------------ #
        # This includes the object identifier columns for both db and CATALOG along with the
        # coordinate columns in ICRS.
        columns = {
            "CATOBJ": "CATALOG.`CATOBJ`",
            f"{self.obj_col}": f"{self.table}.`{self.obj_col}`",
            "CATRA": "CATALOG.`RA`",
            "CATDEC": "CATALOG.`DEC`",
            "DBRA": f"{self.table}.`RA`",  # ALL tables have RA; its a fixed column.
            "DBDEC": f"{self.table}.`DEC`",  # ALL tables have DEC; its a fixed column.
            "NMATCH": f"{self.table}.`CATNMATCH`",
        }
        # ------------------------------------------------ #
        # Adding un-fixed columns to the column list       #
        # ------------------------------------------------ #
        for tbl, error_struct in zip(
            [self.table, "CATALOG"], [self.DATABASE_ERR, self.CATALOG_ERR]
        ):
            # Add coordinate error columns where necessary.
            if error_struct is None or error_struct.mode is None:
                # The error structure is None, we can just skip.
                continue
            elif error_struct.mode == "circular":
                # The error structure is circular. We need to check if the error specifier is a value or a columns.
                if isinstance(error_struct.position_error, TableColumn):
                    # this is a column and needs to be added.
                    columns[
                        "DBRA_ERR" if tbl != "CATALOG" else "CATRA_ERR"
                    ] = f"{tbl}.`{error_struct.position_error.name}`"
                    columns[
                        "DBDEC_ERR" if tbl != "CATALOG" else "CATDEC_ERR"
                    ] = f"{tbl}.`{error_struct.position_error.name}`"
                else:
                    pass
            else:
                # the error structure is axial.
                if isinstance(error_struct.lat_error, TableColumn):
                    columns[
                        "DBDEC_ERR" if tbl != "CATALOG" else "CATDEC_ERR"
                    ] = f"{tbl}.`{error_struct.lat_error.name}`"
                if isinstance(error_struct.lon_error, TableColumn):
                    columns[
                        "DBRA_ERR" if tbl != "CATALOG" else "CATRA_ERR"
                    ] = f"{tbl}.`{error_struct.lon_error.name}`"

        # -- Now construct the SQL query. -- #
        column_string = ", ".join([f"{v} as `{k}`" for k, v in columns.items()])

        SQL_QUERY = (
            "SELECT %(columns)s FROM %(table)s LEFT JOIN CATALOG USING(CATOBJ)"
            % dict(columns=column_string, table=self.table)
        )
        return SQL_QUERY

    def _validate(self):
        # -- Check that at least one of the processes is non-redundant -- #
        if (self.DATABASE_ERR is None) and (self.CATALOG_ERR is None):
            raise ValueError(
                f"{self.debug_flag} Neither component of the match has specified error."
            )

        for err_spec in [self.DATABASE_ERR, self.CATALOG_ERR]:
            if err_spec is None:
                continue

            err_spec.check_empty()

    def _fill_unknown_parameters(self):
        error_specifiers = [self.DATABASE_ERR, self.CATALOG_ERR]
        schemas = [
            self.cross_match_database.schema.table_schema[self.table],
            self.cross_match_database.schema.catalog_schema,
        ]

        for error_specifier, schema in zip(error_specifiers, schemas):
            # -- Skip unspecified error specifiers -- #
            # This way, the user is always trusted to ignore a component if they want to.
            if error_specifier is None:
                continue

            # -- Seek error columns in the schema -- #
            if schema.column_map.RA_ERR is not None:
                mainlog.info(
                    f"{self.debug_flag} Filled parameter lon_err with {schema.column_map.RA_ERR.name}"
                )
                error_specifier.lon_error = schema.column_map.RA_ERR
            if schema.column_map.DEC_ERR is not None:
                mainlog.info(
                    f"{self.debug_flag} Filled parameter lat_err with {schema.column_map.DEC_ERR.name}"
                )
                error_specifier.lat_error = schema.column_map.DEC_ERR


if __name__ == "__main__":
    # load the cross match database from disk
    from pyXMIP.cross_reference import CrossMatchDatabase

    cmd = CrossMatchDatabase("../../docs/source/examples/data/cross_matched.db")
    cmd.schema.catalog_schema.column_map.RA_ERR = {"name": "RA_UPERR", "unit": "arcsec"}
    cmd.schema.catalog_schema.column_map.DEC_ERR = {
        "name": "DEC_UPERR",
        "unit": "arcsec",
    }
    proc = AstrometricReductionProcess(
        cross_match_database=cmd,
        table="NED_STD_MATCH",
        CATALOG_ERR={},
        fill_unknown=True,
    )
    print(proc.CATALOG_ERR)
