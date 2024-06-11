"""
Cross referencing tools and methods for ``pyXMIP`` catalogs.

More details on cross referencing can be found in the user guide: :ref:`cross_Referencing`

Notes
-----

The :py:mod:`cross_reference` module is the core of the user-tools in ``pyXMIP``. This module allows you to cross reference tables,
run queries and reduce results.

"""
import pathlib as pt
import time
from typing import Annotated, Any, Callable, Collection, Self, TypeVar

import numpy as np
import pandas as pd
import sqlalchemy as sql
from astropy.coordinates import ICRS, SkyCoord
from tqdm.auto import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from pyXMIP.schema import CMDSchema
from pyXMIP.structures.databases import (
    DEFAULT_DATABASE_REGISTRY,
    DBRegistry,
    SourceDatabase,
)
from pyXMIP.structures.table import SourceTable
from pyXMIP.utilities.logging import mainlog
from pyXMIP.utilities.types import _CMDTypePydanticAnnotation

# ========================================================= #
# Type hinting constructs                                   #
# ========================================================= #
RProc = TypeVar("RProc")

# ========================================================= #
# Cross Matching Table Class                                #
# ========================================================= #


class CrossMatchDatabase:
    """
    Database of cross-matching results.

    The :py:class:`CrossMatchDatabase` class provides an easy-to-use API for accessing the results of a cross-matching
    procedure (as stored in a SQL database). Furthermore, :py:class:`CrossMatchDatabase` instances provide a variety of
    methods for sorting, processing, and performing other important tasks on the underlying data.

    Notes
    -----

    The :py:class:`CrossMatchDatabase` instance provides access to an underlying SQL database with the following 3 types
    of tables:

    1. **Catalog**: (``CATALOG``) This table contains the original catalog that was used for the cross-matching procedure.
    2. **Match Tables**: These are the results of cross matching the catalog against each of the relevant databases.
    3. **META**: The Meta table contains information on what procedures have already been performed on the database.


    """

    def __init__(
        self,
        path: str | pt.Path,
        schema: CMDSchema = None,
        overwrite_schema: bool = False,
    ):
        """
        Initialize a :py:class:`CrossMatchDatabase` instance.

        Parameters
        ----------
        path: str or :py:class:`pathlib.Path`
            The path to the underlying ``SQL`` database to load.
        schema: :py:class:`schema.CMDSchema`, optional
            The :py:class:`schema.CMDSchema` corresponding to this cross-matching database. If it is not specified, then
            the class will search for a schema file with the correct signature in the same directory as ``path``.

            .. hint::

                For a filename (as specified by ``path``), we expect to locate a schema with path ``path_schema.json``.

        overwrite_schema: bool, optional
            If ``True`` and a schema is specified (``schema != None``), then that schema will be written to disk in the
            standard position even if one already exists.

        """

        # ------------------------------------------------ #
        # Defining attributes                              #
        # ------------------------------------------------ #
        self.path: pt.Path = pt.Path(path)
        """:py:class:`pathlib.Path`: Path to the underlying SQL database.
        """
        assert self.path.exists(), f"The path {self.path.absolute()} doesn't exist."

        self._sql_engine = sql.create_engine(f"sqlite:///{self.path}")

        # -- Hidden attributes for property creation -- #
        self._tables = None

        # -- Seek schema -- #
        self._schema, self.has_schema, self.schema_path = self._get_schema(
            schema=schema, overwrite_schema=overwrite_schema
        )

    def _get_schema(self, schema: CMDSchema = None, overwrite_schema: bool = False):
        """

        Parameters
        ----------
        schema
        overwrite_schema

        Returns
        -------

        """
        _schema_path = pt.Path(str(self.path.with_suffix("")) + "_schema.json")
        if schema is not None:
            # The schema was specified, we don't need to do anything.
            schema: CMDSchema = schema
            has_schema: bool = True

            if (_schema_path.exists() and overwrite_schema) or (
                not _schema_path.exists()
            ):
                # We need to replace the schema path.
                schema.write(
                    _schema_path, overwrite=overwrite_schema, file_format="json"
                )

        else:
            # We need to search for one.
            if _schema_path.exists():
                # the schema file exists and we can load it.
                schema = CMDSchema.read(_schema_path, file_format="json")
                has_schema = True
            else:
                mainlog.warning(
                    f"Failed to load CMDSchema with CrossMatchDatabase instance {self}."
                )
                schema = None
                has_schema = False

        return schema, has_schema, _schema_path

    def __getitem__(self, item: str) -> pd.DataFrame:
        with self.connect() as conn:
            return pd.read_sql_table(item, conn)

    def __str__(self):
        return f"<CrossMatchDatabase @ {self.path}>"

    def __repr__(self):
        return self.__str__()

    def __contains__(self, item: str | SourceDatabase) -> bool:
        if isinstance(item, str):
            return f"{item}_MATCH" in self.match_tables
        else:
            return f"{item.name}_MATCH" in self.match_tables

    @property
    def schema(self) -> CMDSchema | None:
        if self.has_schema:
            return self._schema
        else:
            mainlog.warning(
                f"There is no schema associated with {self}. You may need to set / create one."
            )
            return None

    @schema.setter
    def schema(self, value: CMDSchema):
        self._schema, self.has_schema, _ = self._get_schema(schema=value)

    def overwrite_schema(self):
        self._schema, self.has_schema, self.schema_path = self._get_schema(
            schema=self.schema, overwrite_schema=True
        )

    @property
    def tables(self) -> list[str]:
        """list of str: List of available tables.

        See Also
        --------
        drop_table, match_tables
        """
        if self._tables is None:
            _insp = sql.inspect(self._sql_engine)
            table_names = _insp.get_table_names()
            self._tables = table_names
        return self._tables

    @property
    def meta(self) -> pd.DataFrame:
        """:py:class:`pd.DataFrame`: Dataframe representing the meta-processes that have operated on this instance.

        The ``META`` table (within the larger SQL database) is a repository of the "edit history" of the database. Effectively,
        different processes can check ``META`` to determine if certain tasks have already been performed or if critical
        pre-requisite tasks have actually occurred.

        See Also
        --------
        meta_reset, meta_remove, meta_add, check_meta, build_meta_table

        """
        if "META" in self.tables:
            pass
        else:
            self.build_meta_table(overwrite=True)

        with self._sql_engine.connect() as conn:
            return pd.read_sql_table("META", conn)

    @property
    def has_catalog(self) -> bool:
        """bool: Returns ``True`` if the catalog is already loaded.

        If it's loaded, the original source catalog ends up in the ``CATALOG`` table of the SQL database.

        See Also
        --------
        add_catalog, add_catalog_from_table

        """
        return "CATALOG" in self.tables

    @property
    def match_tables(self) -> list[str]:
        """list of str: List of the available match tables.

        .. hint::

            Match tables are tables in the underlying SQL database which are actually cross-matching database.

        See Also
        --------
        tables, meta

        """
        return [i for i in self.tables if i[-6:] == "_MATCH"]

    def _reset_attributes(self):
        self._tables = None

    def drop_table(self, table_name: str):
        """
        Drop the table from the ``SQL`` database.

        Parameters
        ----------
        table_name: str
            The name of the table to delete.

        Returns
        -------
        None

        See Also
        --------
        tables, match_tables, cross_match
        """
        assert table_name in self.tables, f"Table {table_name} doesn't exist."

        with self._sql_engine.connect() as conn:
            query = sql.text(f"DROP TABLE {table_name}")
            conn.execute(query)

        mainlog.info(f"DELETED table {table_name} from {self.path}.")

    def query(self, query: str) -> pd.DataFrame:
        """
        Query the SQL database.

        Parameters
        ----------
        query: str
            The SQLITE flavor SQL query for this database.

        """
        with self._sql_engine.connect() as conn:
            return pd.read_sql_query(sql.text(query), conn)

    def cross_match(
        self,
        databases: list[SourceDatabase | str],
        registry: DBRegistry = None,
        **kwargs,
    ):
        """
        Add a new cross-matching table to this object.

        Parameters
        ----------
        databases: list of str or :py:class:`structures.databases.SourceDatabase`
            The database to cross-reference against.
        registry: :py:class:`structures.databases.DBRegistry`, optional
            The database registry to use. If unspecified, then the default registry is used.

        """
        from pyXMIP.structures.table import SourceTable

        # ================================================= #
        # Setup: args / kwargs                              #
        # ================================================= #
        if not self.has_schema:
            mainlog.warning(
                f"{self} doesn't have a schema. The catalog will generate it's own (possible erroneous) schema."
            )

        if registry is None:
            registry = DEFAULT_DATABASE_REGISTRY

        for k, database in enumerate(databases):
            if isinstance(database, str):
                # this is a string database and needs to be looked up.
                try:
                    databases[k] = registry[database]
                except KeyError:
                    raise ValueError(
                        f"The database string {database} was not found in the registry. Have you added the database to your registry?"
                    )
            else:
                # This database is a class of its own.
                pass

        if any(database in self for database in databases):
            mainlog.warning(
                f"Databases {[database.name for database in databases if database in self]} are already in this cross-match database. Delete them and re-run to overwrite."
            )
            databases = [database for database in databases if database not in self]

        mainlog.info(f"Adding cross-match results from {databases} to {self.path}.")

        # -- load the catalog as a source table -- #
        assert (
            self.has_catalog
        ), f"There is no CATALOG object in {self}. Add a catalog to proceed."

        with self.connect() as conn:
            _catalog = SourceTable.from_pandas(pd.read_sql_table("CATALOG", conn))

        if self.has_schema:
            _catalog.schema = self.schema.catalog_schema
        else:
            _ = (
                _catalog.schema
            )  # This just base loads the schema to make sure it generates.

        # ================================================= #
        # Cross Referencing                                 #
        # ================================================= #
        cross_match_table(_catalog, self.path, databases=databases, **kwargs)

    def build_meta_table(self, overwrite: bool = False):
        """
        Generate a blank version of the ``META`` table in the database.

        Parameters
        ----------
        overwrite: bool
            If ``True``, overwrites will be allowed.

        See Also
        --------
        meta, check_meta, meta_add, meta_remove, meta_reset
        """
        # ========================================= #
        # Setting up the procedure                  #
        # ========================================= #
        if "META" in self.tables:
            if overwrite:
                with self._sql_engine.connect() as conn:
                    conn.execute(sql.text("DROP TABLE META"))
            else:
                raise ValueError(
                    "Failed to generate new META table because META already exists and overwrite=False."
                )
        else:
            pass

        # ======================================== #
        # Building the table                       #
        # ======================================== #
        tbl = {
            "PROCESS": ["META_GENERATED"],
            "TABLE": ["ALL"],
            "DATE_RUN": [time.asctime()],
            "REDUCTION": [False],
        }

        pd.DataFrame(tbl).to_sql(
            "META", self._sql_engine, if_exists="replace", index=False
        )

    def meta_add(self, process: str, table: str, is_reduction: bool = False):
        """
        Add a record to the ``META`` table.

        Parameters
        ----------
        process: str
            The process to add.
        table: str
            The table to add.
        is_reduction: bool
            Mark this additional process as a reduction process.

            .. note::

                If a process is marked as a reduction process, then it should produce a column in the applied table
                called ``<process_name>_SCORE``.

        See Also
        --------
        meta, check_meta, build_meta_table, meta_remove, meta_reset
        """
        mainlog.debug(f"Added {process} flag to {self.path} for {table}.")
        tbl = {
            "PROCESS": [process],
            "TABLE": [table],
            "DATE_RUN": [time.asctime()],
            "REDUCTION": [is_reduction],
        }

        pd.DataFrame(tbl).to_sql(
            "META", self._sql_engine, index=False, if_exists="append"
        )
        self._reset_attributes()

    def meta_remove(self, process: str, table: str):
        """
        Remove a record from the ``META`` table.

        Parameters
        ----------
        process: str
            The process to remove.
        table: str
            The table to remove.

        See Also
        --------
        meta, check_meta, meta_add, build_meta_table, meta_reset

        """
        _new_meta = self.meta.loc[
            (self.meta["PROCESS"] != process) & (self.meta["TABLE"] != table), :
        ]

        _new_meta.to_sql("META", self._sql_engine, if_exists="replace", index=False)
        self._reset_attributes()

    def meta_reset(self):
        """
        Reset the ``META`` table.

        See Also
        --------
        meta, check_meta, meta_add, meta_remove, build_meta_table
        """
        self.build_meta_table(overwrite=True)
        self._reset_attributes()

    def check_meta(self, process: str, table: str) -> bool:
        """
        Check the ``META`` table for a record.

        Parameters
        ----------
        process: str
            The process to check for.
        table: str
            The table to check for.

        Returns
        -------
        bool

        See Also
        --------
        meta, check_meta, meta_add, meta_remove, meta_reset
        """
        return (
            len(
                self.meta.loc[
                    (self.meta["PROCESS"] == process) & (self.meta["TABLE"] == table), :
                ]
            )
            != 0
        )

    def get_database(
        self, table_name: str, registry: DBRegistry = None
    ) -> SourceDatabase:
        """
        Given a ``table_name`` in this cross-matching database, determines the database from which it originated.

        Parameters
        ----------
        table_name: str
            The table in the SQL file for which to seek a matching database.
        registry: DBRegistry, optional
            The database registry to search for the databases. If ``None`` (default), then the default registry is used.

            .. important::

                If you have custom databases that were used to build the cross-matching database, you'll need to make sure
                they're in the relevant registry or this process will fail.

        Returns
        -------
        SourceDatabase
            The database matching the desired table.

        Raises
        ------
        ValueError
            Raised if the process fails to locate a database with a name matching the table name.
        """
        if registry is None:
            registry = DEFAULT_DATABASE_REGISTRY
        # -- First, check that the table is actually valid -- #
        assert (
            table_name in self.match_tables
        ), f"The table {table_name} is not among the tables of this database or is not a match table."

        # -- Pull the database name -- #
        db_name = table_name[:-6]  # Removes _MATCH from the end of the name string.

        # -- Search for a database -- #
        for name in registry.names:
            if name == db_name:
                return registry[name]

        raise ValueError(
            f"Failed to find a database matching {db_name} in the provided registry. Have you made sure all of the relevant databases are in the registry?"
        )

    @staticmethod
    def chunk_db_operation(
        flag: str,
        afunc: Callable[[Any, str, str], tuple[Collection[Any], dict]],
        allow_overwrite: bool = False,
        **meta_kwargs,
    ) -> Callable[[Callable], Callable]:
        """
        Meta-decorator for wrapping base processes and simplifying runtime context to simple operation on table.

        Parameters
        ----------
        flag: str
            The process flag to assign to the decorated process / method. This should be unique. It is what META uses to
            identify if a process has already been run or not.
        afunc: Callable
            A function (staticmethod) which takes an ``instance``, ``table`` and ``process_flag`` (along with ``**kwargs``) and returns
            ``args`` and ``kwargs`` which can then be passed directly to the decorated method.
        allow_overwrite: bool, optional
            If ``True``, then overwriting is permitted on this operation.
        """
        from pyXMIP.utilities.sql import chunk_sql_query_operation

        def _chunk_db_operation(function: Callable[[pd.DataFrame, ...], pd.DataFrame]):
            """
            Decorator function to perform operations on SQL queries in chunks instead of loading the entire database into
            memory at once.

            Parameters
            ----------
            function: Callable
                The operation to perform on each of the query chunks.

                The signature of ``function`` should be ``func(table: pd.DataFrame, *args, **kwargs)``. ``*args`` and ``**kwargs`` are
                then provided at runtime by the user. (See the notes on the returned function)

            Returns
            -------
            Callable
                The resulting decorator output.

                This is a function with signature ``func(engine: sql.Engine, sql_query: str, otable_name: str,chunksize:int = 1000, *args, **kwargs)``.
                The ``engine`` must be the SQL engine to connect to, the ``sql_query`` is the query to execute and the result is then operated on
                by the wrapped function. Finally, ``otable_name`` specifies what name to provide to the output table. ``*args,**kwargs`` are passed
                directly to ``function``.

            """

            def wrapper(
                self: Self,
                sql_query: str,
                otable_name: str,
                *args,
                **kwargs,
            ) -> None:
                """
                Perform a chunkwise operation on a particular ``sql_query`` and write the resulting output to disk.

                Parameters
                ----------
                self: CrossMatchDatabase
                    The instance to perform this operation.
                sql_query: str
                    The sql query to run.
                otable_name: str
                    The name to give to the output table.
                chunksize: int, optional
                    The maximum allowed chunksize for each operation.
                overwrite: bool, optional
                    If ``True``, then the process will proceed regardless of the META check.
                args:
                    Additional arguments to pass to the underlying function.
                kwargs:
                    Additional keyword arguments to pass to the underlying function.
                """
                # -- Managing the overwrite protocol -- #
                # -- Setup -- #
                # fix the function signature.
                f = lambda *a, **k: function(self, *a, **k)

                # add TQDM params to meta_kwargs if not specified.
                if "tqdm_kwargs" not in meta_kwargs:
                    meta_kwargs["tqdm_kwargs"] = {}

                meta_kwargs["tqdm_kwargs"]["desc"] = meta_kwargs["tqdm_kwargs"].get(
                    "desc", f"[{flag}]"
                )
                _base_wrapper = chunk_sql_query_operation(**meta_kwargs)(f)

                overwrite = kwargs.pop("overwrite", False)
                if overwrite and not allow_overwrite:
                    mainlog.warning(
                        f"[{flag}@{otable_name}] Operation doesn't permit overwriting. Overwrite=False."
                    )
                    overwrite = False

                # -- perform the argument generator process -- #
                # We run the argument fetching function for this process (flag) over unspecified tables (this is a general process).
                if afunc is not None:
                    _args, _kwargs = afunc(self, flag, otable_name)

                    args = list(args) + list(_args)
                    kwargs = {
                        **kwargs,
                        **{k: v for k, v in _kwargs.items() if k not in kwargs},
                    }
                else:
                    pass

                # -- perform operation checks -- #
                operation_check = self.check_meta(flag, otable_name)

                if operation_check and not overwrite:
                    mainlog.warning(
                        f"[{flag}@{otable_name}] Failed to execute because it is already in META and overwrite=False."
                    )
                    return None
                elif operation_check and overwrite:
                    mainlog.info(
                        f"[{flag}@{otable_name}] This process has already been run, but overwrite=True so it has been permitted."
                    )
                else:
                    pass

                engine = self._sql_engine

                _base_wrapper(*(engine, sql_query, otable_name, *args), **kwargs)

                self.meta_add(flag, otable_name)

            return wrapper

        return _chunk_db_operation

    @staticmethod
    def chunk_db_table_operation(
        flag: str,
        afunc: Callable[[Any, str, str], tuple[Collection[Any], dict]],
        allow_overwrite: bool = False,
        **meta_kwargs,
    ) -> Callable[[Callable], Callable]:
        """
        Meta-decorator for wrapping base processes and simplifying runtime context to simple operation on table.

        Parameters
        ----------
        flag: str
            The process flag to assign to the decorated process / method. This should be unique. It is what META uses to
            identify if a process has already been run or not.
        afunc: Callable
            A function (staticmethod) which takes an ``instance``, ``table`` and ``process_flag`` (along with ``**kwargs``) and returns
            ``args`` and ``kwargs`` which can then be passed directly to the decorated method.
        allow_overwrite: bool, optional
            If ``True``, then overwriting is permitted on this operation.
        """
        from pyXMIP.utilities.sql import chunk_sql_table_operation

        def _chunk_db_table_operation(
            function: Callable[[pd.DataFrame, ...], pd.DataFrame]
        ):
            """
            Decorator function to perform operations on SQL queries in chunks instead of loading the entire database into
            memory at once.

            Parameters
            ----------
            function: Callable
                The operation to perform on each of the query chunks.

                The signature of ``function`` should be ``func(table: pd.DataFrame, *args, **kwargs)``. ``*args`` and ``**kwargs`` are
                then provided at runtime by the user. (See the notes on the returned function)

            Returns
            -------
            Callable
                The resulting decorator output.

                This is a function with signature ``func(engine: sql.Engine, sql_query: str, otable_name: str,chunksize:int = 1000, *args, **kwargs)``.
                The ``engine`` must be the SQL engine to connect to, the ``sql_query`` is the query to execute and the result is then operated on
                by the wrapped function. Finally, ``otable_name`` specifies what name to provide to the output table. ``*args,**kwargs`` are passed
                directly to ``function``.

            """

            def wrapper(
                self: Self,
                table: str,
                *args,
                **kwargs,
            ) -> None:
                """
                Perform a chunkwise operation on a particular ``sql_query`` and write the resulting output to disk.

                Parameters
                ----------
                self: CrossMatchDatabase
                    The instance to perform this operation.
                chunksize: int, optional
                    The maximum allowed chunksize for each operation.
                overwrite: bool, optional
                    If ``True``, then the process will proceed regardless of the META check.
                args:
                    Additional arguments to pass to the underlying function.
                kwargs:
                    Additional keyword arguments to pass to the underlying function.
                """
                # -- Setup -- #
                # fix the function signature.
                f = lambda *a, **k: function(self, *a, **k)

                # add TQDM params to meta_kwargs if not specified.
                if "tqdm_kwargs" not in meta_kwargs:
                    meta_kwargs["tqdm_kwargs"] = {}
                meta_kwargs["tqdm_kwargs"]["desc"] = f"[{flag}@{table}]"

                _base_wrapper = chunk_sql_table_operation(**meta_kwargs)(f)

                # -- Managing the overwrite protocols -- #
                overwrite = kwargs.pop("overwrite", False)
                if overwrite and not allow_overwrite:
                    mainlog.warning(
                        f"[{flag}@{table}] Operation doesn't permit overwriting. Overwrite=False."
                    )
                    overwrite = False

                # -- perform the argument generator process -- #
                # We run the argument fetching function for this process (flag) over unspecified tables (this is a general process).
                if afunc is not None:
                    _args, _kwargs = afunc(self, flag, table)

                    args = list(_args) + list(args)
                    kwargs = {
                        **kwargs,
                        **{k: v for k, v in _kwargs.items() if k not in kwargs},
                    }

                else:
                    pass

                # -- perform operation checks -- #
                operation_check = self.check_meta(flag, table)

                if operation_check and not overwrite:
                    mainlog.warning(
                        f"[{flag}@{table}] Failed to execute because it is already in META and overwrite=False."
                    )
                    return None
                elif operation_check and overwrite:
                    mainlog.info(
                        f"[{flag}@{table}] This process has already been run, but overwrite=True so it has been permitted."
                    )
                else:
                    pass

                engine = self._sql_engine

                _base_wrapper(*(engine, table, *args), **kwargs)
                self.meta_add(flag, table)

            return wrapper

        return _chunk_db_table_operation

    @staticmethod
    def _args_gen_OBJECT_TYPES(
        instance: Any, process_flag: str, table: str
    ) -> tuple[list[Any], dict]:
        try:
            return [table, instance.schema.table_schema[table]], {}
        except KeyError as err:
            raise ValueError(
                f"[{process_flag}@{table}] Failed to fetch args / kwargs. MSG={err.__str__()}"
            )

    @chunk_db_table_operation(
        "CORRECT_OBJ_TYPES",
        afunc=_args_gen_OBJECT_TYPES,
        allow_overwrite=False,
        inplace=True,
    )
    def correct_object_types(self, table, *args, **kwargs):
        """
        Standardized process to correct the object type in a cross-matching table.

        Generically, when a :py:class:`CrossMatchDatabase` is generated, the match_tables are not corrected for replace their
        object types with the SIMBAD standard that ``pyXMIP`` recognizes. This process is designed to perform the conversion.

        Parameters
        ----------
        tables: list[str], optional
            The tables to perform the reduction process on. If this parameter is not specified, then all of the tables will
            be processed.
        overwrite: bool, optional
            If ``True``, then the process will be performed regardless of whether or not it has been performed before.

            .. warning::

                This will cause the process to be performed on a table which may already have the process run. This will
                often be either redundant or break your code. Use with caution.
        chunksize: int, optional
            The size of the chunks to pass each table through in. In effect, this controls how memory intensive the
            process is. For each chunk, the table is read in, processed, and written to a temp table in the SQL file. As such,
            smaller chunksizes will lead to longer read / write times but lower memory usage. By default, ``chunksize = 10000``.

        """
        _ = kwargs
        converted_chunk = SourceTable._convert_types(table, args[1])
        return converted_chunk

    @chunk_db_table_operation(
        "STNDIZE_COORDS",
        afunc=_args_gen_OBJECT_TYPES,
        allow_overwrite=False,
        inplace=True,
    )
    def standardize_coordinates(self, table, *args, **kwargs):
        """
        Ensure that all cross-matching tables provided (``tables``) have an RA and DEC column in them. This can then
        be used as the standard backbone coordinate system for all reductions.

        Parameters
        ----------
        tables: list[str], optional
            The tables to perform the reduction process on. If this parameter is not specified, then all of the tables will
            be processed.
        overwrite: bool, optional
            If ``True``, then the process will be performed regardless of whether or not it has been performed before.

            .. warning::

                This will cause the process to be performed on a table which may already have the process run. This will
                often be either redundant or break your code. Use with caution.
        chunksize: int, optional
            The size of the chunks to pass each table through in. In effect, this controls how memory intensive the
            process is. For each chunk, the table is read in, processed, and written to a temp table in the SQL file. As such,
            smaller chunksizes will lead to longer read / write times but lower memory usage. By default, ``chunksize = 10000``.

        """
        table_name, schema = args

        coordinates = SkyCoord(
            *[
                np.array(table[col.name]) * col.unit
                for col in schema.coordinate_columns
            ],
            frame=schema.coordinate_frame,
        )
        icrs_coords = coordinates.transform_to(ICRS)
        table["RA"], table["DEC"] = (
            icrs_coords.frame.spherical.lon.deg,
            icrs_coords.frame.spherical.lat.deg,
        )

        return table

    @chunk_db_table_operation(
        "SCORE",
        afunc=lambda instance, table, flag: ([], {}),
        allow_overwrite=True,
        inplace=True,
    )
    def _score(self, table_chunk, *args, **kwargs):
        _ = kwargs
        weight_dict = args[0]
        table_chunk["SCORE"] = np.sum(
            np.array(
                [
                    v * np.array(table_chunk[f"{k}_SCORE"])
                    for k, v in weight_dict.items()
                ]
            ),
            axis=0,
        )

        return table_chunk

    def add_catalog(self, catalog_path, schema=None, overwrite=False, **kwargs):
        r"""
        Add the catalog to the database.

        Parameters
        ----------
        catalog_path: str
            The path to the catalog data file.
        schema: :py:class:`schema.SourceTableSchema`, optional
            The schema to associate with the catalog. If a schema is not provided then the system will attempt to construct
            one.
        overwrite: bool, optional
            If ``True``, the process will run regardless of whether or not there is an existing CATALOG table.
        kwargs:
            Additional arguments to pass through the method.

        Returns
        -------
        None
        """
        # ========================================= #
        # Setting up the procedure                  #
        # ========================================= #

        # -- looking for the catalog path -- #
        catalog_path = pt.Path(catalog_path)
        assert catalog_path.exists(), f"The catalog at {catalog_path} doesn't exist."

        mainlog.debug(f"Adding {catalog_path} to {self.path}.")
        # Load the catalog into memory.
        catalog = SourceTable.read(catalog_path, format=kwargs.pop("format", None))

        if schema is not None:
            catalog.schema = schema

        self.add_catalog_from_table(catalog, overwrite=overwrite, **kwargs)

    def add_catalog_from_table(self, catalog, overwrite=False, **kwargs):
        """
        Add the catalog to the database.

        Parameters
        ----------
        catalog: :py:class:`structures.table.SourceTable`
            The source table to add.
        overwrite: bool, optional
            If ``True``, the process will run regardless of whether or not there is an existing CATALOG table.
        kwargs:
            Additional arguments to pass through the method.

        Returns
        -------
        None
        """

        _pname = "CATALOG_INCLUDED"
        # -- managing the schema -- #
        schema = catalog.schema

        # We need to assure that we have a NAME column.
        assert (
            schema.NAME is not None
        ), "The schema doesn't have a directive for the NAME column. Try manually providing a schema."
        mainlog.debug(
            f"Schema indicates {schema.NAME} is the object identifier. Renaming to CATOBJ."
        )
        # Coercing format
        # We attach the catalog as-is except for the name column.
        # The name column gets renamed to CATALOG_OBJECT

        catalog.rename_column(schema.NAME, "CATOBJ")
        catalog.remove_columns(kwargs.pop("ignore_columns", []))

        # ========================================= #
        # Run the procedure                         #
        # ========================================= #
        if self.has_catalog:
            if overwrite:
                mainlog.warning(
                    f"Process {_pname} has already been completed. Overwriting..."
                )
            else:
                mainlog.error(
                    f"Process {_pname} has already been completed. Skipping..."
                )
                return None

        # -- fixing datatypes -- #
        # This must be done because strings might come back as bytestrings -> np.objects -> BLOB in sql.
        catalog = catalog.to_pandas()

        for col, dtype in catalog.dtypes.items():
            if dtype == object:
                catalog[col] = catalog.loc[:, col].astype("string")

        # -- passing off to write -- #
        with self._sql_engine.connect() as conn:
            catalog.to_sql("CATALOG", conn, index=False, if_exists="replace")

        self.meta_add(_pname, "all")

    def _run_basic_corrections(self, catalog, tables=None, overwrite=False, **kwargs):
        """
        Run the basic baseline corrections.

        .. warning::

            This method is a standard part of the post-processing procedure. It is unlikely you would
            ever need to run this manually. Proceed with caution!

        Parameters
        ----------
        catalog: :py:class:`structures.table.SourceTable`
            The catalog used for the cross-matching.
        tables: list of str, optional
            The tables to perform this correction on. If unspecified, then all tables are corrected.
        registry: :py:class:`structures.databases.DBRegistry`, optional
            A database registry in which all of the listed databases are found. By default, this will be the
            built-in database registry.
        overwrite: bool, optional
            If ``True``, this will allow the process to run even if it has already been performed.

            .. warning::

                In most cases, this is unadvised because it will most likely fail.

        Returns
        -------
        None
        """
        # --------------------------------------- #
        # Setup processes and determine tables    #
        # --------------------------------------- #
        if tables is None:
            tables = self.match_tables

        # -- Adding the catalog to table -- #
        self.add_catalog_from_table(catalog, overwrite=overwrite, **kwargs)

        # -- Performing by-table corrections -- #
        for table in tables:
            self.correct_object_types(table, overwrite=overwrite, **kwargs)
            self.standardize_coordinates(table, overwrite=overwrite, **kwargs)

    def run_reduction(
        self, reduction_process: RProc, table: str, overwrite: bool = False
    ):
        """
        Run the provided reduction process on the specified table of this cross-matching database.

        Parameters
        ----------
        reduction_process: :py:class:`structures.reduction.ReductionProcess`
            The reduction process to perform on the specified ``table``.
        table: str
            The name of the table on which to perform the reduction process.
        overwrite: bool
            If the reduction process has already been performed on this table, then ``True`` will allow it to proceed while
            ``False`` will cause it to fail.

        """
        # -- check against the meta table -- #
        if self.check_meta(reduction_process.process_name, table):
            # The process has been performed before.
            if not overwrite:
                raise ValueError(
                    f"Reduction process {reduction_process} was already performed on {table} and overwrite = False."
                )
            else:
                mainlog.warning(
                    f"Reduction process {reduction_process} was already performed on {table}. Overwriting results..."
                )
                self.meta_remove(reduction_process.process_name, table)
                # remove from meta so we don't get duplicates later.

        # -- Run the reduction -- #
        reduction_process(self, table)  # --> runs the reduction.

        # -- post-processing -- #
        self.meta_add(
            reduction_process.process_name, table
        )  # Add the new reduction process to the table.

    def score_matches(
        self,
        weight: dict[tuple | str, float],
        tables: Collection[str] | None = None,
        show_summary: bool = False,
    ):
        r"""
        Score the :py:class:`cross_reference.CrossMatchDatabase` instance.

        Parameters
        ----------
        weight: dict of tuple: float or str: float
            The scoring weights for each of the reduction processes. Keys correspond to each of the reduction processes and
            the values correspond to process weight (:math:`\alpha_i \in [0,1]`). If keys are ``str``, then the weight is applied
            universally across tables. If the keys are ``tuple`` of the form ``('<reduction_process>','<table>')`` then the weights
            will be applied differently to different tables.
        tables: list of str, optional
            The tables to perform the scoring procedure on.

            This kwarg is only applied to weights specified with only the reduction process name. If weights are specified
            with both reduction process and table name, then this will overwrite the applied tables.
        show_summary: bool, optional
            If ``True``, then show the summary.
        """
        if not tables:
            tables = self.match_tables
        # ---------------------------------------------- #
        # Enforcing weight restrictions                  #
        # ---------------------------------------------- #
        # restructure the weight dictionary to contain tuples specifically.
        _altered_weights = {}

        for k, v in weight.items():
            if isinstance(k, str):
                _altered_weights = {
                    **_altered_weights,
                    **{(k, _sk): v for _sk in tables},
                }
            else:
                _altered_weights = {**_altered_weights, **{k: v}}

        weights = _altered_weights

        # re-weighting the weights.
        _total_tables = list(set([k[1] for k in _altered_weights]))

        for _tbl in _total_tables:
            _subweight_dict = {
                k: v for k, v in _altered_weights.items() if k[1] == _tbl
            }
            _weight_sum = np.sum(
                list(_subweight_dict.values())
            )  # compute the total weight value.
            # If the weights are not weighted correctly, we need to re-weight. *1/(sum of weights).
            for k, _ in _subweight_dict.items():
                weights[k] *= 1 / _weight_sum

        # ---------------------------------------------- #
        # Debugging summary                              #
        # ---------------------------------------------- #
        if show_summary:
            mainlog.debug(f"SCORING WEIGHT SUMMARY: {self}")
            for table_name in _total_tables:
                _subweights = {k: v for k, v in weights.items() if k[1] == table_name}
                mainlog.debug(f"\t[{table_name}] N reductions = {len(_subweights)}")

                for _k, _v in _subweights.items():
                    mainlog.debug(f"\t\t[{_k}] Weight = {_v}")

        # ----------------------------------------------- #
        # Performing the computations                     #
        # ----------------------------------------------- #
        # This should then be performed in chunked operations on each table.
        for table in tqdm(_total_tables):
            # cycle through all of the tables in the listing.

            # -- Construct the weight dictionary -- #
            _subweights = {k[0]: v for k, v in weights.items() if k[1] == table}

            self._score(table, _subweights)

    def plot_matches(
        self,
        catalog_object,
        table,
        fig=None,
        cmap=None,
        norm=None,
        vmin=None,
        vmax=None,
        resolution=300,
        fov="5 arcmin",
        hips_kwargs=None,
        scatter_kwargs=None,
        **kwargs,
    ):
        """
        Plot the matches to a given catalog source.

        Parameters
        ----------
        catalog_object: str
            The ``CATOBJ`` identifier for the object.
        table: str
            The table to search for matches in.
        fig: :py:class:`plt.axes.Figure`, optional
            The figure to add the axes to.
        cmap: str or :py:class:`plt.cm.ColorMap`
            The colormap to use. Default is viridis.
        norm: :py:class:`plt.colors.Norm`
            The norm to use for the colormap.
        vmin: float
            The minimum value of the background image.
        vmax: float
            The maximum value of the background image.
        resolution: int
            The number of pixels per side of the image (when using HIPS).
        fov: :py:class:`astropy.units.Quantity`
            The field of view for the HIPS image.
        hips_kwargs: dict
            Additional kwargs to pass through the hips generation system.

            +--------------------+-----------+---------------------------------------------------------------+
            | Name               | Type      | Description                                                   |
            +====================+===========+===============================================================+
            | ``enabled``        | `bool`    | [default: ``True``] Use a HIPS image as the background?       |
            +--------------------+-----------+---------------------------------------------------------------+


        scatter_kwargs: dict
            Additional kwargs to pass through scatter.
        **kwargs
            Additional key-word arguments as listed below.

            +--------------------+-----------+---------------------------------------------------------------+
            | Name               | Type      | Description                                                   |
            +====================+===========+===============================================================+
            | ``figsize``        | `tuple`   | [default: ``(10,10)``] The size of the figure to generate.    |
            +--------------------+-----------+---------------------------------------------------------------+


        Returns
        -------

        """
        import matplotlib.pyplot as plt
        from astropy.coordinates import SkyCoord
        from astropy.wcs import WCS

        from pyXMIP.utilities.plot import get_hips_image

        # ============================================================ #
        # Setting up runtime variables                                 #
        # ============================================================ #
        # -- kwargs / args management -- #
        hips_kwargs, scatter_kwargs = (
            {} if hips_kwargs is None else hips_kwargs,
            {} if scatter_kwargs is None else scatter_kwargs,
        )

        if fig is None:
            fig = plt.figure(figsize=kwargs.pop("figsize", (10, 10)))

        # -- Setup central object parameters -- #
        assert (
            self.has_catalog
        ), "Cannot plot matches without a loaded CATALOG table. Try adding the catalog."

        # -- Manage the central object information -- #
        object_data = self.query(
            f"SELECT * FROM CATALOG WHERE CATALOG_OBJECT == '{catalog_object}'"
        )
        RAC, DECC = object_data["RA"], object_data["DEC"]
        central_position = SkyCoord(ra=RAC, dec=DECC, unit="deg")

        # ============================================================ #
        # Pull object positions                                        #
        # ============================================================ #
        data_table = self.query(
            f"SELECT * FROM {table} WHERE CATOBJ == '{catalog_object}'"
        )

        # -- pull data -- #
        _ra, _dec = data_table["RA"], data_table["DEC"]

        # ============================================================ #
        # Producing the Image                                          #
        # ============================================================ #

        if hips_kwargs.pop("enabled", True):
            # We are using a HIPs map.
            hips_image = get_hips_image(
                central_position, fov, (resolution, resolution), **hips_kwargs
            )
            wcs = WCS(hips_image[0].header)
            ax = fig.add_subplot(111, projection=wcs)
            ax.imshow(hips_image[0].data, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax)
        else:
            from astropy.coordinates import Angle

            # We are not using a HIPs map, we need to perform this by hand.
            _fov = Angle(fov)

            _custom_wcs = {
                "CTYPE1": "RA---TAN",
                "CTYPE2": "DEC---TAN",
                "CUNIT1": "deg",
                "CUNIT2": "deg",
                "NAXIS1": resolution,
                "NAXIS2": resolution,
                "CRPIX1": resolution // 2,
                "CRPIX2": resolution // 2,
                "CRVAL1": RAC,
                "CRVAL2": DECC,
                "CDELT1": (1.2 * (2 * _fov) / resolution).to_value("deg"),
                "CDELT2": (1.2 * (2 * _fov) / resolution).to_value("deg"),
            }

            wcs = WCS(_custom_wcs)
            ax = fig.add_subplot(111, projection=wcs)

        # -- Manage the scatter plot -- #

        # Fetch the corrected parameters
        for _sra, _sdec in zip(_ra, _dec):
            ax.scatter(
                _sra,
                _sdec,
                transform=ax.get_transform("world"),
                **scatter_kwargs,
            )

        # adding source position scatter
        ax.scatter(
            central_position.ra,
            central_position.dec,
            transform=ax.get_transform("world"),
            s=(72 * fig.get_size_inches()[0] / 50) ** 2,
            color="red",
            marker="+",
        )
        # ============================================================ #
        # Returning outputs                                            #
        # ============================================================ #
        return ax, fig

    def connect(self):
        return self._sql_engine.connect()

    @classmethod
    def from_file(cls, path):
        """
        Load a :py:class:`CrossMatchDatabase` from file.

        Parameters
        ----------
        path: str
            The path from which to load the database.

        Returns
        -------
        :py:class:`CrossMatchDatabase`
            The returned matching database.

        """
        return cls(path)


PydanticCMD = Annotated[CrossMatchDatabase, _CMDTypePydanticAnnotation]


# ============================================================================================ #
# X-Matching Processes                                                                         #
# ============================================================================================ #
def cross_match(
    input_path: str | pt.Path,
    output_path: str | pt.Path,
    databases: list[SourceDatabase | str] = None,
    registry: DBRegistry = None,
    overwrite: bool = False,
    *args,
    **kwargs,
) -> CrossMatchDatabase:
    r"""
    Cross match a table of known objects against a set of databases and output the result to a path of your choice.

    Parameters
    ----------
    input_path: str
        The path to the input file (must be a table with a readable file format by ``astropy``).

        .. hint::

            This should be the catalog of sources which is to be cross-matched against the databases

    output_path: str, Path
        The ``.db`` file to write the output to. This will be a ``sqlite`` database containing all of the cross-matching
        output data.
    databases: list of :py:class:`structures.databases.SourceDatabase` or list of str, optional
        The databases to cross-match the catalog against. By default, this will be all of the databases currently loaded
        in the ``registry`` :py:class:`structures.databases.DBRegistry` instance provided.

        If entries are ``str``, then they will be looked up in the registry. If they are not, they will be taken as is.

    registry: :py:class:`structures.database.DBRegistry`
        The database registry to lookup databases from. If not provided, defaults to :py:attr:`structures.databases.DEFAULT_DATABASE_REGISTRY`.
    overwrite: bool, optional
        If ``True``, you will be allowed to overwrite a pre-existing ``.db`` file.

    Returns
    -------
    CrossMatchDatabase
        The output matching database.

    See Also
    --------
    :py:func:`cross_match_table`
    """
    mainlog.info(f"X-Matching {input_path} into {output_path}.")

    # ======================================================== #
    # Managing arguments and kwargs. Enforcing types.          #
    # ======================================================== #
    input_path, output_path = pt.Path(input_path), pt.Path(output_path)

    # check the source table reading.
    try:
        source_table = SourceTable.read(input_path)
    except Exception as excep:
        raise ValueError(
            f"Failed to read source table {input_path} because if error {excep.__str__()}."
        )

    # ======================================================== #
    # Running                                                  #
    # ======================================================== #
    return cross_match_table(
        source_table,
        output_path,
        databases=databases,
        registry=registry,
        overwrite=overwrite,
        *args,
        **kwargs,
    )


def cross_match_table(
    table: SourceTable,
    output_path: str | pt.Path,
    databases: list[SourceDatabase | str] = None,
    registry: DBRegistry = None,
    overwrite: bool = False,
    *args,
    **kwargs,
) -> CrossMatchDatabase:
    r"""
    Cross match a table of known objects against a set of databases and output the result to a path of your choice.

    Parameters
    ----------
    table: :py:class:`structures.table.SourceTable`
        The catalog (loaded into memory in a table) to cross-match against databases.
    output_path: str, Path
        The ``.db`` file to write the output to. This will be a ``sqlite`` database containing all of the cross-matching
        output data.
    databases: list of :py:class:`structures.databases.SourceDatabase` or str, optional
        The databases to cross-match the catalog against. By default, this will be all of the databases currently loaded
        in the ``registry`` :py:class:`structures.databases.DBRegistry` instance provided.

        If entries are ``str``, then they will be looked up in the registry. If they are not, they will be taken as is.

    registry: :py:class:`structures.database.DBRegistry`
        The database registry to lookup databases from. If not provided, defaults to :py:attr:`structures.databases.DEFAULT_DATABASE_REGISTRY`.
    overwrite: bool, optional
        If ``True``, you will be allowed to overwrite a pre-existing ``.db`` file.

    Returns
    -------
    CrossMatchDatabase
        The output matching database.

    See Also
    --------
    :py:func:`cross_match`
    """
    import sqlalchemy as sql

    # --------------------------------------- #
    # SETUP
    # --------------------------------------- #
    output_path = pt.Path(output_path)

    # configure the registry
    if registry is None:
        # we default to the standard database registry.
        registry = DEFAULT_DATABASE_REGISTRY

    # configure the databases
    if databases is None:
        databases = list(registry.values())
    else:
        # we need to check the databases.
        databases = [
            database if isinstance(database, SourceDatabase) else registry[database]
            for database in databases
        ]

    mainlog.info(
        f"Cross matching with {len(databases)} databases: {[db.name for db in databases]}."
    )

    # ===================================================== #
    # Setting up the SQL                                    #
    # ===================================================== #
    if output_path.exists():
        # check if the specific table is there
        _tmp_engine = sql.create_engine(f"sqlite:///{output_path}")
        _insp = sql.inspect(_tmp_engine)
        table_names = _insp.get_table_names()

        if "META" in table_names and overwrite:
            mainlog.warning(f"Table META of {output_path} exists. Deleting.")
            with _tmp_engine.connect() as conn:
                conn.execute(sql.text("DROP TABLE META"))

        for tbl in table_names:
            if tbl in [f"{db.name}_MATCH" for db in databases]:
                if not overwrite:
                    raise ValueError(
                        f"Table {tbl} exists in {output_path} and overwrite=False."
                    )
                else:
                    mainlog.warning(
                        f"Table {tbl} exists in {output_path}. Overwrite = True -> deleting."
                    )

                    with _tmp_engine.connect() as conn:
                        _exec = sql.text(f"DROP TABLE '{tbl}'")
                        conn.execute(_exec)

    # ===================================================== #
    # Running                                               #
    # ===================================================== #
    _schema_generator = {"table_schema": table.schema, "db_schema": {}}

    with logging_redirect_tqdm(loggers=[mainlog]):
        # noinspection PyTypeChecker
        for database in tqdm(databases, desc="Cross-Matching"):
            database.source_match(output_path, table, *args, **kwargs)

            # We need to add the database query schema to the schemas collection.
            _table_name = f"{database.name}_MATCH"
            _schema_generator["db_schema"][_table_name] = database.query_schema

    # ===================================================== #
    # Return                                                #
    # ===================================================== #
    mainlog.info("Post-processing the cross-matching database.")

    # -- create the CMDschema -- #
    schema = CMDSchema(
        table_schema=_schema_generator["db_schema"],
        catalog_schema=_schema_generator["table_schema"],
    )

    cmd = CrossMatchDatabase(output_path, schema=schema, overwrite_schema=True)
    cmd._run_basic_corrections(table)

    return CrossMatchDatabase(output_path)
