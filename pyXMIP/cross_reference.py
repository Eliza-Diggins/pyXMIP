"""
``pyXMIP`` cross-referencing toolkit module.

Notes
-----

The :py:mod:`cross_reference` module is the core of the user-tools in ``pyXMIP``. This module allows you to cross reference tables,
run queries and reduce results.
"""
import pathlib as pt
import time

import pandas as pd
import sqlalchemy as sql
from _collections_abc import Collection
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from pyXMIP.structures.databases import DEFAULT_DATABASE_REGISTRY
from pyXMIP.structures.table import SourceTable
from pyXMIP.utilities.core import mainlog

# ============================================================================================ #
# X-Matching Processes                                                                         #
# ============================================================================================ #


def cross_match(
    input_path, output_path, databases=None, overwrite=False, *args, **kwargs
):
    r"""
    Cross match a table of known objects against a set of databases and output the result to a path of your choice.

    Parameters
    ----------
    input_path: str
        The path to the input file (must be a table with a readable file format by ``astropy``).
    output_path: str
        The ``.db`` file to write the output to.
    databases: list, optional
        The databases to include in the cross-matching process. If ``str`` types are listed, then they are assumed to
        be database names, otherwise it is assumed that the database classes are passed. By default, all databases are queried against.
    overwrite: bool, optional
        If ``True``, you will be allowed to overwrite a pre-existing ``.db`` file.
    *args:
        Additional arguments.
    **kwargs:
        Additional kwargs.

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
        overwrite=overwrite,
        *args,
        **kwargs,
    )


def cross_match_table(
    table, output_path, databases=None, overwrite=False, *args, **kwargs
):
    r"""
    Cross match a table of known objects against a set of databases and output the result to a path of your choice.

    Parameters
    ----------
    table: :py:class:`astropy.table.table.Table`
        The ``Table`` instance to cross-reference against.
    output_path: str
        The ``.db`` file to write the output to.
    databases: list, optional
        The databases to include in the cross-matching process. If ``str`` types are listed, then they are assumed to
        be database names, otherwise it is assumed that the database classes are passed. By default, all databases are queried against.
    overwrite: bool, optional
        If ``True``, you will be allowed to overwrite a pre-existing ``.db`` file.
    *args:
        Additional arguments.
    **kwargs:
        Additional kwargs.
    """
    import sqlalchemy as sql

    # ===================================================== #
    # Managing arguments and kwargs                         #
    # ===================================================== #
    output_path = pt.Path(output_path)

    # check if the databases are legit.
    if databases is None:
        # databases become the defaults.
        databases = DEFAULT_DATABASE_REGISTRY.as_list()
    elif not isinstance(databases, Collection):
        databases = [databases]
    else:
        pass

    # make sure they're actually database types.
    for k, v in enumerate(databases):
        if isinstance(v, str):
            databases[k] = DEFAULT_DATABASE_REGISTRY[v]

    mainlog.info(
        f"Cross matching with {len(databases)} databases: {[k.__name__ for k in databases]}."
    )

    # ===================================================== #
    # Setting up the SQL                                    #
    # ===================================================== #
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

    # ===================================================== #
    # Running                                               #
    # ===================================================== #
    with logging_redirect_tqdm(loggers=[mainlog]):
        for database in tqdm(
            databases, desc=f"Matching from {len(databases)} databases"
        ):
            database.source_match(output_path, table, *args, **kwargs)


# ========================================================= #
# Cross Matching Table Class                                #
# ========================================================= #


class CrossMatchDatabase:
    """
    Class wrapper for the underlying ``SQL`` database files where cross-matching data is contained.
    """

    _ppp = ["META_EXISTS", "COLUMN_CORRECT", "OBJECT_CORRECT", "CATALOG_INCLUDED"]
    _resetable_attributes = ["tables", "meta"]

    def __init__(self, path):
        """
        Initialize a :py:class:`CrossMatchDatabase` instance.

        Parameters
        ----------
        path: str
            The path to the underlying ``SQL`` database to load.

        """

        # -- Defining attributes -- #
        self.path = path
        self._sql_engine = sql.create_engine(f"sqlite:///{self.path}")

        # -- Defining property backends -- #
        self._tables = None
        self._meta = None

    @property
    def tables(self):
        """List of available tables."""
        if self._tables is None:
            _insp = sql.inspect(self._sql_engine)
            table_names = _insp.get_table_names()
            self._tables = table_names
        return self._tables

    @property
    def meta(self):
        """The ``META`` table."""
        if self._meta is None:
            with self._sql_engine.connect() as conn:
                self._meta = pd.read_sql_table("META", conn)
        return self._meta

    @property
    def has_catalog(self):
        """Returns ``True`` if the catalog is already loaded."""
        return "CATALOG" in self.tables

    @property
    def match_tables(self):
        """List of the available match tables."""
        return [i for i in self.tables if i[-6:] == "_MATCH"]

    def _reset_attributes(self):
        for u in self.__class__._resetable_attributes:
            setattr(self, f"_{u}", None)

    def drop_table(self, table_name):
        """
        Drop the table from the ``SQL`` database.

        Parameters
        ----------
        table_name: str
            The name of the table to delete.

        Returns
        -------
        None
        """
        assert table_name in self.tables, f"Table {table_name} doesn't exist."

        with self._sql_engine.connect() as conn:
            query = sql.text(f"DROP TABLE {table_name}")
            conn.execute(query)

        mainlog.info(f"DELETED table {table_name} from {self.path}.")

    def query(self, query):
        """
        Query the SQL database.

        Parameters
        ----------
        query: str
            The SQLITE flavor SQL query for this database.

        """
        with self._sql_engine.connect() as conn:
            return pd.read_sql_query(sql.text(query), conn)

    def build_meta_table(self, overwrite=False):
        """
        Generate the post-processing meta table ``META`` inside of the database.

        Parameters
        ----------
        overwrite: bool
            If ``True``, overwrites will be allowed.
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
        }

        pd.DataFrame(tbl).to_sql(
            "META", self._sql_engine, if_exists="replace", index=False
        )
        self.meta_add("META_EXISTS", "ALL")

    def meta_add(self, process, table):
        """
        Add a record to the ``META`` table.

        Parameters
        ----------
        process: str
            The process to add.
        table: str
            The table to add.

        """
        tbl = {
            "PROCESS": [process],
            "TABLE": [table],
            "DATE_RUN": [time.asctime()],
        }

        pd.DataFrame(tbl).to_sql(
            "META", self._sql_engine, index=False, if_exists="append"
        )
        self._reset_attributes()

    def meta_remove(self, process, table):
        """
        Remove a record from the ``META`` table.

        Parameters
        ----------
        process: str
            The process to remove.
        table: str
            The table to remove.

        """
        _new_meta = self.meta.loc[
            (self.meta["PROCESS"] != process) & (self.meta["TABLE"] != table), :
        ]

        _new_meta.to_sql("META", self._sql_engine, if_exists="replace", index=False)
        self._reset_attributes()

    def meta_reset(self):
        """
        Reset the ``META`` table.
        """
        self.build_meta_table(overwrite=True)
        self._reset_attributes()

    def check_meta(self, process, table):
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

        """
        return (
            len(
                self.meta.loc[
                    (self.meta["PROCESS"] == process) & (self.meta["TABLE"] == table), :
                ]
            )
            != 0
        )

    def correct_column_names(self, databases=None, registry=None, overwrite=False):
        """
        Correct the names of table columns to standardize convention.

        .. warning::

            This method is a standard part of the post-processing procedure. It is unlikely you would
            ever need to run this manually. Proceed with caution!

        Parameters
        ----------
        databases: list, optional
            The databases for which this process should be performed.
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
        _pname = "COLUMN_CORRECT"
        # ========================================= #
        # Setting up the procedure                  #
        # ========================================= #
        # -- setting up kwargs and args -- #
        if registry is None:
            registry = DEFAULT_DATABASE_REGISTRY

        if databases is None:
            databases = [
                registry[name.replace("_MATCH", "")] for name in self.match_tables
            ]

        # ========================================= #
        # Run the procedure                         #
        # ========================================= #
        with logging_redirect_tqdm(loggers=[mainlog]):
            for database in tqdm(databases, desc="Fixing column convention."):
                # -- setup and check everything is okay -- #
                table_name = f"{database.__name__}_MATCH"

                if self.check_meta(_pname, table_name):
                    if overwrite:
                        mainlog.warning(
                            f"Process {_pname} has already been completed on {table_name}. Attempting to overwrite [This will probably fail]."
                        )
                    else:
                        mainlog.error(
                            f"Process {_pname} has already been completed on {table_name}. Skipping."
                        )
                        continue

                # -- grab the schema -- #
                _db_schema = database.class_table_schema

                # -- Reduce to only particular columns -- #
                _column_rename_map = {
                    "CATALOG_OBJECT": "CATALOG_OBJECT",
                    "CATALOG_RA": "CATALOG_RA",
                    "CATALOG_DEC": "CATALOG_DEC",
                    "RA": "RA",
                    "DEC": "DEC",
                    "NAME": _db_schema.NAME,
                    "TYPE": _db_schema.TYPE,
                }
                with self._sql_engine.connect() as conn:
                    for k, v in _column_rename_map.items():
                        _query = sql.text(
                            f"ALTER TABLE '{database.__name__}_MATCH' RENAME '{v}' to '{k}'"
                        )
                        conn.execute(_query)

                # -- Add this process to meta -- #
                self.meta_add(_pname, table_name)

    def correct_object_types(self, databases=None, registry=None, overwrite=False):
        """
        Correct the object types in tables.

        .. warning::

            This method is a standard part of the post-processing procedure. It is unlikely you would
            ever need to run this manually. Proceed with caution!

        Parameters
        ----------
        databases: list, optional
            The databases for which this process should be performed.
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
        _pname = "OBJECT_CORRECT"
        # ========================================= #
        # Setting up the procedure                  #
        # ========================================= #
        if registry is None:
            registry = DEFAULT_DATABASE_REGISTRY

        if databases is None:
            databases = [
                registry[name.replace("_MATCH", "")] for name in self.match_tables
            ]

        with logging_redirect_tqdm(loggers=[mainlog]):
            for database in tqdm(
                databases, desc="Converting object types to standards."
            ):
                # -- setup and check everything is okay -- #
                table_name = f"{database.__name__}_MATCH"
                if self.check_meta(_pname, table_name):
                    if overwrite:
                        mainlog.warning(
                            f"Process {_pname} has already been completed on {table_name}. Attempting to overwrite [This will probably fail]."
                        )
                    else:
                        mainlog.error(
                            f"Process {_pname} has already been completed on {table_name}. Skipping."
                        )
                        continue

                # -- grab the schema -- #
                _db_schema = database.class_table_schema

                # -- Reduce to only particular columns -- #
                with self._sql_engine.connect() as conn:
                    for val, table_chunk in enumerate(
                        tqdm(pd.read_sql_table(table_name, conn, chunksize=10000))
                    ):
                        converted_chunk = SourceTable._convert_types(
                            table_chunk, _db_schema
                        )

                        if val == 0:
                            converted_chunk.to_sql(
                                table_name + "_TMP",
                                conn,
                                if_exists="replace",
                                index=False,
                            )
                        else:
                            converted_chunk.to_sql(
                                table_name + "_TMP",
                                conn,
                                if_exists="append",
                                index=False,
                            )

                    # -- replace the tables -- #
                    conn.execute(sql.text(f"DROP TABLE {table_name}"))
                    conn.execute(
                        sql.text(f"ALTER TABLE {table_name}_TMP RENAME TO {table_name}")
                    )

                self.meta_add(_pname, table_name)

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
        _pname = "CATALOG_INCLUDED"
        # ========================================= #
        # Setting up the procedure                  #
        # ========================================= #

        # -- looking for the catalog path -- #
        catalog_path = pt.Path(catalog_path)
        assert catalog_path.exists(), f"The catalog at {catalog_path} doesn't exist."

        mainlog.info(f"Adding {catalog_path} to {self.path}.")
        # Load the catalog into memory.
        catalog = SourceTable.read(catalog_path, format=kwargs.pop("format", None))

        # -- managing the schema -- #
        if schema is None:
            mainlog.warning("No schema provided for catalog. Auto-generating.")
            schema = catalog.schema

        # We need to assure that we have a NAME column.
        assert (
            schema.NAME is not None
        ), "The schema doesn't have a directive for the NAME column. Try manually providing a schema."
        mainlog.info(
            f"Schema indicates {schema.NAME} is the object identifier. Renaming to CATALOG_OBJECT."
        )
        # Coercing format
        # We attach the catalog as-is except for the name column.
        # The name column gets renamed to CATALOG_OBJECT

        catalog.rename_column(schema.NAME, "CATALOG_OBJECT")
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

        with self._sql_engine.connect() as conn:
            catalog.to_pandas().to_sql(
                "CATALOG", conn, index=False, if_exists="replace"
            )

        self.meta_add(_pname, "all")

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
