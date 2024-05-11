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

import pandas as pd
import sqlalchemy as sql
from astropy.coordinates import SkyCoord
from tqdm.auto import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from pyXMIP.structures.databases import (
    DEFAULT_DATABASE_REGISTRY,
    DBRegistry,
    SourceDatabase,
)
from pyXMIP.structures.table import SourceTable
from pyXMIP.utilities.core import mainlog


# ============================================================================================ #
# X-Matching Processes                                                                         #
# ============================================================================================ #
def cross_match(
    input_path,
    output_path,
    databases=None,
    registry=None,
    overwrite=False,
    *args,
    **kwargs,
):
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
    table, output_path, databases=None, registry=None, overwrite=False, *args, **kwargs
):
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
    with logging_redirect_tqdm(loggers=[mainlog]):
        for database in tqdm(
            databases, desc=f"Cross-Matching against {len(databases)} databases"
        ):
            database.source_match(output_path, table, *args, **kwargs)

    # ===================================================== #
    # Return                                                #
    # ===================================================== #
    mainlog.info("Post-processing the cross-matching database.")
    cmd = CrossMatchDatabase(output_path)

    if kwargs.pop("correct", True):
        _tmp_registry = DBRegistry(databases)
        cmd._run_basic_corrections(table, registry=_tmp_registry)
        del _tmp_registry

    return CrossMatchDatabase(output_path)


# ========================================================= #
# Cross Matching Table Class                                #
# ========================================================= #


class CrossMatchDatabase:
    """
    Class representation of a cross-matching database on disk. Allows access to matching tables, catalog data, and
    data reduction protocols through the ``pyXMIP`` interface.

    The :py:class:`CrossMatchDatabase` class provides methods for the following core functionalities:

    1. Interacting with the ``SQL`` database on disk representing all of the cross matching data.
    2. Add, remove, and edit tables within the ``SQL`` database in order to facilitate scientific inquiry.
    3. Allow the performance of **reduction** processes to obtain statistical information about the quality of a given match.

    To achieve these purposes, :py:class:`CrossMatchDatabase` has a number of sub-systems and methods which are summarized in
    this documentation.

    .. rubric:: Components

    :py:class:`CrossMatchDatabase` class instances represent ``SQL`` databases on disk with a few special components:

    - ``MATCH`` tables (named "DATABASE_MATCH" for each database)

      - These tables contain the raw cross-matching data pulled from each of the databases.

    - ``META`` table

      - Contains information about the processes which have been performed on this :py:class:`CrossMatchDatabase`.
      - This may include reduction processes, data post-processing, etc.

    - ``CATALOG`` table

      - The original catalog used to generate the database in the beginning. This allows the user to query additional
        databases and utilize data from the catalog table.


    """

    def __init__(self, path):
        """
        Initialize a :py:class:`CrossMatchDatabase` instance.

        Parameters
        ----------
        path: str or :py:class:`pathlib.Path`
            The path to the underlying ``SQL`` database to load.

        """

        # -- Defining attributes -- #
        self.path = path
        self._sql_engine = sql.create_engine(f"sqlite:///{self.path}")

        # -- Defining property backends -- #
        self._tables = None

    def __getitem__(self, item):
        with self.connect() as conn:
            return pd.read_sql_table(item, conn)

    def __str__(self):
        return f"<CrossMatchDatabase @ {self.path}>"

    def __repr__(self):
        return self.__str__()

    def __contains__(self, item):
        if isinstance(item, str):
            return f"{item}_MATCH" in self.match_tables
        else:
            return f"{item.name}_MATCH" in self.match_tables

    @property
    def tables(self):
        """list of str: List of available tables."""
        if self._tables is None:
            _insp = sql.inspect(self._sql_engine)
            table_names = _insp.get_table_names()
            self._tables = table_names
        return self._tables

    @property
    def meta(self):
        """:py:class:`pandas.DataFrame`: The ``META`` table."""
        if "META" in self.tables:
            pass
        else:
            self.build_meta_table(overwrite=True)

        with self._sql_engine.connect() as conn:
            return pd.read_sql_table("META", conn)

    @property
    def has_catalog(self):
        """bool: Returns ``True`` if the catalog is already loaded."""
        return "CATALOG" in self.tables

    @property
    def match_tables(self):
        """list of str: List of the available match tables."""
        return [i for i in self.tables if i[-6:] == "_MATCH"]

    def _reset_attributes(self):
        self._tables = None

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

    def cross_match(self, databases, schema=None, registry=None, **kwargs):
        """
        Add a new cross-matching table to this object.

        Parameters
        ----------
        databases: list of str or :py:class:`structures.databases.SourceDatabase`
            The database to cross-reference against.
        schema: :py:class:`schema.SourceTableSchema`
            The schema to associated internally with the ``CATALOG`` table.
            If ``None`` (default), then a schema will be deduced.
        registry: :py:class:`structures.databases.DBRegistry`, optional
            The database registry to use. If unspecified, then the default registry is used.
        """
        from pyXMIP.structures.table import SourceTable

        # ================================================= #
        # Setup: args / kwargs                              #
        # ================================================= #
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

        if schema is not None:
            _catalog.schema = schema
        else:
            _catalog.schema.NAME = "CATALOG_OBJECT"

        # ================================================= #
        # Cross Referencing                                 #
        # ================================================= #
        cross_match_table(_catalog, self.path, databases=databases, **kwargs)

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
        mainlog.debug(f"Added {process} flag to {self.path} for {table}.")
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

    def _check_meta_yield_database(
        self, process_flag, table, registry=None, overwrite=False
    ):
        # check that we can locate a database in the registry
        # check if the table already has a meta check for this.
        if self.check_meta(process_flag, table):
            if overwrite:
                mainlog.debug(
                    f"Process {process_flag} has already been completed on {table}. Attempting to overwrite [This will probably fail]."
                )
            else:
                mainlog.debug(
                    f"Process {process_flag} has already been completed on {table}. Skipping."
                )
                return None

        # check that we can locate a database in the registry
        try:
            database = registry[str(table).replace("_MATCH", "")]
        except KeyError:
            mainlog.error(
                f"Could not match database {str(table).replace('_MATCH', '')} to registry with databases {list(registry.keys())}. Skipping."
            )
            return None

        return database

    def correct_coordinates(
        self, tables=None, registry=None, overwrite=False, **kwargs
    ):
        """
        Correct the coordinate columns to assure that there are RA and DEC columns present.

        .. warning::

            This method is a standard part of the post-processing procedure. It is unlikely you would
            ever need to run this manually. Proceed with caution!

        Parameters
        ----------
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
        _process_flag_name = "CORRECT_COORDINATE_COLUMNS"
        # ========================================= #
        # Setting up the procedure                  #
        # ========================================= #
        # -- setting up kwargs and args -- #
        if registry is None:
            registry = DEFAULT_DATABASE_REGISTRY

        if tables is None:
            tables = self.match_tables

        # ========================================= #
        # Run the procedure                         #
        # ========================================= #
        for table in tables:
            # ---------------------------#
            # Table checks
            # ---------------------------#
            # check if the table already has a meta check for this.
            database = self._check_meta_yield_database(
                _process_flag_name, table, registry=registry, overwrite=overwrite
            )
            if database is None:
                continue
            # -- grab the schema -- #
            _db_schema = database.query_schema
            if ("RA" in _db_schema.column_map) and ("DEC" in _db_schema.column_map):
                # correction not needed; RA and DEC are present.
                self.meta_add(_process_flag_name, table)
                continue
            # -- Reduce fix columns -- #
            with self._sql_engine.connect() as conn:
                for val, table_chunk in enumerate(
                    tqdm(
                        pd.read_sql_table(
                            table, conn, chunksize=kwargs.get("chunksize", 10000)
                        ),
                        leave=False,
                    )
                ):
                    # determine the positions of the chunk #
                    positions = SkyCoord(
                        table_chunk[_db_schema.coordinate_columns[0]],
                        table_chunk[_db_schema.coordinate_columns[1]],
                        frame=_db_schema.coordinate_frame,
                        unit=_db_schema.default_angle_units,
                    )
                    new_positions = positions.transform_to("icrs")
                    converted_chunk = table_chunk.iloc[:, :]
                    converted_chunk["RA"] = new_positions.ra.deg
                    converted_chunk["DEC"] = new_positions.dec.deg
                    if val == 0:
                        converted_chunk.to_sql(
                            table + "_TMP",
                            conn,
                            if_exists="replace",
                            index=False,
                        )
                    else:
                        converted_chunk.to_sql(
                            table + "_TMP",
                            conn,
                            if_exists="append",
                            index=False,
                        )
                # -- replace the tables -- #
                conn.execute(sql.text(f"DROP TABLE {table}"))
                conn.execute(sql.text(f"ALTER TABLE {table}_TMP RENAME TO {table}"))
            # -- Add this process to meta -- #
            self.meta_add(_process_flag_name, table)

    def correct_column_names(
        self, tables=None, registry=None, overwrite=False, **kwargs
    ):
        """
        Correct the names of table columns to standardize convention.

        .. warning::

            This method is a standard part of the post-processing procedure. It is unlikely you would
            ever need to run this manually. Proceed with caution!

        Parameters
        ----------
        tables: list of str, optional
            The tables to perform this correction on. If unspecified, then all tables are corrected.
        registry: :py:class:`structures.databases.DBRegistry`, optional
            A database registry in which all of the listed databases are found. By default, this will be the
            built-in database registry.
        overwrite: bool, optional
            If ``True``, this will allow the process to run even if it has already been performed.

            .. warning::

                In most cases, this is unadvised because it will most likely fail.
        remove: list of str, optional
            If desired, these are the column names to remove from the table.

        Returns
        -------
        None

        """
        _ = kwargs
        _process_flag_name = "CORRECT_COLUMN_NAMES"
        # ========================================= #
        # Setting up the procedure                  #
        # ========================================= #
        # -- setting up kwargs and args -- #
        if registry is None:
            registry = DEFAULT_DATABASE_REGISTRY

        if tables is None:
            tables = self.match_tables

        # ========================================= #
        # Run the procedure                         #
        # ========================================= #
        for table in tables:
            # ---------------------------#
            # Table checks
            # ---------------------------#
            # check if the table already has a meta check for this.
            database = self._check_meta_yield_database(
                _process_flag_name, table, registry=registry, overwrite=overwrite
            )
            if database is None:
                continue
            # -- grab the schema -- #
            _db_schema = database.query_schema
            # -- Reduce to only particular columns -- #
            _column_rename_map = {
                "RA": "RA",
                "DEC": "DEC",
                "NAME": _db_schema.NAME,
                "TYPE": _db_schema.TYPE,
            }
            with self._sql_engine.connect() as conn:
                for k, v in _column_rename_map.items():
                    _query = sql.text(
                        f"ALTER TABLE '{database.name}_MATCH' RENAME '{v}' to '{k}'"
                    )
                    conn.execute(_query)
            # -- Add this process to meta -- #
            self.meta_add(_process_flag_name, table)

    def correct_object_types(
        self, tables=None, registry=None, overwrite=False, **kwargs
    ):
        """
        Correct the object types in tables.

        .. warning::

            This method is a standard part of the post-processing procedure. It is unlikely you would
            ever need to run this manually. Proceed with caution!

        Parameters
        ----------
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
        _ = kwargs
        _process_flag_name = "CORRECT_OBJECT_TYPES"
        # ========================================= #
        # Setting up the procedure                  #
        # ========================================= #
        if registry is None:
            registry = DEFAULT_DATABASE_REGISTRY

        if tables is None:
            tables = self.match_tables

        for table in tables:
            # -- setup and check everything is okay -- #
            # check if the table already has a meta check for this.
            database = self._check_meta_yield_database(
                _process_flag_name, table, registry=registry, overwrite=overwrite
            )

            if database is None:
                continue

            # -- grab the schema -- #
            _db_schema = database.query_schema

            if _db_schema.TYPE is None:
                mainlog.error(
                    f"Cannot correct object types for {table} because schema doesn't have TYPE column."
                )
                continue

            # -- Reduce to only particular columns -- #
            with self._sql_engine.connect() as conn:
                for val, table_chunk in enumerate(
                    tqdm(pd.read_sql_table(table, conn, chunksize=10000), leave=False)
                ):
                    converted_chunk = SourceTable._convert_types(
                        table_chunk, _db_schema
                    )

                    if val == 0:
                        converted_chunk.to_sql(
                            table + "_TMP",
                            conn,
                            if_exists="replace",
                            index=False,
                        )
                    else:
                        converted_chunk.to_sql(
                            table + "_TMP",
                            conn,
                            if_exists="append",
                            index=False,
                        )

                    # -- replace the tables -- #
                    conn.execute(sql.text(f"DROP TABLE {table}"))
                    conn.execute(sql.text(f"ALTER TABLE {table}_TMP RENAME TO {table}"))

                self.meta_add(_process_flag_name, table)

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

    def _run_basic_corrections(
        self, catalog, tables=None, registry=None, overwrite=False, **kwargs
    ):
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

        self.add_catalog_from_table(catalog, overwrite=overwrite, **kwargs)
        self.correct_coordinates(
            tables=tables, registry=registry, overwrite=overwrite, **kwargs
        )
        self.correct_object_types(
            tables=tables, registry=registry, overwrite=overwrite, **kwargs
        )
        self.correct_column_names(
            tables=tables, registry=registry, overwrite=overwrite, **kwargs
        )

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
