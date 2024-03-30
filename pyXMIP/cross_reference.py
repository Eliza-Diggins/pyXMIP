"""
Cross-Referencing toolkit
"""
import os
import pathlib as pt

from _collections_abc import Collection
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from pyXMIP.schema import ReductionSchema
from pyXMIP.structures.databases import NED, SIMBAD
from pyXMIP.structures.table import SourceTable
from pyXMIP.utilities.core import mainlog

available_databases = [NED, SIMBAD]


# ============================================================================================ #
# X-Matching Processes                                                                         #
# ============================================================================================ #


def cross_match(
    input_path, output_path, databases=None, overwrite=False, *args, **kwargs
):
    mainlog.info(f"X-Matching {input_path} into {output_path}.")

    # ========================================================#
    # Managing arguments and kwargs. Enforcing types.
    input_path, output_path = pt.Path(input_path), pt.Path(output_path)
    source_table = SourceTable.read(input_path)

    if databases is None:
        databases = available_databases
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
        databases = available_databases
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
