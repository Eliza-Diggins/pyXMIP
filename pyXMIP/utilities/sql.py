"""
SQL interaction module for ``pyXMIP``.
"""
from typing import Callable, Collection

import pandas as pd
import sqlalchemy as sql
from tqdm.asyncio import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from pyXMIP.utilities.logging import mainlog


def _count_table(engine: sql.Engine, table: str) -> int:
    with engine.connect() as conn:
        return int(conn.execute(sql.text(f"SELECT COUNT(*) FROM {table}")).scalar())


def chunk_sql_query_operation(tqdm_kwargs: dict = None):
    """
    Meta-decorator function to perform operations on SQL queries in chunks instead of loading the entire database into
    memory at once.

    Parameters
    ----------
    tqdm_kwargs: dict
        The arguments to pass to tqdm during runtime.

    Returns
    -------
    Callable
        The resulting decorator.

    """

    if tqdm_kwargs is None:
        tqdm_kwargs = {}

    def _chunk_sql_query_operation(
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

        def base_wrapper(
            engine: sql.Engine,
            sql_query: str,
            otable_name: str,
            *args,
            **kwargs,
        ) -> None:
            """
            Perform a chunkwise operation on a particular ``sql_query`` and write the resulting output to disk.

            Parameters
            ----------
            engine: sql.Engine
                The ``sqlalchemy`` engine connection to the underlying database.
            sql_query: str
                The sql query to run.
            otable_name: str
                The name to give to the output table.
            chunksize: int, optional
                The maximum allowed chunksize for each operation.
            args:
                Additional arguments to pass to the underlying function.
            kwargs:
                Additional keyword arguments to pass to the underlying function.
            """
            chunksize = kwargs.pop("chunksize", 1000)

            with engine.connect() as conn, logging_redirect_tqdm(loggers=[mainlog]):
                # Connect to the SQL engine and redirect logging via tqdm.

                for chunk_id, sql_chunk in enumerate(
                    tqdm(
                        pd.read_sql(sql_query, conn, chunksize=chunksize),
                        **tqdm_kwargs,
                    )
                ):
                    # Read the SQL query in chunks of specified chunksize. Operate on those chunks with the specified function.
                    # -- perform the function -- #
                    chunk_image = function(*(sql_chunk, *args), **kwargs)

                    # -- Writing the output to the SQL database -- #
                    # If chunk_id == 0, then we need to create the table,
                    # otherwise, we just append to the existing table.

                    if chunk_id == 0:
                        # -- This database table doesn't exist yet. -- #

                        chunk_image.to_sql(
                            otable_name, conn, if_exists="replace", index=False
                        )
                    else:
                        # -- Write -- #
                        chunk_image.to_sql(
                            otable_name, conn, if_exists="append", index=False
                        )

                    conn.commit()

        return base_wrapper

    return _chunk_sql_query_operation


def chunk_sql_table_operation(
    columns: Collection[str] | str = "*",
    inplace: bool = False,
    tqdm_kwargs: dict = None,
):
    """
    Meta-decorator function to perform operations on SQL tables in chunks instead of loading the entire database into
    memory at once.

    Parameters
    ----------
    tqdm_kwargs: dict
        The arguments to pass to tqdm during runtime.
    columns: str or Collection[str], optional
        The columns to include in the query. By default, this is ``"*"`` which will cause the function to operate
        on the complete table with all columns.
    inplace: bool, optional
        If ``True``, then (instead of creating a table ``"{table}_PROCESSED"``, the original table (``table``) will be replaced
        with the output of this operation.

    Returns
    -------
    Callable
        The resulting decorator.

    """
    if tqdm_kwargs is None:
        tqdm_kwargs = {}

    def _chunk_sql_table_operation(
        function: Callable[[pd.DataFrame, ...], pd.DataFrame]
    ):
        """
        Decorator function to perform operations on SQL table in chunks instead of loading the entire database into
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

            This is a function with signature ``func(engine: sql.Engine, table: str,chunksize:int = 1000, *args, **kwargs)``.
            The ``engine`` must be the SQL engine to connect to. ``*args,**kwargs`` are passed directly to ``function``.

        """
        # The base_wrapper needs to be converted into a new wrapper with the desired (simplified) signature.
        # This means that the sql query must be constructed and the otable generated.

        def wrapper(engine: sql.Engine, table: str, *args, **kwargs):
            """
            Perform a chunkwise operation on a particular ``sql_query`` and write the resulting output to disk.

            Parameters
            ----------
            engine: sql.Engine
                The ``sqlalchemy`` engine connection to the underlying database.
            table: str
                The sql table to run on.
            chunksize: int, optional
                The maximum size of the operation chunks.
            args:
                Additional arguments to pass to the underlying function.
            kwargs:
                Additional keyword arguments to pass to the underlying function.
            """
            # -- Construct the basic wrapper from the general decorator -- #

            N = (_count_table(engine, table) // kwargs.get("chunksize", 1000)) + 1
            tqdm_kwargs["total"] = N

            # Now we can create the wrapper.
            base_wrapper = chunk_sql_query_operation(tqdm_kwargs=tqdm_kwargs)(function)

            # -- Constructing the SQL query for the specified table -- #
            sql_query = "SELECT %(columns)s FROM %(table)s" % dict(
                columns=(columns if isinstance(columns, str) else ", ".join(columns)),
                table=table,
            )
            otable_name = (
                f"{table}_PROCESSED"  # Standard prescription for the output table.
            )

            # -- Execute the base_wrapper -- #
            base_wrapper(*(engine, sql_query, otable_name, *args), **kwargs)

            # -- Cleanup and post-processing -- #
            if inplace:
                with engine.connect() as conn:
                    # -- replace the tables -- #
                    conn.execute(sql.text(f"DROP TABLE {table}"))
                    conn.execute(
                        sql.text(f"ALTER TABLE {otable_name} RENAME TO {table}")
                    )

        return wrapper

    return _chunk_sql_table_operation
