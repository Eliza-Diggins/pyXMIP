import os
import pathlib as pt

from pyXMIP.structures.table import SourceTable


def table_answer_testing(table, filename, answer_store, answer_dir):
    filepath = pt.Path(os.path.join(answer_dir, filename))

    if answer_store:
        # store the data and don't actually test
        table.write(filepath, format="fits", overwrite=True)
    else:
        # check for the available.
        if not filepath.exists():
            table.write(filepath, format="fits")
        else:
            check_table = SourceTable.read(filepath, format="fits")

            print(check_table, table)
            assert len(check_table) == len(
                table
            ), f"Test table had len {len(table)} and check had len {len(check_table)}."
