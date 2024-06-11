"""
Testing utilities for pyXMIP
"""
import os
import pathlib as pt

from astropy.table import Table


def check_astropy_table(table, answer_dir, subpath, answer_store):
    path = pt.Path(os.path.join(answer_dir, subpath))

    if not path.parents[0].exists():
        path.parents[0].mkdir(parents=True)

    if answer_store or not path.exists():
        table.write(path, format="ascii", overwrite=True)
        return None

    else:
        old_table = Table.read(path, format="ascii")

        assert (
            old_table.to_pandas() == table.to_pandas()
        ).all, f"There was a change in the table at {path}."
