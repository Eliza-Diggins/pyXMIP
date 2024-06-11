"""Text utilities for the pyXMIP backend.

Notes
-----

The :py:mod:`utilities.text` module contains various utility functions for text output used elsewhere
in this package.
"""
import os
import pathlib as pt

from pyXMIP.utilities.core import bin_directory

# ------------------------------------------------------------------ #
# Globals                                                            #
# ------------------------------------------------------------------ #
_setup_tools_directory = os.path.join(pt.Path(__file__).parents[2], "setup.py")


# ------------------------------------------------------------------ #
# Simple utility functions                                           #
# ------------------------------------------------------------------ #
def get_package_version() -> str | None:
    """
    Determine the version of the currently installed version.
    """
    try:
        with open(_setup_tools_directory, "r") as setup_file:
            for line in setup_file:
                if line.strip().startswith("version="):
                    # Extract the version string (assuming it's in single or double quotes)
                    version = (
                        line.split("=")[1].strip().strip(",").strip(r"'").strip('"')
                    )
                    return version
    except FileNotFoundError:
        return None


def print_version() -> None:
    print(get_package_version())


def print_cli_header() -> None:
    """
    Print the CLI title.
    """
    with open(os.path.join(bin_directory, "txt", "asciilogo.txt"), "r") as fo:
        header = fo.read()

    print(header % dict(version=get_package_version()))
