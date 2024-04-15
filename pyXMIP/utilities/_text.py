"""Text utilities for the pyXMIP backend.

Notes
-----

The :py:mod:`utilities.text` module contains various utility functions for text output used elsewhere
in this package.
"""
import os
import pathlib as pt

from pyXMIP.utilities.core import xsparams

# -- loading the version dynamically from pyproject.toml -- #
_setup_tools_directory = os.path.join(pt.Path(__file__).parents[2], "setup.py")


def get_package_version():
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


def get_citation():
    return "not yet implemented"


_pxmip_ascii_logo = rf"""
                         __  ____  __ ___ ____
              _ __  _   _\ \/ /  \/  |_ _|  _ \
             | '_ \| | | |\  /| |\/| || || |_) |
             | |_) | |_| |/  \| |  | || ||  __/
             | .__/ \__, /_/\_\_|  |_|___|_|
             |_|    |___/

         ------------------------------------------

      Written by Eliza C. Diggins
          University of Utah Dept. of Physics and Astronomy
          (2020- )

      Version: {get_package_version()}

-------------------------------------------------------------------
"""


def print_cli_header():
    print(_pxmip_ascii_logo)


def show_config():
    try:
        from rich.pretty import pprint

        pprint(xsparams, expand_all=True)
    except ImportError:
        print(xsparams)


if __name__ == "__main__":
    print_cli_header()
