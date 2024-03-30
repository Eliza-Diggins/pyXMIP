"""
Core utilities with ubiquitous use cases in the ``pyXs`` package.
"""
import logging
import operator
import os
import pathlib as pt
import sys
from functools import reduce

import astropy.units as u
import sqlalchemy as sql
import yaml

# -- configuration directory -- #
_bin_directory = os.path.join(pt.Path(__file__).parents[1], "bin")
_config_directory = os.path.join(pt.Path(__file__).parents[1], "bin", "config.yaml")


# ======================================================================================================================#
# YAML loader custom definitions                                                                                       #
# ======================================================================================================================#
def _yaml_unit_constructor(loader: yaml.FullLoader, node: yaml.nodes.MappingNode):
    kw = loader.construct_mapping(node)
    i_s = kw["input_scalar"]
    del kw["input_scalar"]
    return i_s * u.Unit(kw["units"])


def _yaml_sql_type_constructor(loader: yaml.FullLoader, node: yaml.nodes.ScalarNode):
    return getattr(sql.types, loader.construct_scalar(node))


def _get_loader():
    loader = yaml.FullLoader
    loader.add_constructor("!unit", _yaml_unit_constructor)
    loader.add_constructor("!sql", _yaml_sql_type_constructor)
    return loader


# ======================================================================================================================#
# Configuration File                                                                                                    #
# ======================================================================================================================#
try:
    with open(_config_directory, "r+") as config_file:
        xsparams = yaml.load(config_file, _get_loader())

except FileNotFoundError as er:
    raise FileNotFoundError(
        f"Couldn't find the configuration file! Is it at {_config_directory}? Error = {er.__repr__()}"
    )
except yaml.YAMLError as er:
    raise yaml.YAMLError(
        f"The configuration file is corrupted! Error = {er.__repr__()}"
    )


# -- Configuring the loggers -- #
stream = (
    sys.stdout
    if xsparams["system"]["logging"]["main"]["stream"] in ["STDOUT", "stdout"]
    else sys.stderr
)
mainLogger = logging.getLogger("pyXs")

xs_sh = logging.StreamHandler(stream=stream)

# create formatter and add it to the handlers
formatter = logging.Formatter(xsparams["system"]["logging"]["main"]["format"])
xs_sh.setFormatter(formatter)
# add the handler to the logger
mainLogger.addHandler(xs_sh)
mainLogger.setLevel(xsparams["system"]["logging"]["main"]["level"])
mainLogger.propagate = False

mainlog = mainLogger

# -- Setting up the developer debugger -- #
devLogger = logging.getLogger("development_logger")

if xsparams["system"]["logging"]["developer"][
    "enabled"
]:  # --> We do want to use the development logger.
    # -- checking if the user has specified a directory -- #
    if xsparams["system"]["logging"]["developer"]["output_directory"] is not None:
        from datetime import datetime

        dv_fh = logging.FileHandler(
            os.path.join(
                xsparams["system"]["logging"]["developer"]["output_directory"],
                f"{datetime.now().strftime('%m-%d-%y_%H-%M-%S')}.log",
            )
        )

        # adding the formatter
        dv_formatter = logging.Formatter(
            xsparams["system"]["logging"]["main"]["format"]
        )

        dv_fh.setFormatter(dv_formatter)
        devLogger.addHandler(dv_fh)
        devLogger.setLevel("DEBUG")
        devLogger.propagate = False

    else:
        mainlog.warning(
            "User enabled development logger but did not specify output directory. Dev logger will not be used."
        )
else:
    devLogger.propagate = False
    devLogger.disabled = True


def enforce_units(value, preferred_units):
    if isinstance(value, u.Quantity):
        return value.to(preferred_units)
    else:
        return value * u.Unit(preferred_units)


def getFromDict(dataDict, mapList):
    return reduce(operator.getitem, mapList, dataDict)


def setInDict(dataDict, mapList, value):
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value
