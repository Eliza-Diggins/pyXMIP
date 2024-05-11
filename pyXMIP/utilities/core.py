"""
Core utilities with ubiquitous use cases in the ``pyXs`` package.
"""
import logging
import operator
import os
import pathlib as pt
import sys
from contextlib import contextmanager
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


# ======================================================================================================================#
# Logging                                                                                                               #
# ======================================================================================================================#
class PyxmLogger(logging.Logger):
    """custom logging class with customizable verbosity.

    Verbosity levels:

    0 - development logging basically everything.
    1 - high verbosity: lots of info.
    2 - normal output.
    3 - only critical info.

    """

    def __init__(self, name, level=logging.NOTSET, verbosity=0):
        super().__init__(name, level)
        self.verbosity = verbosity
        self._fixed_verb = dict(
            debug=None, info=None, warning=None, error=None, critical=None
        )

    def debug(self, msg, *args, verb=0, **kwargs):
        if self._fixed_verb["debug"] is not None:
            verb = self._fixed_verb["debug"]

        if verb >= self.verbosity:
            super().debug(msg, *args, **kwargs)

    def warning(self, msg, *args, verb=0, **kwargs):
        if self._fixed_verb["warning"] is not None:
            verb = self._fixed_verb["warning"]

        if verb >= self.verbosity:
            super().warning(msg, *args, **kwargs)

    def info(self, msg, *args, verb=0, **kwargs):
        if self._fixed_verb["info"] is not None:
            verb = self._fixed_verb["info"]

        if verb >= self.verbosity:
            super().info(msg, *args, **kwargs)

    def error(self, msg, *args, verb=0, **kwargs):
        if self._fixed_verb["error"] is not None:
            verb = self._fixed_verb["error"]

        if verb >= self.verbosity:
            super().error(msg, *args, **kwargs)

    def critical(self, msg, *args, verb=0, **kwargs):
        if self._fixed_verb["critical"] is not None:
            verb = self._fixed_verb["critical"]

        if verb >= self.verbosity:
            super().critical(msg, *args, **kwargs)

    @contextmanager
    def fixed_verb(self, **kwargs):
        old_fix_verb = {**self._fixed_verb}
        for k, v in self._fixed_verb.items():
            self._fixed_verb[k] = kwargs.get(k, v)

        yield None

        self._fixed_verb = old_fix_verb

    def reset_fixed_verb(self):
        self._fixed_verb = dict(
            debug=None, warning=None, info=None, error=None, critical=None
        )


# -- Setup the logger -- #
stream = (
    sys.stdout
    if xsparams["system"]["logging"]["main"]["stream"] in ["STDOUT", "stdout"]
    else sys.stderr
)
mainLogger = PyxmLogger(
    "pyXMIP", verbosity=xsparams["system"]["logging"]["main"]["verbosity"]
)

xs_sh = logging.StreamHandler(stream=stream)

# create formatter and add it to the handlers
formatter = logging.Formatter(xsparams["system"]["logging"]["main"]["format"])
xs_sh.setFormatter(formatter)
# add the handler to the logger
mainLogger.addHandler(xs_sh)
mainLogger.setLevel(xsparams["system"]["logging"]["main"]["level"])
mainLogger.propagate = False

mainlog = mainLogger


def enforce_units(value, preferred_units):
    """
    Return a version of ``value`` with units of the type preferred or return
    an error if that isn't possible.

    Parameters
    ----------
    value: any
        The value to enforce units on. Can be array or scalar with numerical type with / without units.
    preferred_units: :py:class:`astropy.units.Unit` or str
        The unit to enforce.

    Returns
    -------
    :py:class:`astropy.units.Quantity`
        The output quantity with the preferred units.

    """
    if isinstance(value, u.Quantity):
        return value.to(preferred_units)
    else:
        return value * u.Unit(preferred_units)


def getFromDict(dataDict, mapList):
    return reduce(operator.getitem, mapList, dataDict)


def setInDict(dataDict, mapList, value):
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value
