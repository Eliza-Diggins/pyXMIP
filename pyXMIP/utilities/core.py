"""
Core utilities for many tasks in the ``pyXMIP`` package environment.

Notes
-----

Aside from access to the logging system, this module has very little practical use for the end user.

"""
import functools
import operator
import os
import pathlib as pt
from functools import reduce
from typing import Any, Collection, Iterable, Mapping

import astropy.units as u
import ruamel.yaml
from sqlalchemy.types import TypeEngine

# -- Importing SELF -- #
try:
    from typing import Self  # noqa
except ImportError:
    from typing_extensions import Self as Self  # noqa

# -- configuration directory -- #
bin_directory: pt.Path = pt.Path(os.path.join(pt.Path(__file__).parents[1], "bin"))
# :py:class:`pathlib.Path`: The directory in which the ``pyXMIP`` bin is located.
config_directory: pt.Path = pt.Path(
    os.path.join(pt.Path(__file__).parents[1], "bin", "config.yaml")
)
# :py:class:`pathlib.Path`: The directory in which the ``pyXMIP`` configuration files are located.

yaml = ruamel.yaml.YAML()
yaml.register_class(u.Quantity)
yaml.register_class(TypeEngine)


class AttrDict(dict):
    """
    Attribute accessible dictionary.
    """

    def __init__(self, mapping: Mapping):
        super(AttrDict, self).__init__(mapping)
        self.__dict__ = self

        for key in self.keys():
            self[key] = self.__class__.from_nested_dict(self[key])

    @classmethod
    def from_nested_dict(cls, data: Any) -> Self:
        """Construct nested AttrDicts from nested dictionaries."""
        if not isinstance(data, dict):
            return data
        else:
            return AttrDict({key: cls.from_nested_dict(data[key]) for key in data})

    @classmethod
    def clean_types(cls, _mapping):
        for k, v in _mapping.items():
            if isinstance(v, AttrDict):
                _mapping[k] = cls.clean_types(_mapping[k])
            else:
                pass
        return dict(_mapping)

    def clean(self):
        return self.clean_types(self)


class YAMLConfiguration:
    """
    General class representing a YAML configuration file.
    """

    def __init__(self, path: pt.Path | str):
        # ------------------------------------------------- #
        # Defining parameters                               #
        # ------------------------------------------------- #
        self.path: pt.Path = pt.Path(path)
        # :py:class:`pathlib.Path`: The path to the underlying yaml file.
        self._config: ruamel.yaml.CommentedMap | None = None

    @property
    def config(self):
        if self._config is None:
            self._config = self.load()

        return AttrDict(self._config)

    @classmethod
    def load_from_path(cls, path: pt.Path) -> dict:
        """Read the configuration dictionary from disk."""
        try:
            with open(path, "r+") as cf:
                return yaml.load(cf)

        except FileNotFoundError as er:
            raise FileNotFoundError(
                f"Couldn't find the configuration file! Is it at {config_directory}? Error = {er.__repr__()}"
            )

    def load(self) -> dict:
        return self.__class__.load_from_path(self.path)

    def reload(self):
        """Reload the configuration file from disk."""
        self._config = None

    @classmethod
    def set_on_disk(cls, path: pt.Path | str, name: str | Collection[str], value: Any):
        _old = cls.load_from_path(path)

        if isinstance(name, str):
            name = name.split(".")
        else:
            pass

        setInDict(_old, name, value)

        with open(path, "w") as cf:
            yaml.dump(_old, cf)

    def set_param(self, name: str | Collection[str], value):
        self.__class__.set_on_disk(self.path, name, value)


pxconfig = YAMLConfiguration(config_directory)


def enforce_units(value: Any, preferred_units: u.Unit | str) -> u.Quantity:
    """
    Return a version of ``value`` with units of the type preferred or raise an error if that isn't possible.

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
        try:
            return value.to(preferred_units)
        except u.UnitConversionError as error:
            raise ValueError(
                f"Failed to enforce units {preferred_units} on {value}. MSG: {error.__str__()}"
            )
    else:
        return value * u.Unit(preferred_units)


def getFromDict(dataDict: Mapping, mapList: Iterable[slice]) -> Any:
    """
    Fetch an object from a nested dictionary using a list of keys.

    Parameters
    ----------
    dataDict: dict
        The data dictionary to search.
    mapList: list
        The list of keys to follow.

    Returns
    -------
    Any
        The output value.

    """
    return reduce(operator.getitem, mapList, dataDict)


def setInDict(dataDict: Mapping, mapList: Iterable[slice], value: Any):
    """
    Set the value of an object from a nested dictionary using a list of keys.

    Parameters
    ----------
    dataDict: dict
        The data dictionary to search.
    mapList: list
        The list of keys to follow.
    value: Any
        The value to set the object to
    """
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value


def rsetattr(obj: Any, attr: str, val: Any):
    """
    Recursively set an attribute.

    Parameters
    ----------
    obj: Any
        The object to search.
    attr: str
        The attribute position string.
    val: Any
        The value to set.
    """
    pre, _, post = attr.rpartition(".")
    return setattr(rgetattr(obj, pre) if pre else obj, post, val)


def rgetattr(obj: Any, attr: str, *args) -> Any:
    """
    Recursively get an attribute.

    Parameters
    ----------
    obj: Any
        The object to search.
    attr: str
        The attribute position string.
    """

    def _getattr(obj, attr):
        return getattr(obj, attr, *args)

    return functools.reduce(_getattr, [obj] + attr.split("."))


def produce_dict_heading(dataDict: Mapping, mapList: Iterable[slice]):
    """
    Generate iterative sub-dict in ``dataDict`` to allow setting of a sub-path for which the parent dict doesn't
    yet exist.

    Parameters
    ----------
    dataDict: dict
        The data dictionary to search.
    mapList: list
        The position to generate to.

    Returns
    -------
    None
    """
    for i, _ in enumerate(mapList, start=1):
        try:
            assert getFromDict(dataDict, mapList[:i]) is not None
        except Exception:
            setInDict(dataDict, mapList[:i], {})


def find_descriptors(cls: Any, descriptor_classes: tuple[Any]) -> dict[str, Any]:
    """Searches and returns all descriptors attached to this class or inherited."""
    _descriptor_instances = {
        k: v for k, v in cls.__dict__.items() if isinstance(v, descriptor_classes)
    }
    # iterate through all parent classes.
    _bases = cls.__bases__  # list of all parents.
    while object not in _bases:
        _nbases = []
        for _base in _bases:
            _descriptor_instances = {
                **_descriptor_instances,
                **{
                    k: v
                    for k, v in _base.__dict__.items()
                    if isinstance(v, descriptor_classes)
                },
            }
            _nbases += _base.__bases__
        _bases = _nbases[:]
    return _descriptor_instances


def merge_dicts(template: Mapping, partial: Mapping) -> Mapping:
    _output = {**template}  # Copy template.

    for k, v in template.items():
        if k in partial and v is not None:
            if isinstance(v, dict):
                _output[k] = merge_dicts(v, partial[k])
            else:
                _output[k] = partial[k]

        else:
            pass

    return _output
