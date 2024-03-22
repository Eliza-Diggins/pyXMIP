"""
Schema classes for backend management of various user-facing data classes.
"""
import pathlib as pt
from copy import copy

import astropy.coordinates as astro_coords

from pyXMIP.utilities.core import mainlog


def _get_recursive(key, dictionary):
    o = copy(dictionary)

    for k in key.split("."):
        assert (
                k in o
        ), f"The key {k} failed to appear in the schema at this level: KEYS={o.keys()}."
        o = o[k]

    return o


class Schema(dict):
    """
    The archetypal schema class. This is effectively a wrapper on the :py:class:`dict` type with a few additional
    utility methods.

    .. warning::

        This class is archetypal. Most methods are only functional for sub-classes.
    """
    _required_keys = []

    def __init__(self, mapping, *args, **kwargs):
        super().__init__(mapping)

        assert (
            self._check_valid()
        ), "The mapping provided is not valid for the schema type."

    def _check_valid(self):
        from pyXMIP.utilities.core import getFromDict
        for k in self.__class__._required_keys:
            if '.' in k:
                k = k.split('.')
            else:
                k = [k]

            try:
                _ = getFromDict(self,k)
            except KeyError:
                raise ValueError(f" [{self.__class__.__name__}] Failed to find entry {k} on schema check.")

        return True

    def isvalid(self):
        """Returns ``True`` if the schema is a valid format; ``False`` otherwise."""
        return self._check_valid()

    @classmethod
    def from_file(cls, path, format=None, **kwargs):
        """
        Load a :py:class:`schema.Schema` instance from file.

        Parameters
        ----------
        path: str
            The path to the file to load as a schema.
        format: str, optional
            The format of the file being read. If ``None`` (default), then the format is determined based on the file extension.
            Currently, ``yaml``, ``toml``, and ``json`` are all supported formats.
        kwargs
            Additional kwargs to pass to the ``__init__`` call.

        Returns
        -------

        """
        if format is None:
            # determine the format from the extension.
            _suffix = pt.Path(path).suffix

            if _suffix in [".yaml", ".yml", ".YAML", ".YML", ".schema"]:
                # the format is yaml
                format = "yaml"
            elif _suffix in [".tml", ".toml", ".TOML", ".TML", ".config", ".cfg"]:
                format = "toml"
            elif _suffix in ["json", "JSON"]:
                format = "json"
            else:
                raise ValueError(
                    f"The format {_suffix} is not recognized. Schema files must be .yaml, .json, or .toml."
                )
        else:
            pass

        try:
            return getattr(cls, f"_from_{format}")(path, **kwargs)
        except AttributeError:
            raise ValueError(
                f"The format {format} is not recognized. Schema files must be .yaml, .json, or .toml."
            )

    @classmethod
    def _from_yaml(cls, path, **kwargs):
        import yaml

        with open(path, "r") as f:
            return cls(yaml.load(f, yaml.FullLoader))

    @classmethod
    def _from_toml(cls, path, **kwargs):
        import tomllib as toml

        return cls(toml.load(path))

    @classmethod
    def _from_json(cls, path, **kwargs):
        import json

        with open(path, "r") as f:
            return cls(json.load(f))

class _SpecialColumn:
    def __set_name__(self, owner, name):
        self._name = name

    def __get__(self,instance,owner):
        try:
            return instance.colmap[self._name]
        except KeyError:
            return None

    def __set__(self,instance,value):
        instance.colmap[self._name] = value


class SourceTableSchema(Schema):
    r"""
    Schema class for representing table naming conventions for :py:class:`table.SourceTable` objects.

    The :py:class:`SourceTableSchema` informs :py:class:`table.SourceTable` of special roles for particular columns in the dataset.
    Additionally, the schema controls the native coordinate system, and the source type mapping from the native convention of the
    source table to the SIMBAD convention adopted in pyXMIP.

    .. admonition:: Format

        The format of the :py:class:`SourceTableSchema` contains 3 top-level headers:

        - ``column_map``: contains mappings from special ``pyXMIP`` column names to their equivalent in the native table.
          - Special keys are ``Z``, ``NAME``, ``TYPE``, ``RA``, ``DEC``, ``L``, ``B``.
        - ``settings``: contains the default settings for the specified schema.
          - ``default_coord_system``: The default coordinate frame (represented as a string).
        - ``object_map``: contains each of the possible object types from the native dataset and corresponds to the SIMBAD convention used
          in ``pyXMIP``.
    """

    _required_keys = ["column_map", "settings.default_coord_system", "object_map"]

    _coordinate_requirements = {
        astro_coords.ICRS    : ["RA", "DEC"],
        astro_coords.Galactic: ["L", "B"]
    }
    _expected_special_columns = {"Z"   : ['Z', 'z', 'REDSHIFT', 'redshift', 'BEST_Z'],
                                "TYPE": ['TYPE', 'type', 'object_type', 'OBJECT_TYPE', 'OTYPES'],
                                "NAME": ['NAME', 'name', 'OBJECT', 'object', 'ID', 'OBJECT_ID','IAUNAME'],
                                "RA": ['RA','ra'],
                                 'DEC': ['DEC','dec'],
                                 'L': ['l','lii','L','LII'],
                                 'B': ['B','b','bii','BII']
                                }

    # -- Constructing the schema column map accessors -- #
    Z = _SpecialColumn()
    TYPE = _SpecialColumn()
    NAME = _SpecialColumn()
    RA = _SpecialColumn()
    DEC = _SpecialColumn()
    B = _SpecialColumn()
    L = _SpecialColumn()


    @property
    def colmap(self):
        """
        The column map containing mappings between the "special" columns (using the native naming convention of pyXMIP) and the
        true column names in the native dataset.

        .. note::

            The column map is a dictionary with ``keys`` corresponding to the pyXMIP designation for the column and the ``values``
            corresponding to the native designation. **Not all columns must be mapped in** ``colmap``. ``colmap`` contains only those column
            mappings with are "special": coordinate columns, object type columns, redshift columns, and object designation columns.

        Returns
        -------
        dict
            The column map.
        """
        return self["column_map"]

    @property
    def inv_colmap(self):
        """
        The inverse column map, now with keys selected from the native columns and keys for the special columns of pyXMIP.

        Returns
        -------
        dict
            The inverse column map.
        """
        return {v: k for k, v in self["column_map"].items()}

    @property
    def coordinate_system(self):
        """
        The specific coordinate system of the :py:class:`schema.SourceTableSchema` instance.

        If :py:meth:`schema.SourceTableSchema.coordinate_system` is set to another value, it may be either ``str`` or
        :py:class:`astropy.coordinates.GenericFrame`, and must be consistent with the available columns as specified by the ``colmap``.

        Returns
        -------
        :py:class:`astropy.coordinates.GenericFrame`
            The correct coordinate frame object for the coordinate system used by the schema.
        list
            The table column names corresponding to the coordinate system.
        """
        _cs = getattr(astro_coords, self['settings']["default_coord_system"])
        return _cs, self.available_coordinate_frames()[_cs]

    def set_special_column(self,key,value):
        assert key in self.__class__._expected_special_columns

    @property
    def coordinate_columns(self):
        """The columns representing the current coordinate system."""
        return self.coordinate_system[1]

    @property
    def coordinate_frame(self):
        """The specific coordinate system currently corresponding to this schema."""
        return self.coordinate_system[0]

    @coordinate_system.setter
    def coordinate_system(self, value):
        # check the typing.
        if isinstance(value, str):
            pass
        elif isinstance(value, astro_coords.GenericFrame):
            value = value.__name__
        else:
            raise TypeError(f"Cannot set coordinate system to type {type(value)}.")

        # check validity
        assert hasattr(
            astro_coords, value
        ), f"The coordinate frame {value} is not in astropy's coordinate systems."
        assert (
                getattr(astro_coords, value) in self.available_coordinate_frames()
        ), f"The coordinate frame {value} is not a valid coordinate frame for this schema."

        self['settings']["default_coord_system"] = value

    def available_coordinate_frames(self):
        """
        List the available :py:class:`astropy.coordinates.GenericFrame` types which can be associated with this schema given
        the available data.

        Returns
        -------
        list
            The available coordinate frames.
        """
        _available_coordinate_frames = self._get_coordinate_frames(
            self.colmap
        )

        if _available_coordinate_frames is not None:
            _available_coordinate_frames = {
                k: [self.colmap[u] for u in v]
                for k, v in _available_coordinate_frames.items()
            }
        return _available_coordinate_frames

    @classmethod
    def construct(cls, table):
        """
        Construct a "best-guess" schema from a table by deducing the meaning of it's columns.

        Parameters
        ----------
        table: :py:class:`structures.table.SourceTable`
            The table on which to determine the best guess schema.

        Returns
        -------
        :py:class:`SourceTableSchema`
            The finalized guess schema.
        """
        mainlog.debug(f" [SourceTableSchema] Constructing {cls.__name__} from fits table.")
        # ----------------------------------------------#
        # setup components
        _table_columns = table.columns
        _constructed_schema = cls._build_schema_template()

        _default_coordinate_frame = None

        # ----------------------------------------------------------------------#
        # Managing special columns
        # ----------------------------------------------------------------------#
        _constructed_schema['column_map'] = {}
        for k,v in cls._expected_special_columns.items():
            _matched_keys = [u for u in v if u in _table_columns]

            if len(_matched_keys):
                mainlog.debug(f" [SourceTableSchema] Identified special key {k} with column {_matched_keys[0]} of the table.")
                _constructed_schema['column_map'][k] = _matched_keys[0]
            else:
                mainlog.debug(f" [SourceTableSchema] Failed to identify automatic match for special column {k}.")

        # ----------------------------------------------------------------------#
        # Managing object maps
        # ----------------------------------------------------------------------#
        _constructed_schema['object_map'] = {}

        # ----------------------------------------------------------------------#
        # Managing coordinates
        # ----------------------------------------------------------------------#
        _construction_coords_avail = cls._get_coordinate_frames(_constructed_schema['column_map'])
        assert len(
            _construction_coords_avail
        ), " [SourceTableSchema] No valid coordinate system could be determined."

        _default_coordinate_frame = list(_construction_coords_avail.keys())[0].__name__
        _constructed_schema['settings']["default_coord_system"] = _default_coordinate_frame

        mainlog.debug(
            f" [SourceTableSchema] Located {len(_construction_coords_avail)} possible coordinate frames. Selected {_default_coordinate_frame} as default."
        )

        for _, v in _construction_coords_avail.items():
            for h in v:
                _constructed_schema["column_map"][h] = h

        return cls(_constructed_schema)

    @classmethod
    def _get_coordinate_frames(cls, column_map):
        """searches the column map for allowed values"""
        _allowed_coordinate_frames = {}
        for k, v in cls._coordinate_requirements.items():
            if all(u in column_map for u in v):
                # this is an allowed coordinate mapping
                _allowed_coordinate_frames[k] = [column_map[_v] for _v in v]
        if len(_allowed_coordinate_frames):
            return _allowed_coordinate_frames
        else:
            return None

    @classmethod
    def _build_schema_template(cls):
        """construct a schema from the template"""
        from pyXMIP.utilities.core import getFromDict, setInDict
        _output_schema = {}

        for key in cls._required_keys:
            # parsing the keys and subkeys.
            if "." in key:
                subkeys = key.split(".")

                for k_id, key in enumerate(subkeys[:-1], 1):
                    try:
                        _ = getFromDict(_output_schema, subkeys[:k_id])
                    except KeyError:
                        setInDict(_output_schema, subkeys[:k_id], {})

                setInDict(_output_schema, subkeys, None)
            else:
                _output_schema[key] = None

        return _output_schema

