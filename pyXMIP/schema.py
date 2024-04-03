"""
Schema classes for backend management of various user-facing data classes.
"""
import pathlib as pt
from copy import copy

import astropy.coordinates as astro_coords
from astropy.units import Quantity

from pyXMIP.utilities.core import getFromDict, mainlog, setInDict


class _SCHEMAENTRY:
    def __init__(self, dict_location, required=False, dtype=None, allow_set=True):
        self.dict_location = dict_location
        self.required = required
        self.dtype = dtype
        self.allow_set = allow_set

    def __set_name__(self, owner, name):
        self._name = name
        self.dict_location += [self._name]

    def __get__(self, instance, owner):
        try:
            return getFromDict(instance, self.dict_location)
        except KeyError:
            return None
        except TypeError:
            return None

    def __set__(self, instance, value):
        setInDict(instance, self.dict_location, value)


class _SpecialColumn:
    def __set_name__(self, owner, name):
        self._name = name

    def __get__(self, instance, owner):
        try:
            return instance.column_map[self._name]
        except KeyError:
            return None

    def __set__(self, instance, value):
        instance.column_map[self._name] = value


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

    def __init__(self, mapping):
        """
        Initialize the generalized schema class.

        Parameters
        ----------
        mapping: dict
            The dictionary underlying this schema.
        """
        super().__init__(mapping)

        assert (
            self._check_valid()
        ), "The mapping provided is not valid for the schema type."

    @classmethod
    def _get_descriptors(cls):
        # fetching descriptors.
        return [m for m, v in cls.__dict__.items() if isinstance(v, _SCHEMAENTRY)]

    @classmethod
    def _get_required_descriptors(cls):
        return [
            m
            for m, v in cls.__dict__.items()
            if isinstance(v, _SCHEMAENTRY) and v.required is True
        ]

    def _check_valid(self):
        for schema_entry in self.__class__._get_descriptors():
            # Iterate through the schema entries associated with this object.
            _ref = self.__class__.__dict__[schema_entry]

            if isinstance(_ref.required, bool) and _ref.required:
                assert (
                    getattr(self, schema_entry) is not None
                ), f"Failed to find required entry {schema_entry} at {_ref.dict_location} in the provided schema."
            elif isinstance(_ref.required, str):
                _parent_val = getattr(self, _ref.required)

                if _parent_val is not None:
                    assert (
                        getattr(self, schema_entry) is not None
                    ), f"Failed to find required (BY PARENT {_ref.required}) entry {schema_entry} at {_ref.dict_location} in the provided schema."
            _val = getattr(self, schema_entry)

            if _val is not None:
                if _ref.dtype is not None:
                    assert isinstance(
                        _val, _ref.dtype
                    ), f"Found {schema_entry}, but it had type {type(_val)}, not {_ref.dtype}."

        return True

    def isvalid(self):
        """Returns ``True`` if the schema is a valid format; ``False`` otherwise."""
        return self._check_valid()

    @classmethod
    def from_file(cls, path, format=None):
        """
        Load a :py:class:`schema.Schema` instance from file.

        Parameters
        ----------
        path: str
            The path to the file to load as a schema.
        format: str, optional
            The format of the file being read. If ``None`` (default), then the format is determined based on the file extension.
            Currently, ``yaml``, ``toml``, and ``json`` are all supported formats.
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
            return getattr(cls, f"_from_{format}")(path)
        except AttributeError:
            raise ValueError(
                f"The format {format} is not recognized. Schema files must be .yaml, .json, or .toml."
            )

    @classmethod
    def _from_yaml(cls, path):
        import yaml

        with open(path, "r") as f:
            return cls(yaml.load(f, yaml.FullLoader))

    @classmethod
    def _from_toml(cls, path):
        import tomllib as toml

        return cls(toml.load(path))

    @classmethod
    def _from_json(cls, path):
        import json

        with open(path, "r") as f:
            return cls(json.load(f))

    @classmethod
    def _build_schema_template(cls):
        """construct a schema from the template"""

        _output_schema = {}

        for schema_entry in cls._get_required_descriptors():
            _loc = cls.__dict__[schema_entry].dict_location
            for _l in range(1, len(_loc)):
                setInDict(_output_schema, _loc[:_l], {})

            setInDict(_output_schema, cls.__dict__[schema_entry].dict_location, None)

        return _output_schema


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
    _coordinate_requirements = {
        astro_coords.ICRS: ["RA", "DEC"],
        astro_coords.Galactic: ["L", "B"],
    }
    _expected_special_columns = {
        "Z": ["Z", "z", "REDSHIFT", "redshift", "BEST_Z"],
        "TYPE": ["TYPE", "type", "object_type", "OBJECT_TYPE", "OTYPES"],
        "NAME": ["NAME", "name", "OBJECT", "object", "ID", "OBJECT_ID", "IAUNAME"],
        "RA": ["RA", "ra"],
        "DEC": ["DEC", "dec"],
        "L": ["l", "lii", "L", "LII"],
        "B": ["B", "b", "bii", "BII"],
    }

    # -- Construct the attributes -- #
    column_map = _SCHEMAENTRY([], required=True, dtype=dict, allow_set=True)
    default_coord_system = _SCHEMAENTRY(
        ["settings"], required=True, dtype=str, allow_set=True
    )
    object_map = _SCHEMAENTRY([], required=True, dtype=dict, allow_set=True)

    # -- Constructing the schema column map accessors -- #
    Z = _SpecialColumn()
    TYPE = _SpecialColumn()
    NAME = _SpecialColumn()
    RA = _SpecialColumn()
    DEC = _SpecialColumn()
    B = _SpecialColumn()
    L = _SpecialColumn()

    @property
    def inv_column_map(self):
        """
        The inverse column map, now with keys selected from the native columns and keys for the special columns of pyXMIP.

        Returns
        -------
        dict
            The inverse column map.
        """
        return {v: k for k, v in self.column_map.items()}

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
        _cs = getattr(astro_coords, self.default_coord_system)
        return _cs, self.available_coordinate_frames()[_cs]

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

        self.default_coord_system = value

    def available_coordinate_frames(self):
        """
        List the available :py:class:`astropy.coordinates.GenericFrame` types which can be associated with this schema given
        the available data.

        Returns
        -------
        list
            The available coordinate frames.
        """
        _available_coordinate_frames = self._get_coordinate_frames(self.column_map)

        if _available_coordinate_frames is not None:
            _available_coordinate_frames = {
                k: [self.column_map[u] for u in v]
                for k, v in _available_coordinate_frames.items()
            }
        return _available_coordinate_frames

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
        mainlog.debug(
            f" [SourceTableSchema] Constructing {cls.__name__} from fits table."
        )
        # ----------------------------------------------#
        # setup components
        _table_columns = table.columns
        _constructed_schema = cls._build_schema_template()

        _default_coordinate_frame = None

        # ----------------------------------------------------------------------#
        # Managing special columns
        # ----------------------------------------------------------------------#
        _constructed_schema["column_map"] = {}
        for k, v in cls._expected_special_columns.items():
            _matched_keys = [u for u in v if u in _table_columns]

            if len(_matched_keys):
                mainlog.debug(
                    f" [SourceTableSchema] Identified special key {k} with column {_matched_keys[0]} of the table."
                )
                _constructed_schema["column_map"][k] = _matched_keys[0]
            else:
                mainlog.debug(
                    f" [SourceTableSchema] Failed to identify automatic match for special column {k}."
                )

        # ----------------------------------------------------------------------#
        # Managing object maps
        # ----------------------------------------------------------------------#
        _constructed_schema["object_map"] = {}

        # ----------------------------------------------------------------------#
        # Managing coordinates
        # ----------------------------------------------------------------------#
        _construction_coords_avail = cls._get_coordinate_frames(
            _constructed_schema["column_map"]
        )
        assert len(
            _construction_coords_avail
        ), " [SourceTableSchema] No valid coordinate system could be determined."

        _default_coordinate_frame = list(_construction_coords_avail.keys())[0].__name__
        _constructed_schema["settings"][
            "default_coord_system"
        ] = _default_coordinate_frame

        mainlog.debug(
            f" [SourceTableSchema] Located {len(_construction_coords_avail)} possible coordinate frames. Selected {_default_coordinate_frame} as default."
        )

        for _, v in _construction_coords_avail.items():
            for h in v:
                _constructed_schema["column_map"][h] = h

        return cls(_constructed_schema)


class ReductionSchema(Schema):
    # -- RUN_PARAMS -- #
    RUN_PARAMS = _SCHEMAENTRY([], required=True, dtype=dict)
    REFDBS = _SCHEMAENTRY(["RUN_PARAMS"], required=True, dtype=(list, str))
    POISSON = _SCHEMAENTRY(["RUN_PARAMS"], required=False, dtype=bool)
    INSTRUMENTAL = _SCHEMAENTRY(["RUN_PARAMS"], required=False, dtype=bool)
    OTYPES = _SCHEMAENTRY(["RUN_PARAMS"], required=False, dtype=bool)

    # -- IO PARAMS -- #
    IO_PARAMS = _SCHEMAENTRY([], required=True, dtype=dict)
    DBPATH = _SCHEMAENTRY(["IO_PARAMS"], required=True, dtype=str)

    # -- POISSON PARAMS -- #
    POISSON_PARAMS = _SCHEMAENTRY([], required="POISSON", dtype=dict)

    # -- INSTRUMENT PARAMS -- #
    INSTRUMENT_PARAMS = _SCHEMAENTRY([], required="INSTRUMENTAL", dtype=dict)
    PSF = _SCHEMAENTRY(
        ["INSTRUMENT_PARAMS"], required="INSTRUMENT_PARAMS", dtype=Quantity
    )

    # -- OBJECT PARAMS -- #
    OBJECT_PARAMS = _SCHEMAENTRY([], required="OTYPES", dtype=dict)


if __name__ == "__main__":
    u = ReductionSchema(
        {
            "RUN_PARAMS": {"INSTRUMENTAL": True},
            "IO_PARAMS": {"DBPATH": ""},
            "INSTRUMENT_PARAMS": {},
        }
    )
