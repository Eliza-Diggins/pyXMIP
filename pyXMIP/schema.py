"""
Schema classes for instance specific meta data specification in various pyXMIP classes.

For a more complete discussion of schema, visit :ref:`schema`.
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
    _required_keys = []

    def __init__(self, mapping, *args, **kwargs):
        super().__init__(mapping)

        assert (
            self._check_valid()
        ), "The mapping provided is not valid for the schema type."

    def _check_valid(self):
        if all(k in self.keys() for k in self.__class__._required_keys):
            return True
        else:
            return False

    def isvalid(self):
        return self._check_valid()

    @classmethod
    def from_file(cls, path, format=None, **kwargs):
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


class SkyCollectionSchema(Schema):
    _required_keys = ["column_map", "default_coord_system"]

    _coordinate_requirements = {
        astro_coords.ICRS: ["RA", "DEC"],
        astro_coords.Galactic: ["l", "b"],
    }

    @property
    def colmap(self):
        return self["column_map"]

    @property
    def inv_colmap(self):
        return {v: k for k, v in self["column_map"].items()}

    @property
    def coordinate_system(self):
        _cs = getattr(astro_coords, self["default_coord_system"])
        return _cs, self.available_coordinate_frames()[_cs]

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

        self["default_coord_system"] = value

    @property
    def redshift_column(self):
        if "Z" in self.colmap:
            return self.colmap["Z"]
        else:
            return None

    @redshift_column.setter
    def redshift_column(self, value):
        self.colmap["Z"] = value

    def available_coordinate_frames(self):
        _available_coordinate_frames = self._get_coordinate_frames(
            list(self.colmap.keys())
        )

        if _available_coordinate_frames is not None:
            _available_coordinate_frames = {
                k: [self.colmap[u] for u in v]
                for k, v in _available_coordinate_frames.items()
            }

        return _available_coordinate_frames

    @classmethod
    def construct(cls, table):
        mainlog.debug(f"Constructing {cls.__name__} from fits table.")
        # ----------------------------------------------#
        # setup components
        _table_columns = table.columns
        _constructed_schema = {k: None for k in cls._required_keys}
        _default_coordinate_frame = None

        # ---------------------------#
        # Coordinate Determination
        _construction_coords_avail = cls._get_coordinate_frames(_table_columns)

        assert len(
            _construction_coords_avail
        ), "No valid coordinate system could be determined."
        _default_coordinate_frame = list(_construction_coords_avail.keys())[0].__name__
        _constructed_schema["default_coord_system"] = _default_coordinate_frame
        mainlog.debug(
            f"Located {len(_construction_coords_avail)} possible coordinate frames. Selected {_default_coordinate_frame} as default."
        )

        _constructed_schema["column_map"] = {}
        for _, v in _construction_coords_avail.items():
            for h in v:
                _constructed_schema["column_map"][h] = h

        return cls(_constructed_schema)

    @classmethod
    def _get_coordinate_frames(cls, columns):
        _allowed_coordinate_frames = {}
        for k, v in cls._coordinate_requirements.items():
            if all(u in columns for u in v):
                # this is an allowed coordinate mapping
                _allowed_coordinate_frames[k] = v
        if len(_allowed_coordinate_frames):
            return _allowed_coordinate_frames
        else:
            return None


if __name__ == "__main__":
    from astropy.table import Table

    t = Table.read("/home/ediggins/pyROSITA_test/eRASS1_Hard.v1.0.fits", format="fits")
    t["l"] = t["LII"]
    t["b"] = t["BII"]
    q = SkyCollectionSchema.construct(t)
    print(q)
    print(q.coordinate_system)
    q.coordinate_system = "Galactic"
    print(q.coordinate_system)
    q.coordinate_system = "Ecliptic"
