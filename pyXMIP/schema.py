"""
Schema management module for ``pyXMIP``

The :py:mod:`schema` module provides structure and data validation processes for user specified objects. In most cases,
this functions to effectively allow ``pyXMIP`` to convert from the user's data convention to those used internally.

Notes
-----

For a long-form discussion about schema, see :ref:`schema`.

"""
import os
import pathlib as pt
from typing import ClassVar, Type, Union

import astropy.coordinates as astro_coords
import pydantic as pyn

from pyXMIP.utilities.core import bin_directory, find_descriptors, rgetattr, rsetattr
from pyXMIP.utilities.logging import mainlog
from pyXMIP.utilities.types import (
    ColumnMap,
    ICRSCoordinateStdErrorSpecifier,
    Registry,
    SourceTableSchemaSettings,
    TableColumn,
    construct_template,
)

#: directory containing the built-in schema.
builtin_schema_directory = os.path.join(bin_directory, "builtin_schema")
#: directory containing the built-in :py:class:`SourceTableSchema` instances.
builtin_source_table_schema_directory = os.path.join(
    builtin_schema_directory, "source_table"
)
#: directory containing the built-in :py:class:`ReductionSchema` instances.
builtin_reduction_schema_directory = os.path.join(builtin_schema_directory, "reduction")


class SchemaEntry:
    """
    Descriptor class for accessing variables in the :py:class:`Schema` class.
    """

    def __init__(self, dict_location: str = None):
        """
        Initialize the entry.

        Parameters
        ----------
        dict_location: str
            The position of the underlying attribute in the :py:class:`Schema` instance. For example,
            the redshift column is ``'schema.column_map.Z'``, so ``dict_location = 'column_map'``.
        """
        self.dict_location = dict_location

    def __set_name__(self, owner, name):
        self._name = name

        if self.dict_location is None:
            self.dict_location = self._name
        else:
            self.dict_location += f".{self._name}"

    def __get__(self, instance, owner):
        try:
            return rgetattr(instance, self.dict_location)
        except AttributeError:
            return None

    def __set__(self, instance, value):
        rsetattr(instance, self.dict_location, value)


class ColumnEntry:
    """
    Specialized descriptor for accessing :py:class:`Schema` columns.
    """

    def __init__(self, auto_id=None, default_unit=None):
        self.auto_id = auto_id if auto_id is not None else []
        self.default_unit = default_unit

    def __set_name__(self, owner, name):
        self._name = name

    def __get__(self, instance, owner):
        try:
            return getattr(instance.column_map, self._name).name
        except AttributeError:
            return None

    def __set__(self, instance, value):
        if isinstance(value, TableColumn):
            setattr(instance.column_map, self._name, value)
        elif isinstance(value, dict):
            setattr(instance.column_map, self._name, TableColumn(**value))
        else:
            setattr(instance.column_map, self._name, TableColumn(name=value))


class Schema(pyn.BaseModel):
    """
    The generic schema class from which all other ``pyXMIP`` schemas are derived.

    .. warning::

        This is an abstract class.

    Notes
    -----

    The :py:class:`Schema` class is simply a wrapper on the PyDantic :py:class:`pydantic.BaseModel` class. This allows
    for the validation of user provided data.

    """

    def __init__(self, **mapping):
        """
        Initialize the generalized schema class.

        Parameters
        ----------
        mapping:
            Keyword arguments to pass into the schema.
        """
        super().__init__(**mapping)

    def __repr__(self):
        return f"<{self.__class__.__name__} instance>"

    def write(
        self,
        filename: Union[str, pt.Path],
        overwrite: bool = False,
        file_format: str = "yaml",
    ):
        """
        Write this schema instance to file. Various data formats are permitted.

        Parameters
        ----------
        filename:str
            The path at which to write the schema.
        overwrite: bool, optional
            If ``True``, then write will occur regardless of the prior existence of the file. Otherwise, an error is
            raised.
        file_format: str, optional
            The format to select for outputting the model.

        Raises
        ------
        IOError
            Raised if the specified ``filename`` cannot be reached.
        AssertionError
            Raised if the file already exists and ``overwrite=False``.

        """
        # ----------------------------------------------------- #
        # Arg / kwarg validation                                #
        # ----------------------------------------------------- #
        path = pt.Path(filename)

        # check the path for validity. The upstream directory must already exist.
        if not path.parents[0].exists():
            raise IOError(f"Parent directory(ies) of {path} don't exist.")

        # Check if the file already exists and has data.
        if not overwrite:
            # We don't have overwrite. Check if we can proceed.
            assert (
                not path.exists()
            ), f"Cannot write to {filename} because it already exists and overwrite is False."
        elif overwrite and path.exists():
            # We have overwrite = True and the file exists, so we delete it.
            path.unlink()

        # ----------------------------------------------------- #
        # Formatting Determination                              #
        # ----------------------------------------------------- #
        match file_format:  # Check the various formats.
            case "yaml":
                # This will be output to yaml.
                import yaml

                _model_dict = self.dict()

                with open(filename, "w") as f:
                    yaml.dump(_model_dict, f)

            case "json":
                import json

                _model_str = self.model_dump(mode="json")

                with open(filename, "w") as f:
                    json.dump(_model_str, f)

            case "toml":
                import toml

                _model_dict = self.dict()

                with open(filename, "w") as f:
                    toml.dump(_model_dict, f)
            case "txt":
                _model_str = str(self.dict())

                with open(filename, "w") as f:
                    f.write(_model_str)

    @classmethod
    def read(
        cls, filename: Union[str, pt.Path], file_format: str = None, *args, **kwargs
    ):
        """
        Read schema from file. Various data formats are permitted.

        Parameters
        ----------
        filename:str
            The path at which to read the schema.
        file_format: str, optional
            The intended format for the schema being read. If ``None``, then ``pyXMIP`` will attempt
            to determine it based on the suffix.

        Raises
        ------
        IOError
            Raised if the specified ``filename`` cannot be reached.
        AssertionError
            Raised if the file already exists and ``overwrite=False``.

        """
        # ----------------------------------------------------- #
        # Arg / kwarg validation                                #
        # ----------------------------------------------------- #
        path = pt.Path(filename)

        # check the path for validity. The upstream directory must already exist.
        if not path.exists():
            raise IOError(f"{path} doesn't exist.")

        # ----------------------------------------------------- #
        # Formatting Determination                              #
        # ----------------------------------------------------- #
        if file_format is None:
            file_format = path.suffix[1:]

        match file_format:  # Check the various formats.
            case "yaml":
                # This will be output to yaml.
                import yaml

                with open(filename, "r") as f:
                    if not len(args):
                        # the args should have specified the loader. Because they're not provided, we assume full loader.
                        args = [yaml.FullLoader]

                    _r = yaml.load(f, *args, **kwargs)

            case "json":
                import json

                with open(filename, "r") as f:
                    _r = json.load(f, *args, **kwargs)
            case "toml":
                import toml

                with open(filename, "r") as f:
                    _r = toml.load(f, *args, **kwargs)

            case "txt":
                with open(filename, "r") as f:
                    _r = {k: v for k, v in [u.split(":") for u in f.read().split(",")]}
            case _:
                raise ValueError(f"Failed to recognize format {file_format}.")

        return cls(**_r)


class SourceTableSchema(Schema):
    """
    Schema for interactions with :py:class:`structures.table.SourceTable`.
    """

    # =========================================== #
    # Validated parameters                        #
    # =========================================== #
    column_map: ColumnMap
    """:py:class:`utilities.types.ColumnMap`: Mapping between table columns and standard recognized by ``pyXMIP``.

    The column map instructs ``pyXMIP`` about what different columns in a table mean.
    """
    object_map: dict | None = None
    """dict: The mapping between the native source object types and the standardized types.

    This parameter is used to determine how to map the typing from the catalog / table into the standardized
    SIMBAD types used by ``pyXMIP``.

    .. hint::

        You can find a complete list of SIMBAD object types `here <https://simbad.cds.unistra.fr/guide/otypes.htx>`_.
    """
    settings: SourceTableSchemaSettings = SourceTableSchemaSettings()
    """:py:class:`utilities.types.SourceTableSchemaSettings`: The schema settings for this instance."""

    # =========================================== #
    # Columns                                     #
    # =========================================== #
    # Coordinate Columns #
    RA: ClassVar[str] = ColumnEntry(auto_id=["RA", "ra", "Ra"], default_unit="deg")
    """The column containing the object RA coordinate."""
    DEC: ClassVar[str] = ColumnEntry(auto_id=["dec", "DEC", "Dec"], default_unit="deg")
    """The column containing the object DEC coordinate."""
    RA_ERR: ClassVar[str] = ColumnEntry(auto_id=[], default_unit="arcsec")
    """The column containing the object RA coordinate error."""
    DEC_ERR: ClassVar[str] = ColumnEntry(auto_id=[], default_unit="arcsec")
    """The column containing the object DEC coordinate error."""
    GAL_L: ClassVar[str] = ColumnEntry(
        auto_id=["Gal_L", "GL", "LII"], default_unit="deg"
    )
    """The column containing the object's galactic longitude."""
    GAL_B: ClassVar[str] = ColumnEntry(
        auto_id=["GAL_B", "GB", "BII"], default_unit="deg"
    )
    """The column containing the object's galactic latitude."""
    GAL_L_ERR: ClassVar[str] = ColumnEntry(auto_id=[], default_unit="arcsec")
    """The column containing the object's galactic longitude error."""
    GAL_B_ERR: ClassVar[str] = ColumnEntry(auto_id=[], default_unit="arcsec")
    """The column containing the object's galactic latitude error."""
    SEP: ClassVar[str] = ColumnEntry(auto_id=["Separation"], default_unit="arcsec")
    # Name columns #
    NAME: ClassVar[str] = ColumnEntry(
        auto_id=["IAUNAME", "NAME", "MAIN_ID", "OBJECT_ID"]
    )
    """The column containing the object's main identifier."""
    # Extra Columns #
    Z: ClassVar[str] = ColumnEntry(auto_id=["REDSHIFT", "redshift", "Z", "z"])
    """The column containing the object's redshift."""
    TYPE: ClassVar[str] = ColumnEntry(auto_id=["TYPE"])
    """The column containing the object type of the catalog object."""

    # ============================================ #
    # Settings                                     #
    # ============================================ #
    default_coord_system: ClassVar[str] = SchemaEntry(dict_location="settings")
    """str: The default coordinate system to use.

    Use this to select the native coordinate system for the schema. Whichever coordinate system is selected, there must
    be sufficient identifiable columns to actually obtain those coordinates.

    """
    object_type_separator: ClassVar[str] = SchemaEntry(dict_location="settings")
    """str: The separator character / string between different object types.

    In many databases, returned object types may have multiple types separated by a delimiter. This delimiter should
    be specified here if it is needed.

    .. important::

        In all native ``pyXMIP`` data objects, a "|" is used for the delimiter. Even if the object type is a single
        value, all object types in ``pyXMIP`` tables are formatted ``"|entry_1|entry_2|...|"``.

    """

    _coordinate_requirements: ClassVar[dict] = {
        astro_coords.ICRS: [RA, DEC],
        astro_coords.Galactic: [GAL_L, GAL_B],
    }

    def available_coordinate_frames(self):
        """
        List the available astropy coordinate frames which can be used to determine positions of the data based on
        columns provided to this schema.

        Returns
        -------
        dict
            The available coordinate frames.
        """
        # -- determine which coordinate frames are available -- #
        available_frames = [
            frame
            for frame, cols in self._coordinate_requirements.items()
            if all(getattr(self, col._name) is not None for col in cols)
        ]

        return {
            frame: [u._name for u in self._coordinate_requirements[frame]]
            for frame in available_frames
        }

    @property
    def coordinate_system(self) -> tuple[astro_coords.GenericFrame, list[str]]:
        """
        The coordinate system and associated columns for this schema. Unless this attribute is set after initialization, the
        coordinate system will be the same as determined by :py:attr:`SourceTableSchema.default_coord_system`.

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
    def coordinate_frame(self) -> astro_coords.GenericFrame:
        """The specific coordinate system currently corresponding to this schema."""
        return self.coordinate_system[0]

    @property
    def coordinate_columns(self) -> list[TableColumn]:
        """The columns representing the current coordinate system."""
        # The setup returns the actual descriptor object, we need the underlying parameter (hence the call to .__get__)
        return [
            getattr(self.column_map, coord_col)
            for coord_col in self.coordinate_system[1]
        ]

    @coordinate_system.setter
    def coordinate_system(self, value: Union[str, astro_coords.GenericFrame]):
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

        self.settings.default_coord_system = value

    @classmethod
    def from_table(cls, table):
        r"""
        Construct a "best-guess" schema from a table by deducing the meaning of its columns.

        Parameters
        ----------
        table: :py:class:`structures.table.SourceTable` or :py:class:`pd.DataFrame`
            The table on which to determine the best guess schema.

        Returns
        -------
        :py:class:`SourceTableSchema`
            The finalized guess schema.
        """
        mainlog.debug(f" [SourceTableSchema] Constructing {cls.__name__} from table.")

        # ----------------------------------------------#
        # setup components
        _table_columns = table.columns  # isomorphic call between table types.
        _constructed_schema = construct_template(cls)

        _default_coordinate_frame = None

        # ----------------------------------------------------------------------#
        # Managing special columns
        # ----------------------------------------------------------------------#
        _constructed_schema["column_map"] = {}
        for k, v in find_descriptors(cls, ColumnEntry).items():
            # search for matching keys
            _matched_keys = [key for key in _table_columns if key in v.auto_id]

            if len(_matched_keys):
                # A matching key was found.
                _constructed_schema["column_map"][k] = {"name": _matched_keys[0]}

                # try to fetch a unit.
                if hasattr(table[_matched_keys[0]], "unit"):
                    _constructed_schema["column_map"][k]["unit"] = table[
                        _matched_keys[0]
                    ].unit
                else:
                    _constructed_schema["column_map"][k]["unit"] = v.default_unit

                if _constructed_schema["column_map"][k]["unit"] is None:
                    del _constructed_schema["column_map"][k]["unit"]

                mainlog.debug(
                    f" [SourceTableSchema] Identified special column {k} with column {_matched_keys[0]} [unit={_constructed_schema['column_map'][k].get('unit','')}] of the table."
                )

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
        _construction_coords_avail = [
            frame
            for frame, cols in cls._coordinate_requirements.items()
            if all(col._name in _constructed_schema["column_map"] for col in cols)
        ]
        assert len(
            _construction_coords_avail
        ), " [SourceTableSchema] No valid coordinate system could be determined."

        _default_coordinate_frame = _construction_coords_avail[0].__name__
        _constructed_schema["settings"].default_coord_system = _default_coordinate_frame

        mainlog.debug(
            f" [SourceTableSchema] Located {len(_construction_coords_avail)} possible coordinate frames. Selected {_default_coordinate_frame} as default."
        )
        return cls(**_constructed_schema)

    def guess_coordinate_standard_error_struct(self):
        # -- build the base template -- #
        template = construct_template(ICRSCoordinateStdErrorSpecifier)

        # --------------------------- #
        # Look for valid coordinates  #
        # --------------------------- #
        _available_coordinate_frames = self.available_coordinate_frames()

        # For each of the available frames, we need to determine if the coordinate errors are available.
        _frames_with_error_cols = {}
        for frame, required_columns in _available_coordinate_frames.items():
            required_error_columns = [
                getattr(self.column_map, u + "_ERR")
                if hasattr(self.column_map, u + "_ERR")
                else None
                for u in required_columns
            ]
            required_columns = [
                getattr(self.column_map, u) for u in required_columns
            ]  # Warning! note the change in type.

            if any(u is None for u in required_error_columns):
                continue
            else:
                _frames_with_error_cols[frame] = (
                    required_error_columns,
                    required_columns,
                )

        if not len(_frames_with_error_cols):
            # We can't find a way to improve on the base template. Return it blank.
            return template
        else:
            if self.coordinate_system[0] in _frames_with_error_cols:
                err_cols = (
                    self.coordinate_system[0],
                    _frames_with_error_cols[self.coordinate_system[0]],
                )
            else:
                err_cols = [(k, v) for k, v in _frames_with_error_cols.items()][0]

        # ------------------------------------ #
        # Adding error columns to the template #
        # ------------------------------------ #
        template["longitude_column"], template["latitude_column"] = (
            err_cols[1][1][0].dict(),
            err_cols[1][1][1].dict(),
        )
        (
            template["longitude_column_error"]["base_column"],
            template["latitude_column_error"]["base_column"],
        ) = (err_cols[1][1][0].dict(), err_cols[1][1][1].dict())
        (
            template["longitude_column_error"]["error_column"],
            template["latitude_column_error"]["error_column"],
        ) = (err_cols[1][0][0].dict(), err_cols[1][0][1].dict())

        return template


class CMDSchema(Schema):
    """
    Schema class for interaction with :py:class:`cross_reference.CrossMatchDatabase`.

    Notes
    -----

    This class is effectively a collection of :py:class:`SourceTableSchema` instances compiled together to
    represent the entire :py:class:`cross_reference.CrossMatchDatabase` structure.

    """

    table_schema: dict[str, SourceTableSchema]
    """dict: Dictionary of schema specific to each of the cross matching databases.

    For each of the tables in the underlying ``SQL`` database (excluding ``CATALOG``), there should be a schema corresponding
    to the table format for the source. For the most part, these are the :py:attr:`structures.databases.RemoteDatabase.default_query_schema` of
    the queried databases.

    .. hint::

        For the most part, these are generated automatically when the cross matching occurs and you don't need to worry about
        generating one of these on your own. If you do, the dictionary has a simple format: keys are the NAMES of the tables in
        the :py:class:`cross_reference.CrossMatchDatabase`, the values should be the corresponding schema.
    """
    catalog_schema: SourceTableSchema
    """:py:class:`SourceTableSchema` The schema for the catalog.

    This is the schema corresponding to the original catalog. In almost all cases, this is done automatically.
    """


class ReductionSchema(Schema):
    pass


class SchemaRegistry(Registry):
    """
    Registry class for containing collections of :py:class:`Schema` classes.
    """

    def __init__(self, mapping: dict[str, str | pt.Path], schema_class: Type[Schema]):
        """
        Initialize the :py:class:`SchemaRegistry`

        Parameters
        ----------
        mapping: dict
            The dictionary of ``name``,``path`` pairs to use as the registry.
        schema_class:
            The schema class to initialize objects with.
        """
        super().__init__(mapping)

        self.schema_class: Type[Schema] = schema_class

    def __getitem__(self, item: str) -> Schema:
        return self.schema_class.read(super().__getitem__(item))

    @classmethod
    def from_directory(cls, directory, schema_class: Type[Schema]):
        """
        Generate a registry of :py:class:`Schema` instances from a directory of files.

        Parameters
        ----------
        directory: str
            The directory to load from.
        schema_class:
            The schema class for this registry.

        Returns
        -------
        :py:class:`SchemaRegistry`
            The corresponding registry.
        """
        import pathlib as pt

        directory = pt.Path(directory)
        assert directory.exists(), f"The directory {directory} doesn't exist."
        assert directory.is_dir(), f"The path {directory} isn't a directory."

        sub_objects = os.listdir(directory)
        _dict = {}

        for o in sub_objects:
            pth = pt.Path(os.path.join(directory, o))
            if pth.is_file():
                suffix = pth.suffix

                if suffix in [".toml", ".json", ".yaml"]:
                    _dict[pth.name.replace(pth.suffix, "")] = str(pth)

        return cls(_dict, schema_class)


#: The default :py:class:`SourceTableSchema` registry.
DEFAULT_SOURCE_SCHEMA_REGISTRY = SchemaRegistry.from_directory(
    builtin_source_table_schema_directory, SourceTableSchema
)

#: The default :py:class:`ReductionSchema` registry.
DEFAULT_REDUCTION_SCHEMA_REGISTRY = SchemaRegistry.from_directory(
    builtin_reduction_schema_directory, ReductionSchema
)


if __name__ == "__main__":
    NED_schema = DEFAULT_SOURCE_SCHEMA_REGISTRY["NED"]
    print(type(NED_schema.model_dump(mode="json")))
