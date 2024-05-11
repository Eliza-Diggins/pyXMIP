"""
``pyXMIP`` schema module for managing conversion between data conventions and ``pyXMIP`` standards.

What is a Schema
----------------

There are a variety of different types of schema in ``pyXMIP``, the most common of which is :py:class:`schema.SourceTableSchema`.
Schema have a variety of purposes, but they all do effectively the same thing: they allow ``pyXMIP`` to convert your data
to a form it understands.

Under-the-hood, schema are *fancy dictionaries* with various properties and methods specific to their particular use-case.
Let's begin by listing some of the most important schema:

- :py:class:`schema.SourceTableSchema` - By far the most important of the schemas.

  - Everytime the user loads a table in ``pyXMIP``, it get an associated :py:class:`schema.SourceTableSchema` instance.

    - Sometimes, this can be deduced from the available column names (i.e. RA and DEC columns); however, in other cases
      it cannot.
    - ``pyXMIP`` will attempt to generate a schema with access to all of the correct parameters; however, if this fails it may
      be necessary to write your own.

  - The corresponding :py:class:`schema.SourceTableSchema` allows ``pyXMIP`` to convert column names to standardized forms, understand
    different object types, and manage a variety of settings.

- :py:class:`schema.ReductionSchema`

  - The reduction process converts raw match data to a statistically optimized cross-reference catalogs. To do this, the user generally
    writes a **reduction file** (``yaml``) containing all of the settings for that reduction process. This is read into ``pyXMIP`` and
    becomes an instance of the :py:class:`schema.ReductionSchema` class.

Notes
-----

For a complete guide on ``Schema`` classes, check out the reference page: :ref:`schema`.



"""
import os
import pathlib as pt
from copy import copy

import astropy.coordinates as astro_coords

from pyXMIP.utilities._registries import _Registry
from pyXMIP.utilities.core import _bin_directory, getFromDict, mainlog, setInDict

#: directory containing the built-in schema.
builtin_schema_directory = os.path.join(_bin_directory, "builtin_schema")
#: directory containing the built-in :py:class:`SourceTableSchema` instances.
builtin_source_table_schema_directory = os.path.join(
    builtin_schema_directory, "source_table"
)
#: directory containing the built-in :py:class:`ReductionSchema` instances.
builtin_reduction_schema_directory = os.path.join(builtin_schema_directory, "reduction")


class SchemaError(Exception):
    """
    Collective error type for issues with schema.
    """

    def __init__(self, message=None):
        self.message = message
        super().__init__(self.message)


class SchemaParameterError(SchemaError):
    """
    Exception raised when a parameter in a schema is invalid.
    """

    def __init__(self, parameter_name=None, instance=None, reason="NA"):
        self.message = f"PARAMETER {parameter_name} associated with SCHEMA {instance} is not valid. REASON: {reason}."
        super().__init__(self.message)


class SchemaEntry:
    """
    Descriptor class representing entries in a schema.
    """

    def __init__(
        self,
        dict_location,
        required=False,
        dtype=None,
        default=None,
        allowed_values=None,
    ):
        """
        Initialize the :py:class:`SchemaEntry` instance.

        Parameters
        ----------
        dict_location: list of str
            The location of this parameter in the underlying schema `dict`. This should **not include** the name of the instance. For
            example, if the parameter dictionary is ``{"Header1":{"example_a":0,"example_b":1}}``, then ```example_a`` should have ``dict_location``
            set to ``Header1``.
        required: callable, optional
            Function of the ``instance`` associated with the parameter which returns a bool representing if a parameter is actually required. By default, ``lambda x: False``.
        dtype: tuple or type or None
            The allowed data types for this parameter. By default ``None``, meaning no-restrictions.
        default: any
            The default value of this parameter if it is not provided. By default ``None``, meaning that it has no value.
        """
        self.dict_location = dict_location
        self.required = required
        self.dtype = dtype
        self.default = default
        self.allowed_values = allowed_values

    def __set_name__(self, owner, name):
        self._name = name
        self.dict_location += [self._name]

    def __get__(self, instance, owner):
        try:
            return getFromDict(instance, self.dict_location)
        except KeyError:
            return self.default
        except TypeError:
            return self.default

    def __set__(self, instance, value):
        setInDict(instance, self.dict_location, value)

        self.validate(instance)

    def validate(self, instance):
        """
        Validate that the instantiation of this setting in a given schema ``instance`` is valid.

        Parameters
        ----------
        instance: any
            The instance implementing this setting to check.

        Raises
        ------
        SchemaParameterError
            Raised if not valid.

        """
        value = getattr(instance, self._name)

        if value is not None:
            # the value is set, assure that it is set to something legitimate.
            if self.dtype is not None:
                if not isinstance(value, self.dtype):
                    raise SchemaParameterError(
                        parameter_name=self._name,
                        instance=instance,
                        reason=f"Type {type(value)} is not among {self.dtype}.",
                    )

            if self.allowed_values is not None:
                if value not in self.allowed_values:
                    raise SchemaParameterError(
                        parameter_name=self._name,
                        instance=instance,
                        reason=f"Value {value} is not in permitted values {self.allowed_values}.",
                    )

        else:
            # The value isn't set. Is it required?
            if self.required(instance):
                raise SchemaParameterError(
                    parameter_name=self._name,
                    instance=instance,
                    reason="Parameter is required for this instance, but not specified.",
                )

    @staticmethod
    def true(instance):
        """Callable to arbitrarily return ``True`` for required parameters."""
        return True

    @staticmethod
    def false(instance):
        """Callable to arbitrarily return ``False`` for required parameters."""
        return False

    @staticmethod
    def has(entry):
        def _func(instance):
            return getattr(instance, entry) is not None

        return _func

    @staticmethod
    def enabled(entry):
        def _func(instance):
            return getattr(instance, entry) is True

        return _func


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
    The generic schema class from which all other ``pyXMIP`` schemas are derived.

    .. warning::

        This is an abstract class.

    Notes
    -----

    Schema classes are wrappers of the built-in :py:class:`dict` type with additional structure. Generically,
    they can be read from file and checked for standardized structures.

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

    def __repr__(self):
        return f"<{self.__class__.__name__} instance>"

    @classmethod
    def _get_descriptors(cls):
        # fetching descriptors.
        return [m for m, v in cls.__dict__.items() if isinstance(v, SchemaEntry)]

    @classmethod
    def _get_required_descriptors(cls):
        return [
            m
            for m, v in cls.__dict__.items()
            if isinstance(v, SchemaEntry) and v.required is True
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
    Schema class for characterizing the structure of :py:class:`structures.table.SourceTable` objects.

    The :py:class:`SourceTableSchema` informs :py:class:`structures.table.SourceTable` of special roles for particular columns in the dataset.
    Additionally, the schema controls the native coordinate system, and the source type mapping from the native convention of the
    source table to the SIMBAD convention adopted in pyXMIP.


    Attributes
    ----------

    column_map:
        Mapping between special columns (keys) and their column names in the resulting tables (values). Only special columns
        will appear in the ``column_map``. It is the backbone of the table standardization approach used in pyXMIP.
    default_coord_system:
        The default coordinate system for this schema. Specificed as a string representation.
    object_map:
        Mapping between source table object types (keys) and pyXMIP convention [SIMBAD] object types (values).
    Z:
        The column name in the table which contains redshift information.
    TYPE:
        The column name for the object type.
    NAME:
        The column denoting the name of each object in the table.
    RA:
        The name of the column containing the RA coordinate.
    DEC:
        The name of the column containing the DEC coordinate.
    L:
        The name of the column containing galactic longitude.
    B:
        The name of the column containing galactic latitude.

    Notes
    -----

    **Schema Structure**:

    The :py:class:`SourceTableSchema` is composed of 3 components split across 3 headings in its ``yaml`` representation:

    1. ``column_map``: A dictionary mapping special / important table columns to their corresponding name in the table.
    2. ``object_map``: A dictionary mapping object types (star, etc.) to their SIMBAD / ``pyXMIP`` standard.
    3. ``settings``: A general purpose set of settings and other miscellaneous configuration options.

    In a ``.yaml`` file, the schema should have the following general layout:

    .. code-block:: yaml

        column_map:
            special_column: column_in_table
        object_map:
            type_in_table: standard_type
        settings:
            setting: value

    **Column Map**:

    The column map may contain any of the following special keys:

    +----------------------+-----------+-------------------------------------------------------------------------------+---------------------------------------------------------+
    | Special Column       | Required? | Description                                                                   | Identifiers                                             |
    +======================+===========+===============================================================================+=========================================================+
    | ``TYPE``             | ``False`` | Column containing the  object type of each entry (star, galaxy, etc.)         | TYPE, type, object_type, OBJECT_TYPE, OTYPES            |
    |                      |           | The column is not strictly necessary; however, any and all type reliant       |                                                         |
    |                      |           | processes will require that this is specified.                                |                                                         |
    +----------------------+-----------+-------------------------------------------------------------------------------+---------------------------------------------------------+
    | ``Z``                | ``False`` | Column containing the object redshifts. [Currently unused]                    | Z, z, REDSHIFT, redshift, BEST_Z                        |
    +----------------------+-----------+-------------------------------------------------------------------------------+---------------------------------------------------------+
    | ``NAME``             | ``True``  | Column containing the object identifier.                                      | NAME, name, OBJECT, object, ID, OBJECT_ID,              |
    |                      |           |                                                                               | IAUNAME, CATALOG_OBJECT                                 |
    +----------------------+-----------+-------------------------------------------------------------------------------+---------------------------------------------------------+
    | ``RA``               | ``False`` | Column containing the object coordinates.                                     | RA, ra                                                  |
    +----------------------+-----------+-------------------------------------------------------------------------------+---------------------------------------------------------+
    | ``DEC``              | ``False`` | Column containing the object coordinates.                                     | DEC, dec                                                |
    +----------------------+-----------+-------------------------------------------------------------------------------+---------------------------------------------------------+
    | ``B``                | ``False`` | Column containing the object coordinates.                                     | B, b, bii, BII                                          |
    +----------------------+-----------+-------------------------------------------------------------------------------+---------------------------------------------------------+
    | ``L``                | ``False`` | Column containing the object coordinates.                                     | l, lii, L, LII                                          |
    +----------------------+-----------+-------------------------------------------------------------------------------+---------------------------------------------------------+

    .. hint::

        Although there is no requirement for any given coordinate to be present, **at least 1** complete sky coordinate system must be present.

    **Object Map**:

    The ``object_map`` is not strictly required, but is necessary to do any object type conversions. This should be structured as a standard dictionary
    with the object-types present in the :py:class:`structures.table.SourceTable` as keys and the corresponding SIMBAD standard object types as values.

    Duplicate values are permitted for mapping, as are missing SIMBAD types. **Every** type in the table must be present.

    .. note::

        In some databases, such as SIMBAD, multiple object types may be present.

    **Settings**:

    The following table contains all of the available settings keys and values.

    +--------------------------+-----------+-------------------------------------------------------------------------------+
    | Setting                  | Required? | Description                                                                   |
    +==========================+===========+===============================================================================+
    | ``default_coord_system`` | ``True``  | The default coordinate system to use for the table. This must be the class    |
    |                          |           | name of an :py:class:`astropy.coordinates.GenericFrame` class. In order to be |
    |                          |           | permissible, the necessary columns must be present in ``column_map``.         |
    +--------------------------+-----------+-------------------------------------------------------------------------------+
    | ``default_angle_units``  | ``False`` | The default units for angular measures. If this is **already specified** in   |
    |                          |           | your data (say in the ``fits`` header), then this setting will only convert   |
    |                          |           | those angles to these units. If not, it is assumed these units are the ones   |
    |                          |           | present in the data. By default, ``'deg'``.                                   |
    +--------------------------+-----------+-------------------------------------------------------------------------------+
    | ``object_type_separator``| ``False`` | If the ``TYPE`` column for the table contains multiple types separated by a   |
    |                          |           | delimiter, specify that delimeter here. By default, it is assumed to be ``|`` |
    +--------------------------+-----------+-------------------------------------------------------------------------------+

    """
    _coordinate_requirements = {
        astro_coords.ICRS: ["RA", "DEC"],
        astro_coords.Galactic: ["L", "B"],
    }
    _expected_special_columns = {
        "Z": ["Z", "z", "REDSHIFT", "redshift", "BEST_Z"],
        "TYPE": ["TYPE", "type", "object_type", "OBJECT_TYPE", "OTYPES"],
        "NAME": [
            "NAME",
            "name",
            "OBJECT",
            "object",
            "ID",
            "OBJECT_ID",
            "IAUNAME",
            "CATALOG_OBJECT",
        ],
        "RA": ["RA", "ra"],
        "DEC": ["DEC", "dec"],
        "L": ["l", "lii", "L", "LII"],
        "B": ["B", "b", "bii", "BII"],
    }

    # -- Construct the attributes -- #
    column_map = SchemaEntry([], required=SchemaEntry.true, dtype=dict)
    default_coord_system = SchemaEntry(
        ["settings"], required=SchemaEntry.true, dtype=str
    )
    default_angle_units = SchemaEntry(
        ["settings"], required=SchemaEntry.false, dtype=str, default="deg"
    )
    object_type_separator = SchemaEntry(
        ["settings"], required=SchemaEntry.false, dtype=str, default="|"
    )
    object_map = SchemaEntry([], required=SchemaEntry.true, dtype=dict)

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
    """
    Schema class for representing reduction pipelines for cross-matching databases.
    """

    # ------------------------------------------------- #
    # Runtime parameter definitions                     #
    # ------------------------------------------------- #
    enable_astrometry = SchemaEntry(
        ["RUNTIME_PARAMS"], required=SchemaEntry.false, dtype=bool, default=True
    )
    """bool: ``True`` to enable instrumental reduction pipeline."""
    enable_poisson = SchemaEntry(
        ["RUNTIME_PARAMS"], required=SchemaEntry.false, dtype=bool, default=True
    )
    """bool: ``True`` to enable the Poisson reduction pipeline."""
    # -------------------------------------------------- #
    # Astrometry parameters                              #
    # -------------------------------------------------- #
    astrometry_weight = SchemaEntry(
        ["ASTROMETRY_PARAMS"],
        required=SchemaEntry.enabled("enable_astrometry"),
        dtype=bool,
        default=1,
    )
    """float: The weight of this process relative to others."""
    astrometry_mode_cat = SchemaEntry(
        ["ASTROMETRY_PARAMS"],
        required=SchemaEntry.enabled("enable_astrometry"),
        dtype=str,
        default="circular",
    )
    """str: The astrometric error mode to use for the catalog.

    Options are ``None``, ``circular`` or ``axial``. If ``circular``, a circular gaussian error is utilized, otherwise an axial (RA/DEC) gaussian is used.
    """
    astrometry_mode_db = SchemaEntry(
        ["ASTROMETRY_PARAMS"],
        required=SchemaEntry.enabled("enable_astrometry"),
        dtype=(str, dict),
        default="circular",
    )
    """dict or str: The astrometric error mode to use for the reference databases.

    Options are ``None``,``circular`` or ``axial``. If ``circular``, a circular gaussian error is utilized, otherwise an axial (RA/DEC) gaussian is used. If a ``str`` is provided,
    the same approach is used for all of the databases. Otherwise, dictionary should have format ``{table_name: mode}`` for each table in the cross-match database.
    """
    astrometry_err_cat = SchemaEntry(
        ["ASTROMETRY_PARAMS"],
        required=SchemaEntry.enabled("enable_astrometry"),
        dtype=None,
        default="circular",
    )
    r"""list: The error measure (1 :math:`\sigma`) for catalog sources.

    Must be formatted as a list. Each entry corresponds to **1 dimension** of the error. If :py:attr:`~ReductionSchema.astrometry_mode_cat` is ``'circular'``, then
    the list must have length 1; if ``'axial'``, then list is length 2 with RA error and DEC error respectively.

    For each error that needs to be specified, either a ``list`` or :py:class:`astropy.units.Quantity` may be supplied.

    - ``list``: 2-components: ``value``, ``unit``.

      - If ``value`` is a ``str``: Corresponds to a column in the ``CATALOG`` table of the cross-match database which contains the relevant error.
      - If ``value`` is ``numeric``: Assumed to be the actual value of the relevant error in the specified units.

    - :py:class:`astropy.units.Quantity`: Set both ``value`` and ``unit`` to those of the quantity.

    Examples
    --------

    .. code-block:: yaml

        ASTROMETRY_PARAMS:
          astrometry_mode_cat: 'circular'
          astrometry_err_cat:
            - ['POSERR','arcsec'] # --> Use the POSERR column for the 1-sigma circular precision.
          astrometry_mode_db: 'axial'
          astrometry_err_db:
            SIMBAD_STD_MATCH:
              - ['RA_ERR','deg'] # --> use column RA_ERR (assumed that units are deg).
              - ['DEC_ERR','arcmin'] # --> use column DEC_ERR (assumed that units are arcmin).

    """
    astrometry_err_db = SchemaEntry(
        ["ASTROMETRY_PARAMS"],
        required=SchemaEntry.enabled("enable_astrometry"),
        dtype=(str, dict),
        default="circular",
    )
    r"""dict or list: The error measure (1 :math:`\sigma`) for database sources.

    If the error parameters are the same for all relevant databases, then only one entry is required (``list``). Otherwise, a ``dict`` must be provided
    containing the table names and the relevant parameters.


    Each entry (for each database) must be formatted as a list. Each entry corresponds to **1 dimension** of the error.
    If :py:attr:`~ReductionSchema.astrometry_mode_cat` is ``'circular'``, then the list must have length 1; if ``'axial'``,
    then list is length 2 with RA error and DEC error respectively.

    For each error that needs to be specified, either a ``list`` or :py:class:`astropy.units.Quantity` may be supplied.

    - ``list``: 2-components: ``value``, ``unit``.

      - If ``value`` is a ``str``: Corresponds to a column in the ``CATALOG`` table of the cross-match database which contains the relevant error.
      - If ``value`` is ``numeric``: Assumed to be the actual value of the relevant error in the specified units.

    - :py:class:`astropy.units.Quantity`: Set both ``value`` and ``unit`` to those of the quantity.

    Examples
    --------

    .. code-block:: yaml

        ASTROMETRY_PARAMS:
          astrometry_mode_cat: 'circular'
          astrometry_err_cat:
            - ['POSERR','arcsec'] # --> Use the POSERR column for the 1-sigma circular precision.
          astrometry_mode_db:
            SIMBAD_STD_MATCH: 'axial'
            NED_STD_MATCH: 'circular'
          astrometry_err_db:
            SIMBAD_STD_MATCH:
              - ['RA_ERR','deg'] # --> use column RA_ERR (assumed that units are deg).
              - ['DEC_ERR','arcmin'] # --> use column DEC_ERR (assumed that units are arcmin).
            NED_STD_MATCH:
              - ['ASTERR','deg'] # --> circular error provided by the NED_STD ASTERR column.

    """
    # -------------------------------------------------- #
    # Poisson Parameters                                 #
    # -------------------------------------------------- #
    poisson_weight = SchemaEntry(
        ["POISSON_PARAMS"],
        required=SchemaEntry.enabled("enable_poisson"),
        dtype=bool,
        default=1,
    )
    """float: The weight of this process relative to others."""


class SchemaRegistry(_Registry):
    """
    Registry class for containing collections of :py:class:`Schema` classes.
    """

    def __init__(self, mapping, schema_class):
        """
        Initialize the :py:class:`SchemaRegistry`

        Parameters
        ----------
        mapping: dict
            The dictionary of ``name``,``path`` pairs to use as the registry.
        schema_class: :py:class:`Schema`
            The schema class to initialize objects with.
        """
        super().__init__(mapping)

        self.schema_class = schema_class

    def __getitem__(self, item):
        return self.schema_class.from_file(super().__getitem__(item))

    @classmethod
    def from_directory(cls, directory, schema_class):
        """
        Generate a registry of :py:class:`Schema` instances from a directory of files.

        Parameters
        ----------
        directory: str
            The directory to load from.

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
