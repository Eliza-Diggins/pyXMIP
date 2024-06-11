"""
Module containing core data structures and validation procedures for interacting with user data.

Notes
-----

``pyXMIP`` uses the ``pydantic`` library to validate user data. To do so, all of the data structures used in
the package require some boilerplate code to provide for validation. To do this without losing the readability of the
main package, this module contains all of that boilerplate code.

"""
from typing import Any, Callable, Mapping, Type

import astropy.coordinates as astro_coords
import numpy as np
from astropy.units import Quantity, Unit
from pydantic import (
    AfterValidator,
    BaseModel,
    ConfigDict,
    Field,
    GetCoreSchemaHandler,
    GetJsonSchemaHandler,
    model_validator,
)
from pydantic.json_schema import JsonSchemaValue
from pydantic_core import PydanticUndefinedType, core_schema
from sqlalchemy.types import (
    BOOLEAN,
    DATETIME,
    FLOAT,
    INTEGER,
    NULLTYPE,
    TEXT,
    Interval,
    TypeEngine,
)
from typing_extensions import Annotated, Literal

from pyXMIP.utilities.logging import mainlog

# -- Importing SELF -- #
try:
    from typing import Self  # noqa
except ImportError:
    from typing_extensions import Self as Self  # noqa


def convert_np_type_to_sql(dtype: np.dtype) -> Type[TypeEngine]:
    """
    Convert a standard ``numpy`` data type into its corresponding ``sqlalchemy`` dtype.

    Parameters
    ----------
    dtype: :py:class:`numpy.dtype`
        The data type (numpy) to convert.

    Returns
    -------
    TypeEngine
        The class representing that data type in sqlalchemy.

    """

    _kind = dtype.kind

    if _kind == "f":
        return FLOAT
    elif _kind in ["i", "u"]:
        return INTEGER
    elif _kind in ["m"]:
        return Interval
    elif _kind in ["M"]:
        return DATETIME
    elif _kind in ["U", "S", "O"]:
        return TEXT
    elif _kind in ["b"]:
        return BOOLEAN
    else:
        return NULLTYPE


# ============================================================= #
# Custom Annotations
# ============================================================= #
class _UnitTypePydanticAnnotation:
    """
    This is a PyDantic annotation to process unit-typed objects.
    """

    @classmethod
    def __get_pydantic_core_schema__(
        cls,
        _source_type: Any,
        _handler: GetCoreSchemaHandler,
    ) -> core_schema.CoreSchema:
        """
        We return a pydantic_core.CoreSchema that behaves in the following ways:

        * strings will be parsed as `Unit` instances with the int as the x attribute
        * `Unit` instances will be parsed as `Unit` instances without any changes
        * Nothing else will pass validation
        * Serialization will always return just the string
        """

        def validate_from_str(value: str) -> Unit:
            result = Unit(value)
            return result

        from_str_schema = core_schema.chain_schema(
            [
                core_schema.str_schema(),
                core_schema.no_info_plain_validator_function(validate_from_str),
            ]
        )

        return core_schema.json_or_python_schema(
            json_schema=from_str_schema,
            python_schema=core_schema.union_schema(
                [
                    # check if it's an instance first before doing any further work
                    core_schema.is_instance_schema(Unit),
                    from_str_schema,
                ]
            ),
            serialization=core_schema.plain_serializer_function_ser_schema(
                lambda instance: str(instance)
            ),
        )

    @classmethod
    def __get_pydantic_json_schema__(
        cls, _core_schema: core_schema.CoreSchema, handler: GetJsonSchemaHandler
    ) -> JsonSchemaValue:
        # Use the same schema that would be used for `str`
        return handler(core_schema.str_schema())


class _QuantityTypePydanticAnnotation:
    """
    This is a PyDantic annotation to process quantity-typed objects.
    """

    @classmethod
    def __get_pydantic_core_schema__(
        cls,
        _source_type: Any,
        _handler: GetCoreSchemaHandler,
    ) -> core_schema.CoreSchema:
        """
        We return a pydantic_core.CoreSchema that behaves in the following ways:

        * strings will be parsed as `Unit` instances with the int as the x attribute
        * `Unit` instances will be parsed as `Unit` instances without any changes
        * Nothing else will pass validation
        * Serialization will always return just the string
        """

        def validate_from_str(value: str) -> Quantity:
            result = Quantity(value)
            return result

        from_str_schema = core_schema.chain_schema(
            [
                core_schema.str_schema(),
                core_schema.no_info_plain_validator_function(validate_from_str),
            ]
        )

        return core_schema.json_or_python_schema(
            json_schema=from_str_schema,
            python_schema=core_schema.union_schema(
                [
                    # check if it's an instance first before doing any further work
                    core_schema.is_instance_schema(Quantity),
                    from_str_schema,
                ]
            ),
            serialization=core_schema.plain_serializer_function_ser_schema(
                lambda instance: str(instance)
            ),
        )

    @classmethod
    def __get_pydantic_json_schema__(
        cls, _core_schema: core_schema.CoreSchema, handler: GetJsonSchemaHandler
    ) -> JsonSchemaValue:
        # Use the same schema that would be used for `str`
        return handler(core_schema.str_schema())


class _FrameTypePydanticAnnotation:
    """
    This is a PyDantic annotation to process coordinate frame objects.
    """

    @classmethod
    def __get_pydantic_core_schema__(
        cls,
        _source_type: Any,
        _handler: GetCoreSchemaHandler,
    ) -> core_schema.CoreSchema:
        def validate_from_str(value: str) -> astro_coords.GenericFrame:
            result = getattr(astro_coords, value)

            assert type(result) is type(
                astro_coords.GenericFrame
            ), f"Astropy coordinates {value} attribute is not a frame."
            return result

        from_str_schema = core_schema.chain_schema(
            [
                core_schema.str_schema(),
                core_schema.no_info_plain_validator_function(validate_from_str),
            ]
        )

        return core_schema.json_or_python_schema(
            json_schema=from_str_schema,
            python_schema=core_schema.union_schema(
                [
                    # check if it's an instance first before doing any further work
                    core_schema.is_instance_schema(astro_coords.GenericFrame),
                    from_str_schema,
                ]
            ),
            serialization=core_schema.plain_serializer_function_ser_schema(
                lambda instance: instance.__name__
            ),
        )

    @classmethod
    def __get_pydantic_json_schema__(
        cls, _core_schema: core_schema.CoreSchema, handler: GetJsonSchemaHandler
    ) -> JsonSchemaValue:
        # Use the same schema that would be used for `str`
        return handler(core_schema.str_schema())


class _CMDTypePydanticAnnotation:
    """
    This is a PyDantic annotation to process :py:class:`CrossMatchDatabase` objects.
    """

    @classmethod
    def __get_pydantic_core_schema__(
        cls,
        _source_type: Any,
        _handler: GetCoreSchemaHandler,
    ) -> core_schema.CoreSchema:
        from pyXMIP.cross_reference import CrossMatchDatabase

        def validate_from_str(value: str) -> CrossMatchDatabase:
            result = CrossMatchDatabase(value)
            return result

        from_str_schema = core_schema.chain_schema(
            [
                core_schema.str_schema(),
                core_schema.no_info_plain_validator_function(validate_from_str),
            ]
        )

        return core_schema.json_or_python_schema(
            json_schema=from_str_schema,
            python_schema=core_schema.union_schema(
                [
                    # check if it's an instance first before doing any further work
                    core_schema.is_instance_schema(CrossMatchDatabase),
                    from_str_schema,
                ]
            ),
            serialization=core_schema.plain_serializer_function_ser_schema(
                lambda instance: instance.__name__
            ),
        )

    @classmethod
    def __get_pydantic_json_schema__(
        cls, _core_schema: core_schema.CoreSchema, handler: GetJsonSchemaHandler
    ) -> JsonSchemaValue:
        # Use the same schema that would be used for `str`
        return handler(core_schema.str_schema())


# ============================================================= #
# Custom Validators
# ============================================================= #
PydanticUnit = Annotated[Unit, _UnitTypePydanticAnnotation]
# Any: The Pydantic type representing the :py:class:`astropy.units.core.Unit` class.

PydanticQuantity = Annotated[Quantity, _QuantityTypePydanticAnnotation]
# Any: The Pydantic type representing the :py:class:`astropy.units.core.Quantity` class.

PydanticFrame = Annotated[astro_coords.GenericFrame, _FrameTypePydanticAnnotation]
# Any: The Pydantic type representing the :py:class:`astropy.units.core.Quantity` class.


def construct_template(datamodel: Type[BaseModel]) -> dict:
    """
    Construct a template version of a known data model. This is largely useless for external users,
    but is used in the backend to produce model schema which can then be filled in with information from objects like tables.

    Parameters
    ----------
    datamodel: BaseModel
        The ``pydantic`` model to produce a template version of.

    Returns
    -------
    dict
        The output template.
    """
    structure = {
        k: None for k, v in datamodel.model_fields.items()
    }  # Grab all of the model fields.

    for k, v in datamodel.model_fields.items():
        if type(v.default) is not PydanticUndefinedType and v.default is not None:
            structure[k] = v.default
        else:
            _underlying_type = v.annotation

            if type(_underlying_type) is type(BaseModel):
                structure[k] = construct_template(_underlying_type)
            else:
                pass

    return structure


# ============================================================= #
# Custom Structures
# ============================================================= #
#   Custom structures are Pydantic recognized data structures for validation
#   and eventual usage in schema objects.


class Registry(dict):
    """
    Bare-bones registry class.
    """

    def __init__(self, mapping: Mapping):
        """
        Initialize the generalized registry class.

        Parameters
        ----------
        mapping: dict
            The dictionary underlying this registry.
        """
        super().__init__(mapping)

    def as_list(self) -> list:
        return list(self.values())

    def add(self, key: str, value: Any):
        assert key not in self, f"Key {key} already in registry."
        self[key] = value


class TableColumn(BaseModel):
    """
    Pydantic data class representing a single column in a table.

    The :py:class:`TableColumn` class is utilized by :py:class:`ColumnMap` to construct the column structure for
    :py:class:`schema.SourceTableSchema`.
    """

    # -- Managing configuration -- #
    model_config = ConfigDict(validate_assignment=True)

    name: str
    """str: The column name corresponding to this column."""
    unit: PydanticUnit | None = None
    """The unit to attribute to this column."""


class SourceTableSchemaSettings(BaseModel):
    """
    Pydantic data class representing the settings for a :py:class:`schema.SourceTableSchema`.
    """

    # -- Managing configuration -- #
    model_config = ConfigDict(validate_assignment=True)

    # ============================= #
    # Settings                      #
    # ============================= #
    default_coord_system: Literal["ICRS", "Galactic"] = "ICRS"
    """str: The default coordinate system to use.

    Use this to select the native coordinate system for the schema. Whichever coordinate system is selected, there must
    be sufficient identifiable columns to actually obtain those coordinates.

    """
    object_type_separator: str = Field(default="|", min_length=1, max_length=1)
    """str: The separator character / string between different object types.

        In many databases, returned object types may have multiple types separated by a delimiter. This delimiter should
        be specified here if it is needed.

        .. important::

            In all native ``pyXMIP`` data objects, a "|" is used for the delimiter. Even if the object type is a single
            value, all object types in ``pyXMIP`` tables are formatted ``"|entry_1|entry_2|...|"``.
    """


class ColumnMap(BaseModel):
    """
    Pydantic data class representing a mapping from columns in a table to standard pyxmip columns.

    The :py:class:`ColumnMap` class is used to construct the column structure for :py:class:`schema.SourceTableSchema`.
    """

    # ========================================================================== #
    # Validators                                                                 #
    # ========================================================================== #
    # -- Managing configuration -- #
    model_config = ConfigDict(validate_assignment=True)

    @staticmethod
    def validate_table_column_unit(
        expected_unit: Unit, equivalencies: list[tuple] | None = None
    ):
        """
        PyDactic validator method to enforce unit norms on :py:class:`TableColumn` instances in the :py:class:`ColumnMap`.

        Parameters
        ----------
        expected_unit: Unit
            The expected unit. The unit specified in the data must be convertible to this unit.
        equivalencies
            The equivalencies to allow when checking for type consistency. By default, ``None``.

        Raises
        ------
        ValidationError
            [UNIT CONSISTENCY FAULT] Failed to convert ``value`` to units ``expected_unit``.

        """

        def _check_unit_validation(value: TableColumn) -> TableColumn:
            # -- Check if NONE -- #
            # If the unit isn't provided, we can just assign the expected unit.
            if value.unit is None:
                value.unit = expected_unit

            # -- Check unit type -- #
            # The user might have specified a string, convert to unit.
            if isinstance(value.unit, str):
                value.unit = Unit(value.unit)

            try:
                _ = value.unit.to(expected_unit, equivalencies=equivalencies)
            except Exception:
                raise AssertionError(
                    f"[UNIT CONSISTENCY FAULT] Failed to convert {value} to units {expected_unit}. Raised by ColumnMap.validate_table_column_unit."
                )
            return value

        return _check_unit_validation

    @staticmethod
    def validate_table_column_no_unit(value: TableColumn) -> TableColumn:
        """
        PyDactic validator method to enforce unit norms on :py:class:`TableColumn` instances in the :py:class:`ColumnMap`.

        Parameters
        ----------
        value
            The :py:class:`TableColumn` being validated.


        Raises
        ------
        ValidationError
            [UNIT CONSISTENCY FAULT] Failed to convert ``value`` to units ``expected_unit``.

        """
        assert (
            value.unit is None
        ), f"[UNIT CONSISTENCY FAULT] This TableColumn should not have units but is {value}. Raised by ColumnMap.validate_table_column_no_unit."
        return value

    # Coordinate Columns #
    RA: Annotated[
        TableColumn, AfterValidator(validate_table_column_unit(Unit("deg")))
    ] | None = None
    """The column containing the object RA coordinate."""
    DEC: Annotated[
        TableColumn, AfterValidator(validate_table_column_unit(Unit("deg")))
    ] | None = None
    """The column containing the object DEC coordinate."""
    RA_ERR: Annotated[
        TableColumn, AfterValidator(validate_table_column_unit(Unit("deg")))
    ] | None = None
    """The column containing the object RA coordinate error."""
    DEC_ERR: Annotated[
        TableColumn, AfterValidator(validate_table_column_unit(Unit("deg")))
    ] | None = None
    """The column containing the object DEC coordinate error."""
    GAL_L: Annotated[
        TableColumn, AfterValidator(validate_table_column_unit(Unit("deg")))
    ] | None = None
    """The column containing the object's galactic longitude."""
    GAL_B: Annotated[
        TableColumn, AfterValidator(validate_table_column_unit(Unit("deg")))
    ] | None = None
    """The column containing the object's galactic latitude."""
    GAL_L_ERR: Annotated[
        TableColumn, AfterValidator(validate_table_column_unit(Unit("deg")))
    ] | None = None
    """The column containing the object's galactic longitude error."""
    GAL_B_ERR: Annotated[
        TableColumn, AfterValidator(validate_table_column_unit(Unit("deg")))
    ] | None = None
    """The column containing the object's galactic latitude error."""
    SEP: Annotated[
        TableColumn, AfterValidator(validate_table_column_unit(Unit("deg")))
    ] | None = None
    # Name columns #
    NAME: Annotated[TableColumn, AfterValidator(validate_table_column_no_unit)]
    """The column containing the object's main identifier."""
    # Extra Columns #
    Z: Annotated[
        TableColumn, AfterValidator(validate_table_column_no_unit)
    ] | None = None
    """The column containing the object's redshift."""
    TYPE: Annotated[
        TableColumn, AfterValidator(validate_table_column_no_unit)
    ] | None = None
    """The column containing the object type of the catalog object."""

    def search_inverse(self, key):
        return {v: k for k, v in self.model_fields.items()}[key]


class ColumnStdErrorSpecifier(BaseModel):
    """
    Pydantic data class specifying the error on a given column of a table.

    .. note::

        This class is table-agnostic, it won't know anything about the data it's describing until runtime. Thus,
        failure to create the class accurately may lead to late appearing validation errors because i.e. your table doesn't
        have the specified column.

    Notes
    -----

    The :py:class:`ColumnErrorSpecifier` should be understood as an "archetypal" structure for any reduction process where
    the precision in a particular field of the matching database is used. Users must specify a ``base_column``, corresponding to
    the column (i.e. RA) for which this specifier is determining the error (i.e. RA_PREC). Error can either by dynamic, as determined by
    providing an ``error_column``, or by providing an ``error_value``.
    """

    base_column: str | None = None
    """:py:class:`TableColumn`: The name of the column who's error is being specified.

    This should be the name of a column in the table which is being described by this data object.
    """
    error_column: TableColumn | None = None
    """:py:class:`TableColumn` optional: The column which specifies the error in :py:attr:`ColumnErrorSpecifier.base_column`

    At least one of :py:attr:`ColumnErrorSpecifier.error_column` and :py:attr:`ColumnErrorSpecifier.error_value` is required
    for validation to pass without errors.
    """
    error_value: PydanticQuantity | None = None
    """:py:class:`astropy.units.core.Quantity`, optional: The global value of the error in :py:attr:`ColumnErrorSpecifier.base_column`

    .. hint::

        If a table doesn't have different errors for different entries but instead has a flat error profile, this is the
        option you should use.

    At least one of :py:attr:`ColumnErrorSpecifier.error_column` and :py:attr:`ColumnErrorSpecifier.error_value` is required
    for validation to pass without errors.
    """
    error_convert_to_std: Callable = lambda x: x
    """Callable, optional: Optional parameter to provide a function to convert error to standard deviation.

    This class **assumes** that the provided error is the standard deviation of a normally distributed error. If the available
    error column is not the standard deviation, specify this function to take ``x`` (the uncorrected error) as an argument and
    return the standard deviation form of the error.
    """

    def validate_on_table(self, table):
        """
        Validate this column error specifier against a table.

        Parameters
        ----------
        table: :py:class:`structures.table.SourceTable`
            The table to check this data structure against.

        Notes
        -----

        This validator performs 3 checks:

        1. Assert that ``self.base_column`` is actually a column in the table.
        2. If ``self.error_column`` is specified, then it must also be in the table.
        3. At least one of ``self.error_column`` and ``self.error_value`` must be specified.

        """

        if self.error_column is not None:
            assert (
                self.error_column.name in table.columns
            ), f"The base column {self.error_column.name} is not present in the validating table."

        assert not (
            self.error_column is None and self.error_value is None
        ), "At least one of error_column or error_value must be specified."

    @model_validator(mode="after")
    def validator(self):
        """
        Validate this column error specifier.


        Notes
        -----

        This validator performs 1 check:

        1. At least one of ``self.error_column`` and ``self.error_value`` must be specified.

        """
        # Check that at least one of the needed columns is specified.
        assert not (
            self.error_column is None and self.error_value is None
        ), "At least one of error_column or error_value must be specified."

        return self


class ICRSCoordinateStdErrorSpecifier(BaseModel):
    """
    Pydantic data class specifying the astrometric errors in a table.

    Notes
    -----

    .. hint::

        Like many of the structures in :py:mod:`utilities.types`, this class is composed of various subtypes. When documentation
        says that something must have the "format" of another class, that implies that it may **either** be an instance of the
        actual class, or the ``dict`` provided must be interpretable as that class.

    To initialize a :py:class:`CoordinateErrorSpecifier`, you need to provide the following:

    - :py:attr:`CoordinateErrorSpecifier.longitude_column`: The column in the table representing the longitude (i.e. RA).
    - :py:attr:`CoordinateErrorSpecifier.latitude_column`: The column in the table representing the latitude (i.e. DEC).

    Both of these should have the native format of :py:class:`TableColumn`.

    .. warning::

        The two coordinate columns must be from the same coordinate system as ``pyXMIP`` treats them as generic lat/lon instead
        of trying to interpret



    """

    # Error specifiers #
    lon_error: TableColumn | PydanticQuantity | None = None
    """:py:class:`TableColumn` or :py:class:`astropy.units.core.Quantity`: The error specification in the longitude."""
    lat_error: TableColumn | PydanticQuantity | None = None
    """:py:class:`TableColumn` or :py:class:`astropy.units.core.Quantity`: The error specification in the latitude."""
    position_error: TableColumn | PydanticQuantity | None = None
    """:py:class:`TableColumn` or :py:class:`astropy.units.core.Quantity`: The circular error specification."""

    @property
    def mode(self):
        if self.position_error is not None:
            return "circular"
        elif all([cerr is not None for cerr in [self.lat_error, self.lon_error]]):
            return "axial"
        else:
            return None

    @model_validator(mode="after")
    def validator(self):
        """
        Validate the :py:class:`CoordinateErrorSpecifier`.

        Notes
        -----

        The only check performed by this validator is that either the lat/lon error is specified or the circular is specified.
        """
        if all(
            [i is None for i in [self.lat_error, self.lon_error, self.position_error]]
        ):
            return self

        if self.position_error is None:
            assert not any(
                cerr is None for cerr in [self.lat_error, self.lon_error]
            ), f"CoordinateErrorSpecifier {self} is not circular but doesn't have 2 axes for position error."

        else:
            for cerr in [self.lat_error, self.lon_error]:
                if cerr is not None:
                    mainlog.warning(
                        f"CoordinateErrorSpecifier {self} is in MODE=circular, but axial errors are specified: {cerr}."
                    )

        return self

    def check_empty(self):
        assert any(
            [
                i is not None
                for i in [self.lat_error, self.lon_error, self.position_error]
            ]
        ), "At least one error specifier is needed."
