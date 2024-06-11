"""
Module for interacting with tabulated data in ``pyXMIP``.

Notes
-----

For general usage, :py:class:`SourceTable` operates exactly the same as :py:class:`astropy.table.table.Table` with only a few
key exceptions:

1. :py:class:`SourceTable` implements a :py:class:`schema.SourceTableSchema` on top of the existing astropy structure.
2. :py:class:`SourceTable` utilizes its schema to deduce additional information not inherently included in the table.

For general usage information, particularly regarding interacting with, slicing from, and altering tables, we recommend
reading the associated ``astropy`` documentation `here <https://docs.astropy.org/en/stable/table/index.html>`_.
"""
import re
from typing import Generic, Type, TypeVar

import numpy as np
import pandas as pd
from astropy import units as units
from astropy.coordinates import GenericFrame, SkyCoord
from astropy.coordinates.angles.core import Latitude, Longitude
from astropy.io import fits
from astropy.table import Column, Table, TableAttribute

from pyXMIP.schema import SourceTableSchema
from pyXMIP.utilities.core import enforce_units

Instance = TypeVar("Instance")
Value = TypeVar("Value")
Attribute = TypeVar("Attribute")


class SchemaColumn(Generic[Instance, Attribute, Value]):
    """
    Special column of the table specified by the schema.
    """

    def __set_name__(self, owner: Type[Instance], name: str) -> None:
        self._name = name

    def __get__(self, instance: Instance, owner: Type[Instance]) -> Column:
        try:
            return instance[getattr(instance.schema, self._name)]
        except (KeyError, ValueError):
            raise ValueError(f"Schema has no {self._name} column identified.")


class CoordinateColumn(SchemaColumn):
    def __init__(self, latlon: str, frame: GenericFrame):
        self.latlon: str = latlon
        self.frame: GenericFrame = frame

    def __get__(
        self, instance: Instance, owner: Type[Instance]
    ) -> Latitude | Longitude:
        positions = instance.get_coordinates()
        positions = positions.transform_to(self.frame)

        return getattr(positions.frame.spherical, self.latlon)


class SourceTable(Table):
    """
    ``pyXMIP`` specific form of :py:class:`astropy.table.table.Table` for representing source catalogs and interacting
    with :py:class:`schema.SourceTableSchema`.
    """

    _schema: SourceTableSchema = TableAttribute()

    # =========================================== #
    # Columns                                     #
    # =========================================== #
    # Coordinate Columns #
    RA: Column = SchemaColumn()
    """The column containing the object RA coordinate."""
    DEC: Column = SchemaColumn()
    """The column containing the object DEC coordinate."""
    RA_ERR: Column = SchemaColumn()
    """The column containing the object RA coordinate error."""
    DEC_ERR: Column = SchemaColumn()
    """The column containing the object DEC coordinate error."""
    GAL_L: Column = SchemaColumn()
    """The column containing the object's galactic longitude."""
    GAL_B: Column = SchemaColumn()
    """The column containing the object's galactic latitude."""
    GAL_L_ERR: Column = SchemaColumn()
    """The column containing the object's galactic longitude error."""
    GAL_B_ERR: Column = SchemaColumn()
    """The column containing the object's galactic latitude error."""
    SEP: Column = SchemaColumn()
    """The column containing the separation between the search position and the output object."""
    # Name columns #
    NAME: Column = SchemaColumn()
    """The column containing the object's main identifier."""
    # Extra Columns #
    Z: Column = SchemaColumn()
    """The column containing the object's redshift."""
    TYPE: Column = SchemaColumn()
    """The column containing the object type of the catalog object."""

    @property
    def schema(self) -> SourceTableSchema:
        """
        The :py:class:`schema.SourceTableSchema` associated with this table.
        """
        if self._schema is None:
            self.generate_schema()
        return self._schema

    @schema.setter
    def schema(self, value: SourceTableSchema) -> None:
        self._schema = value

    @property
    def lon(self) -> units.Quantity:
        """The native longitude (default coordinate frame) of the catalog.

        Notes
        -----

        This property accesses the :py:attr:`schema.SourceTableSchema.coordinate_columns` attribute to determine what the
        default coordinate columns are. It then accesses them and returns the corresponding column of the table.

        If the column has a ``.unit`` attribute, that will be utilized; otherwise, the column is assumed to have the
        default schema units.

        See Also
        --------

        lat, lon, get_coordinates

        """
        unit = self.schema.coordinate_columns[0].unit
        if unit is None:
            unit = self.schema.default_angle_units

        return enforce_units(
            units.Quantity(self[self.schema.coordinate_columns[0].name]),
            unit,
        )

    @property
    def lat(self) -> units.Quantity:
        """The native latitude (default coordinate frame) of the catalog.

        Notes
        -----

        This property accesses the :py:attr:`schema.SourceTableSchema.coordinate_columns` attribute to determine what the
        default coordinate columns are. It then accesses them and returns the corresponding column of the table.

        If the column has a ``.unit`` attribute, that will be utilized; otherwise, the column is assumed to have the
        default schema units.

        See Also
        --------

        lat, lon, get_coordinates

        """
        unit = self.schema.coordinate_columns[1].unit
        if unit is None:
            unit = self.schema.default_angle_units

        return enforce_units(
            units.Quantity(self[self.schema.coordinate_columns[1].name]),
            unit,
        )

    def get_coordinates(self) -> SkyCoord:
        """
        Obtain the default coordinate positions for the catalog objects.

        Returns
        -------
        :py:class:`astropy.coordinates.SkyCoord`
            The coordinates of the sources in default frame.

        See Also
        --------

        lat, lon, get_coordinates
        """
        return SkyCoord(self.lon, self.lat, frame=self.schema.coordinate_frame)

    def generate_schema(self) -> None:
        """
        Generate and assign a default / automated schema to this table.

        Notes
        -----

        See :py:meth:`schema.SourceTableSchema.construct`

        """
        self._schema = SourceTableSchema.from_table(self)

    def get_formatted_types(self):
        """
        Returns an altered version of the :py:attr:`SourceTable.TYPE` column with standardized formatting.

        Returns
        -------
        Column
            The corrected formatting.

        Notes
        -----

        What this method does is detect the :py:attr:`schema.SourceTableSchema.object_type_separator` in each row,
        then reformat so that all of the types have the same standard formatting.

        """
        return format_table_types(self, self.schema)

    def count_types(self) -> Table:
        """
        Count the number of each object type in the catalog.

        Returns
        -------
        Table
            The count table. Formatted as the type followed by the number of such types.
        """
        assert (
            self.schema.TYPE in self.columns
        ), f"Cannot count types because there is no TYPE column {self.schema.TYPE}."

        # ------------------------------------------------- #
        # Correcting the formatting
        # ------------------------------------------------- #
        if len(self) == 0:
            # This is operational redundancy. If the table has no entries, we can just return 0 for each object type.
            return Table({k: [0] for k in self.schema.object_map})

        # Collect the separator #
        _sep = self.schema.object_type_separator

        # -- Managing types for non-zero length -- #
        _types = pd.Series(
            self.get_formatted_types()
        )  # --> this does all the reformatting.
        return Table(
            {
                k: [
                    len(_types[_types.str.contains(f"{re.escape(f'{_sep}{k}{_sep}')}")])
                ]
                for k in self.schema.object_map
            }
        )

    def append_to_fits(self, path, hdu):
        """
        Append this table to an existing fits table.

        Parameters
        ----------
        path: str
            The path to the fits file to append to.
        hdu: str
            The name of the HDU to append to.

        """
        _self_hdu = fits.table_to_hdu(self)

        try:
            with fits.open(path, "update") as hdu_list:
                if hdu in hdu_list:
                    _hdul, _len_hdul = hdu_list[hdu], len(hdu_list[hdu].data)
                    new_hudl = fits.BinTableHDU.from_columns(
                        _hdul.columns, nrows=_len_hdul + len(self)
                    )
                    for colname in hdu_list[hdu].columns.names:
                        new_hudl.data[colname][_len_hdul:] = _self_hdu.data[colname]

                    del hdu_list[hdu]
                else:
                    new_hudl = fits.table_to_hdu(self)

                new_hudl.name = hdu
                hdu_list.append(new_hudl)
                hdu_list.flush()
        except Warning:
            raise ValueError()

    def append_to_sql(self, table, conn, **kwargs):
        self.to_pandas().to_sql(
            str(table).upper(), con=conn, if_exists="append", index=False, **kwargs
        )

    @classmethod
    def _convert_types(cls, table, schema):
        # -- omap -- #
        omap = schema.object_map

        # -- pull the table's types -- #
        _sep = schema.object_type_separator
        _types = pd.Series(format_table_types(table, schema))
        _types = _types.str[1:-1].str.split(_sep)  # --> convert to the string subsets.

        # -- Fixing the table types -- #
        _types = _types.apply(
            lambda x: "|" + "|".join([omap.get(l, f"NSM_{l}") for l in x]) + "|"
        )

        table[schema.TYPE] = _types

        return table


def format_table_types(table: Table, schema: SourceTableSchema = None) -> list:
    """
    Returns an altered version of the :py:attr:`SourceTable.TYPE` column with standardized formatting.

    Parameters
    ----------
    table: Table
        The table from which the types should be formatted.
    schema: SourceTableSchema, optional
        The :py:class:`schema.SourceTableSchema` associated with the table. If ``table`` is a :py:class:`SourceTable`, then
        this will be determined automatically. Otherwise, the schema must be provied.

    Returns
    -------
    Column
        The corrected formatting.

    Notes
    -----

    What this method does is detect the :py:attr:`schema.SourceTableSchema.object_type_separator` in each row,
    then reformat so that all of the types have the same standard formatting.

    """
    if schema is None:
        if isinstance(table, SourceTable):
            schema = table.schema
        else:
            raise ValueError(
                "The table you provided doesn't have a built-in schema. One must be provided explicitly."
            )

    # grab the separator for the type column
    _sep = schema.object_type_separator

    # iterate through and correct.
    _types = list(table[schema.TYPE])

    for i, val in enumerate(_types):
        if not len(val):
            # This is a non-specified type. We set it to the ? type.
            _types[i] = f"{_sep}?{_sep}"
            continue

        if val[0] != _sep:
            # Correct the first column to enforce the separator.
            _types[i] = _sep + _types[i]

        if val[-1] != _sep:
            # Correct the final column.
            _types[i] = _types[i] + _sep

    return _types


def load(path, *args, **kwargs):
    """
    Load a catalog into ``pyXMIP``.

    Parameters
    ----------
    path: str
        The path to the catalog to load.


    Returns
    -------
    :py:class:`SourceTable`
        The resulting table.
    """
    return SourceTable.read(path, *args, **kwargs)


def correct_column_types(table):
    # -- Fix object column types -- #
    for col in table.columns:
        if table[col].dtype == "object":
            table[col] = np.array([str(j) for j in table[col]], dtype="<U64")

    return table
