"""
Module for managing tabular data in ``pyXMIP``.

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

import numpy as np
import pandas as pd
from astropy import units as units
from astropy.coordinates import ICRS, Galactic, SkyCoord
from astropy.io import fits
from astropy.table import Column, Table, TableAttribute

from pyXMIP.schema import SourceTableSchema
from pyXMIP.utilities.core import enforce_units


class _SchemaColumn:
    def __set_name__(self, owner, name):
        self._name = name

    def __get__(self, instance, owner):
        try:
            return instance[getattr(instance.schema, self._name)]
        except KeyError:
            raise ValueError(f"Schema has no {self._name} column identified.")

    def __set__(self, instance, value):
        setattr(instance.schema, self._name, value)


class _CoordinateColumn(_SchemaColumn):
    def __init__(self, latlon, frame):
        self.latlon = latlon
        self.frame = frame

    def __get__(self, instance, owner):
        positions = instance.get_coordinates()
        positions = positions.transform_to(self.frame)

        return getattr(positions.frame.spherical, self.latlon)


class SourceTable(Table):
    """
    ``pyXMIP`` specific form of :py:class:`astropy.table.table.Table` for representing source catalogs and interacting
    with :py:class:`schema.SourceTableSchema`.

    Attributes
    ----------
    TYPE:
        The object type column of the table.
    NAME:
        The name column of the table.
    Z:
        The redshift column (if it exists).
    RA:
        The RA of the table objects.
    DEC:
        The DEC of the table objects.
    L:
        The galactic longitude of the table objects.
    B:
        The galactic latitude of the table objects.
    """

    _schema = TableAttribute()
    TYPE = _SchemaColumn()
    NAME = _SchemaColumn()
    Z = _SchemaColumn()
    RA = _CoordinateColumn("lon", ICRS)
    DEC = _CoordinateColumn("lat", ICRS)
    L = _CoordinateColumn("lon", Galactic)
    B = _CoordinateColumn("lat", Galactic)

    @property
    def schema(self):
        """
        The :py:class:`schema.SourceTableSchema` associated with this table.
        """
        if self._schema is None:
            self.generate_schema()
        return self._schema

    @property
    def lon(self):
        """The native longitude (default coordinate frame) of the catalog."""
        return enforce_units(
            units.Quantity(self[self.schema.coordinate_columns[0]]),
            self.schema.default_angle_units,
        )

    @property
    def lat(self):
        """The native latitude (default coordinate frame) of the catalog."""
        return enforce_units(
            units.Quantity(self[self.schema.coordinate_columns[1]]),
            self.schema.default_angle_units,
        )

    @property
    def type_native(self):
        """TYPE column converted to ``pyXMIP`` standard (SIMBAD standard) conventions."""
        _base = self.TYPE
        return Column(
            data=[self.schema["object_map"][l] for l in _base], name="Native Type"
        )

    @schema.setter
    def schema(self, value):
        self._schema = value

    def get_coordinates(self):
        """
        Obtain the default coordinate positions for the catalog objects.

        Returns
        -------
        :py:class:`astropy.coordinates.SkyCoord`
            The coordinates of the souces in default frame.
        """
        return SkyCoord(self.lon, self.lat, frame=self.schema.coordinate_frame)

    def generate_schema(self):
        """
        Automatically generate a schema for this table.

        .. hint::

            This is also done automatically by simply calling ``self.schema``.

        Returns
        -------
        None

        """
        self._schema = SourceTableSchema.construct(self)

    def get_formatted_types(self):
        # grab the separator for the type column
        _sep = self.schema.object_type_separator

        # iterate through and correct.
        _types = list(self.TYPE)

        for i, val in enumerate(_types):
            if not len(val):
                _types[i] = f"{_sep}?{_sep}"
                continue

            if val[0] != _sep:
                _types[i] = _sep + _types[i]

            if val[-1] != _sep:
                _types[i] = _types[i] + _sep

        return _types

    def count_types(self):
        """
        Count the number of each object type in the catalog.

        Returns
        -------

        """
        assert (
            self.schema.TYPE in self.columns
        ), f"Cannot count types because there is no TYPE column {self.schema.TYPE}."

        # ------------------------------------------------- #
        # Correcting the formatting
        # ------------------------------------------------- #
        if len(self) == 0:
            return Table({k: [0] for k in self.schema["object_map"]})

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
                for k in self.schema["object_map"]
            }
        )

    def correct_type_formatting(self):
        self.TYPE = pd.Series(self.get_formatted_types())

    def append_to_fits(self, path, hudl):
        _self_hudl = fits.table_to_hdu(self)

        try:
            with fits.open(path, "update") as hudl_list:
                if hudl in hudl_list:
                    _hudl, _len_hudl = hudl_list[hudl], len(hudl_list[hudl].data)
                    new_hudl = fits.BinTableHDU.from_columns(
                        _hudl.columns, nrows=_len_hudl + len(self)
                    )
                    for colname in hudl_list[hudl].columns.names:
                        new_hudl.data[colname][_len_hudl:] = _self_hudl.data[colname]

                    del hudl_list[hudl]
                else:
                    new_hudl = fits.table_to_hdu(self)

                new_hudl.name = hudl
                hudl_list.append(new_hudl)
                hudl_list.flush()
        except Warning:
            raise ValueError()

    def append_to_sql(self, table, engine):
        with engine.connect() as conn:
            self.to_pandas().to_sql(table, con=conn, if_exists="append", index=False)

    @classmethod
    def _convert_types(cls, table, schema):
        # -- omap -- #
        omap = schema.object_map

        # -- pull the table's types -- #
        _sep = schema.object_type_separator
        _types = pd.Series(table[schema.TYPE])
        _types = _types.str[1:-1].str.split(_sep)  # --> convert to the string subsets.

        # -- Fixing the table types -- #
        _types = _types.apply(
            lambda x: "|" + "|".join([omap.get(l, f"NSM_{l}") for l in x]) + "|"
        )

        table[schema.TYPE] = _types

        return table


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
