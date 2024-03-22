"""
Custom table module for specialized table types.
"""
from astropy.table import Table, TableAttribute, Column
from astropy.io import fits
from pyXMIP.schema import SourceTableSchema
from astropy.coordinates import SkyCoord
import os
import numpy as np

class SourceTable(Table):
    _schema = TableAttribute()

    @property
    def schema(self):
        if self._schema is None:
            self.generate_schema()
        return self._schema

    @property
    def lon(self):
        return self[self.schema.coordinate_columns[0]]

    @property
    def lat(self):
        return self[self.schema.coordinate_columns[1]]

    def get_coordinates(self):
        coord_unit = self.schema['settings'].get('coordinate_units',None) #! This is necessary when tables have units that aren't recognized.

        if coord_unit is None:
            return SkyCoord(self.lon,self.lat,frame=self.schema.coordinate_frame)
        else:
            return SkyCoord(list(self.lon.value),list(self.lat.value),frame=self.schema.coordinate_frame,unit=coord_unit)
    @property
    def redshift(self):
        if self.schema.Z is not None:
            return self[self.schema.Z]
        else:
            raise ValueError(f"The schema associated with this table doesn't have a mapping to a redshift column.")

    @property
    def name(self):
        if self.schema.NAME is not None:
            return self[self.schema.NAME]
        else:
            raise ValueError(f"The schema associated with this table doesn't have a mapping to a redshift column.")

    @property
    def type(self):
        if self.schema.TYPE is not None:
            return self[self.schema.TYPE]
        else:
            raise ValueError(f"The schema associated with this table doesn't have a mapping to a redshift column.")

    @property
    def type_native(self):
        _base = self.type
        return  Column(data = [self.schema['object_map'][l] for l in _base],name='Native Type')

    @schema.setter
    def schema(self,value):
        self._schema = value

    def generate_schema(self):
        self._schema = SourceTableSchema.construct(self)

    def count_types(self):
        return Table({
            k:[len(self[self[self.schema.TYPE] == k])] for k in self.schema['object_map']
        })

    def _enforce_fits_compliant_dtypes(self):
        # -- remove schema and desc headers -- #
        self.meta = {}

        # -- Fixing column datatypes -- #
        for col in self.columns:
            if self[col].dtype == 'object':
                self[col] = np.array([str(j) for j in self[col]],dtype='<U8')

            if self[col].unit == 'degrees':
                self[col].unit = 'deg'
                self[col].format = None

    def append_to_fits(self,path,hudl):
        self._enforce_fits_compliant_dtypes()
        _self_hudl = fits.table_to_hdu(self)


        with fits.open(path,'update') as hudl_list:
            if hudl in hudl_list:
                _hudl, _len_hudl = hudl_list[hudl], len(hudl_list[hudl].data)
                new_hudl = fits.BinTableHDU.from_columns(_hudl.columns,nrows=_len_hudl+len(self))
                for colname in hudl_list[hudl].columns.names:
                    new_hudl.data[colname][_len_hudl:] = _self_hudl.data[colname]

                del hudl_list[hudl]
            else:
                new_hudl = fits.table_to_hdu(self)

            new_hudl.name = hudl
            hudl_list.append(new_hudl)
            hudl_list.flush()




if __name__ == '__main__':
    from databases import NED
    from astropy.coordinates import SkyCoord
    from astropy import units
    # = NED.query_radius(SkyCoord(ra=1,dec=1,unit='deg'),2*units.arcmin)
    from astropy.io import fits

    with fits.open('test.fits') as hudl:
        print(hudl.info())
