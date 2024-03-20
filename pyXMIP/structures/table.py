"""

"""
from astropy.table import Table, TableAttribute
from pyXMIP.schema import SkyCollectionSchema
class SourceTable(Table):
    _schema = TableAttribute()

    @property
    def schema(self):
        if self._schema is None:
            self.generate_schema()
        return self._schema

    @schema.setter
    def schema(self,value):
        self._schema = value

    def generate_schema(self):
        self._schema = SkyCollectionSchema.construct(self)

if __name__ == '__main__':

    q = SourceTable.read('/home/ediggins/pyROSITA_test/eRASS1_Hard.v1.0.fits',format='fits')
    print(q.schema)
    print(SourceTable.schema)