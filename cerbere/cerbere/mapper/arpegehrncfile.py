#mapper arpege hr 
import netCDF4
import numpy
import ast
from datetime import datetime
from collections import OrderedDict
#from .ncfile import NCFile
from .. import READ_ONLY, WRITE_NEW
from cerbere.mapper.ncfile import NCFile
from cerbere.datamodel.grid import Grid
from ..datamodel.field import Field
from cerbereutils.science.wind import uv2dir
from ..datamodel.variable import Variable

VIRTUALFIELD_DESCR = {'wind_speed': 'wind speed',
                      'wind_direction': 'wind direction (where the wind is blowing)'
                      }
VIRTUALFIELD_STDNAME = {'wind_speed': 'wind_speed',
                      'wind_direction': 'wind_to_direction'
                      }
VIRTUALFIELD_UNITS = {
                'wind_speed': 'm s-1',
                'wind_direction': 'degrees'
                }
class arpegeHRNCFile(NCFile):
    
    def __init__(self, url=None, mode=READ_ONLY,ncformat='NETCDF3_CLASSIC', **kwargs):
        super(arpegeHRNCFile,self).__init__(url=url,mode=mode,ncformat=ncformat,depths=True,center_on_greenwhich=False,**kwargs)
        if self._handler is None:
                self._handler = self.get_handler()
        return  

    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature

        Excludes the geolocation fields in the returned list
        '''
        geophyvars = super(arpegeHRNCFile, self).get_fieldnames()
        geophyvars.extend(["wind_speed", "wind_direction"])
        #geophyvars.remove('height')   
        return geophyvars

        
    def read_values(self, fieldname, slices=None):
        if fieldname == "wind_speed":
            u10 = super(arpegeHRNCFile, self).read_values('u10m', slices)
            v10 = super(arpegeHRNCFile, self).read_values('v10m', slices)
            u = numpy.sqrt(u10 * u10 + v10 * v10)
            shaap=numpy.shape(u)
            u=numpy.reshape(u,(1,shaap[2],shaap[3]))
            return u
        elif fieldname == "wind_direction":
            u10 = super(arpegeHRNCFile, self).read_values('u10m', slices)
            v10 = super(arpegeHRNCFile, self).read_values('v10m', slices)
            tmp=uv2dir(u10, v10)
            shaap=numpy.shape(tmp)
            tmp=numpy.reshape(tmp,(1,shaap[2],shaap[3]))
            return tmp
        elif fieldname not in ['lat','lon','time','height']:
            tmp= super(arpegeHRNCFile, self).read_values(fieldname, slices)
            shaap=numpy.shape(tmp)
            tmp=numpy.reshape(tmp,(1,shaap[2],shaap[3]))
            return tmp
        else: #case geolocation fields
            return super(arpegeHRNCFile, self).read_values(fieldname, slices)
            
    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        if fieldname in ['wind_speed', 'wind_direction']:
            # create a virtual field
            variable = Variable(
                    shortname=fieldname,
                    description=VIRTUALFIELD_DESCR[fieldname],
                    authority=self.get_naming_authority(),
                    standardname=VIRTUALFIELD_STDNAME[fieldname]
                    )
            field = Field(
                    variable,
                    OrderedDict([('time', 1),
                                 #('z', 1),
                                 ('y', self.get_dimsize('y')),
                                 ('x', self.get_dimsize('x'))
                                 ]),
                    datatype=numpy.dtype(numpy.float32),
                    units=VIRTUALFIELD_UNITS[fieldname]
                    )
            field.attach_storage(self.get_field_handler(fieldname))
        else:
            field = NCFile.read_field(self, fieldname)
            vava=field.variable
            uni=field.units
            if 'z' in field.get_dimnames():
                field = Field(
                    vava,
                    OrderedDict([('time', 1),
                                 #('z', self.get_dimsize('z')),
                                 ('y', self.get_dimsize('y')),
                                 ('x', self.get_dimsize('x'))
                                 ]),
                    datatype=numpy.dtype(numpy.float32),
                    units=uni
                    )
                field.attach_storage(self.get_field_handler(fieldname))
        return field
        