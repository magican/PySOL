'''
Created on 6 dec. 2012

@author: jfpiolle
'''

import numpy
from collections import OrderedDict

from .ncfile import NCFile
from cerbereutils.science.wind import uv2dir
from ..datamodel.field import Field
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


class ECMWFNCFile(NCFile):
    '''
    storage class for ECMWF files
    '''
    def __init__(self, url=None, **kwargs):
        """
        """
        super(ECMWFNCFile, self).__init__(url=url, **kwargs)
        return

    def get_fieldnames(self):
        # add two virtual variables for wind speed and direction
        fieldnames = super(ECMWFNCFile, self).get_fieldnames()
        fieldnames.extend(["wind_speed", "wind_direction"])
        return fieldnames

    def read_values(self, fieldname, slices=None):
        if fieldname == "wind_speed":
            u10 = super(ECMWFNCFile, self).read_values('10u', slices)
            v10 = super(ECMWFNCFile, self).read_values('10v', slices)
            u = numpy.sqrt(u10 * u10 + v10 * v10)
            return u
        elif fieldname == "wind_direction":
            u10 = super(ECMWFNCFile, self).read_values('10u', slices)
            v10 = super(ECMWFNCFile, self).read_values('10v', slices)
            return uv2dir(u10, v10)
        else:
            return super(ECMWFNCFile, self).read_values(fieldname, slices)

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
                                 ('y', self.get_dimsize('y')),
                                 ('x', self.get_dimsize('x'))
                                 ]),
                    datatype=numpy.dtype(numpy.float32),
                    units=VIRTUALFIELD_UNITS[fieldname]
                    )
            field.attach_storage(self.get_field_handler(fieldname))
        else:
            field = NCFile.read_field(self, fieldname)
        return field
