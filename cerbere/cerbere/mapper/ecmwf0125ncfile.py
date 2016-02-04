'''
Created on 6 dec. 2012

@author: jfpiolle
'''

import numpy
import logging
from collections import OrderedDict
from .ncfile import NCFile
try:
    from cerform.wind import uv2dir
except ImportError as err:
    try:
        from cerbereutils.science.wind import uv2dir
        logging.warning('cerbereutils has been split into multiple packages '
                        'and no longer exists. Please add the new cer* '
                        'packages to your python environment.')
    except ImportError:
        raise err
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


class ECMWF0125NCFile(NCFile):
    '''
    storage class for ECMWF files
    '''
    def __init__(self, url=None, **kwargs):
        """
        """
        super(ECMWF0125NCFile, self).__init__(url=url, **kwargs)
        return

    def get_fieldnames(self):
        # add two virtual variables for wind speed and direction
        fieldnames = super(ECMWF0125NCFile, self).get_fieldnames()
        fieldnames.extend(["wind_speed", "wind_direction"])
        return fieldnames

    def read_values(self, fieldname, slices=None):
        '''return the values of a given geophysical field
        the height dimensions is removed
        here we considere that the user doesnt have to know the netcdf dimensions
        and therefore would provide slice as a dictionnary standard x:slice y:slice'''
        
        if slices != None:
            if isinstance(slices,list) and len(slices) == 3:
                slices.insert(0,slice(0,1,1))
            elif isinstance(slices,dict):
                slices['z'] = slice(0,1,1) #add a slice for height dim
        if fieldname == "wind_speed":
            u10 = super(ECMWF0125NCFile, self).read_values('u10m',slices=slices)
            v10 = super(ECMWF0125NCFile, self).read_values('v10m',slices=slices)
            val = numpy.sqrt(u10 * u10 + v10 * v10)
        elif fieldname == "wind_direction":
            u10 = super(ECMWF0125NCFile, self).read_values('u10m',slices=slices)
            v10 = super(ECMWF0125NCFile, self).read_values('v10m',slices=slices)
            val = uv2dir(u10, v10)
        elif fieldname in ['u10m','v10m']:
            val = super(ECMWF0125NCFile, self).read_values(fieldname,slices=slices)
        else:
            val = super(ECMWF0125NCFile, self).read_values(fieldname, slices)
        if len(val.shape) == 4: #remove de height dimension
            val = numpy.reshape(val,[1,val.shape[2],val.shape[3]])
        return val

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
#                                 ('z', self.get_dimsize('z')),
                                 ('y', self.get_dimsize('y')),
                                 ('x', self.get_dimsize('x'))
                                 ]),
                    datatype=numpy.dtype(numpy.float32),
                    units=VIRTUALFIELD_UNITS[fieldname]
                    )
            field.attach_storage(self.get_field_handler(fieldname))
        elif fieldname in ['u10m','v10m']:
            field = super(ECMWF0125NCFile, self).read_field(fieldname)
            vava=field.variable
            uni=field.units
            field = Field(
                vava,
                OrderedDict([('time', 1),
#                             ('z', self.get_dimsize('z')),
                             ('y', self.get_dimsize('y')),
                             ('x', self.get_dimsize('x'))
                             ]),
                datatype=numpy.dtype(numpy.float32),
                units=uni
                )
            field.attach_storage(self.get_field_handler(fieldname))
        else:
            field = super(ECMWF0125NCFile, self).read_field( fieldname)
        return field
    
    def get_dimsize(self, dimname):
        res = super(ECMWF0125NCFile, self).get_dimsize(dimname)
        return res
    
    def get_matching_dimname(self, geodimname):
        """Return the equivalent name in the native format for a standard
        dimension.

        This is a translation of the standard names to native ones. It is used
        for internal purpose only and should not be called directly.

        The standard dimension names are:

        * x, y, time for :class:`~cerbere.datamodel.grid.Grid`
        * row, cell, time for :class:`~cerbere.datamodel.swath.Swath` or
          :class:`~cerbere.datamodel.image.Image`

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        Args:
            dimname (str): standard dimension name.

        Returns:
            str: return the native name for the dimension. Return `dimname` if
                the input dimension has no standard name.

        See Also:
            see :func:`get_standard_dimname` for the reverse operation
        """

        dim_matching = {
                'time': 'time',
                'x':'longitude',
                'y':'latitude',
                'z':'height'
                }
        return dim_matching[geodimname]
    
    def get_geolocation_field(self, fieldname):
        if fieldname == 'z':
            res = 'height'
        else:
            res = super(ECMWF0125NCFile, self).get_geolocation_field(fieldname)
        return  res

    def get_standard_dimname(self, geodimname):
        '''here we keep the height dim to enable the formatting of slices'''
        dim_std = {
                'time': 'time',
                'longitude':'x',
                'latitude':'y',
                'height':'z',
                }
        if geodimname in dim_std:
            res = dim_std[geodimname]
        else:
            res = geodimname
        return res