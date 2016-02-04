'''
cerbere.mapper.rapidscatpodaacncfile
'''
import netCDF4
import numpy
import ast
import os
import logging
import pdb
from datetime import datetime
from collections import OrderedDict
#from .ncfile import NCFile
from .. import READ_ONLY, WRITE_NEW
from cerbere.mapper.ncfile import NCFile
from cerbere.datamodel.grid import Grid
from ..datamodel.field import Field
# from ..science.wind import uv2dir
from cerform.wind import uv2dir
# from wind import uv2dir
from ..datamodel.variable import Variable

# VIRTUALFIELD_DESCR = {'time': 'time',
#                       'wind_speed': 'wind speed',
#                       'wind_direction': 'wind direction (where the wind is blowing)'
#                       }
# VIRTUALFIELD_STDNAME = {'time': 'time',
#                         'wind_speed': 'wind_speed',
#                       'wind_direction': 'wind_to_direction'
#                       }
VIRTUALFIELD_UNITS = {
#                     'time':'seconds since 1980-01-01 00:00:00',
                'wind_speed': 'm s-1',
                'wind_direction': 'degrees'
                }
class RapidScatPODAACNCFile(NCFile):
    
    def __init__(self, url=None, mode=READ_ONLY,ncformat='NETCDF3_CLASSIC', **kwargs):
        super(RapidScatPODAACNCFile,self).__init__(url=url,mode=mode,ncformat=ncformat,depths=True,center_on_greenwhich=False,**kwargs)
        if self._handler is None:
                self._handler = self.get_handler()
        return  

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
                    'row':'along_track',
                    'cell':'cross_track',
                    }
        if geodimname in dim_matching:
            res= dim_matching[geodimname]
        else:
            res = geodimname
        return res

    def get_standard_dimname(self, geodimname):
        """
        Returns the equivalent standard dimension name for a
        dimension in the native format.

        This is a translation of the native names to standard ones. It is used
        for internal purpose and should not be called directly.

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        Args:
            dimname (string): native dimension name

        Return:
            str: the (translated) standard name for the dimension. Return
            `dimname` if the input dimension has no standard name.

        See Also:
            see :func:`get_matching_dimname` for the reverse operation
        """
        dim_matching = {
                'time': 'time',
                'along_track':'row',
                'cross_track':'cell',
                }
        if geodimname in dim_matching:
            res= dim_matching[geodimname]
        else:
            res = geodimname

        return res
        
    def read_field(self, fieldname):
        """
        Return the field, without its values.
  
        Actual values can be retrieved with read_values() method.
        """
#         if fieldname in ['wind_speed', 'wind_direction']:
#             # create a virtual field
#             variable = Variable(
#                     shortname=fieldname,
#                     description=VIRTUALFIELD_DESCR[fieldname],
#                     authority=self.get_naming_authority(),
#                     standardname=VIRTUALFIELD_STDNAME[fieldname]
#                     )
#             field = Field(
#                     variable,
#                     OrderedDict([
#                                  #('z', 1),
#                                  ('row', self.get_dimsize('row')),
#                                  ('cell', self.get_dimsize('cell'))
#                                  ]),
#                     datatype=numpy.dtype(numpy.float32),
#                     units=VIRTUALFIELD_UNITS[fieldname]
#                     )
#             field.attach_storage(self.get_field_handler(fieldname))
        
#             variable = Variable(
#                     shortname=fieldname,
#                     description=fieldname,
#                     authority=self.get_naming_authority(),
#                     standardname=fieldname
#                     )
#             field = Field(
#                     variable,
#                     OrderedDict([
#                                 ('row', self.get_dimsize('row')),
#                                 ('cell', self.get_dimsize('cell'))
#                                 ]),
# #                     OrderedDict([('time',self.get_dimsize('row'))]),
#                     datatype=numpy.dtype(numpy.float32),
#                     units=VIRTUALFIELD_UNITS[fieldname]
#                     )
#             field.attach_storage(self.get_field_handler(fieldname))
            

        field = super(RapidScatPODAACNCFile,self).read_field(fieldname)
        if fieldname == 'time':
            dims = OrderedDict([
                                ('row', self.get_dimsize('row')),
                                ('cell', self.get_dimsize('cell'))
                                ])
            field.dimensions = dims
        return field

    def read_values(self, fieldname, slices=None):
#         if fieldname == 'time':
#             val = super(RapidScatPODAACNCFile, self).read_values('time')
#             res = numpy.repeat(val,self.get_dimsize('cell'))
#             res = numpy.array([val,]*self.get_dimsize('cell')).transpose()
#         else:
        res = super(RapidScatPODAACNCFile,self).read_values(fieldname)
        return res

        
    def get_start_time(self):
        '''minimal time of the file'''
        date = self.read_global_attribute('RangeBeginningDate')
        times = self.read_global_attribute('RangeBeginningTime')
        resf = datetime.strptime(date+times,'%Y-%j%H:%M:%S.%f')
        return resf
#     
    def get_end_time(self):
        '''maximal time of the file'''
        date = self.read_global_attribute('RangeEndingDate')
        times = self.read_global_attribute('RangeEndingTime')
        resf = datetime.strptime(date+times,'%Y-%j%H:%M:%S.%f')
        return resf
    
    def get_collection_id(self):
        """return the identifier of the product collection"""
        return 'PODAAC-ISS-RAPIDSCAT-L2B-0125'

    