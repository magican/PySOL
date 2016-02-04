# -*- coding: utf-8 -*-
"""
.. mod:: cerbere.mapper.mwferslandmaskncfile

Mapper class for MWF-ERS Land mask in netCDF

:copyright: Copyright 2014 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: agrouaze@ifremer.fr
"""

import collections

import netCDF4
import numpy
from datetime import datetime
import logging
from .ncfile import NCFile
from ..datamodel.field import Field
from ..datamodel.variable import Variable


class MWFERSLandmaskNCFile(NCFile):
    """Mapper class for the MWF-ERS Land mask in netCDF

    The MWF-ERS mask is a vectorized 2-dimensional array.
    in the netcdf file
    sea is 0
    and land is 2
    in the mapper it is
    sea : -1
    land: 0
    """
    def get_geolocation_field(self, fieldname):
#         logging.debug('get_geolocation_field mwfers fieldname %s',fieldname)
        if fieldname == 'lon':
            res = 'longitude'
        elif fieldname == 'lat':
            res = 'latitude'
        else:
            res= fieldname
        return  res

    def get_matching_dimname(self, geodimname):
        if geodimname == 'y':
            res = 'latitude'
        elif geodimname == 'x':
            res = 'longitude'
        else:
            res = geodimname
        return res

    def get_standard_dimname(self, geodimname):
        return geodimname

    def get_fieldnames(self):
        return ['lon','lat','landmask','time']

    def read_field_attributes(self, fieldname):
        """
        """
        if fieldname == 'landmask':
            res = 'quality_flag'
        elif fieldname == 'lon':
            res = 'longitude'
        elif fieldname == 'lat':
            res = 'latitude'
        else:
            res = fieldname
        return super(MWFERSLandmaskNCFile, self).read_field_attributes(res)

    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        # special implementation case for time field which is not
        # available as a variable a land mask
        logging.debug('fieldname %s',fieldname)
        if fieldname == 'lon':
            variable = Variable(shortname=fieldname,
                                description='lon',
                                authority=self.get_naming_authority(),
                                standardname='longitude'
                                )
            dimensions = self.get_full_dimensions('longitude')
            field = Field(variable,
                          dimensions=dimensions,
                          datatype=numpy.dtype(numpy.float),
                          units='degree'
                          )
            field.attach_storage(self.get_field_handler(fieldname))
            field.dimensions = dimensions
        elif fieldname == 'lat':
            variable = Variable(shortname=fieldname,
                                description='lat',
                                authority=self.get_naming_authority(),
                                standardname='latitude'
                                )
            dimensions = self.get_full_dimensions('latitude')
            field = Field(variable,
                          dimensions=dimensions,
                          datatype=numpy.dtype(numpy.float),
                          units='degree'
                          )
            field.attach_storage(self.get_field_handler(fieldname))
            field.dimensions = dimensions
        elif fieldname == 'landmask':
            # create a field for time
            variable = Variable(shortname=fieldname,
                                description='landmask',
                                authority=self.get_naming_authority(),
                                standardname='landmask'
                                )
            dimensions_landmask = self.get_full_dimensions('quality_flag')
            logging.debug('dimensions_landmask %s %s',dimensions_landmask,type(dimensions_landmask))
            field = Field(variable,
                          dimensions=dimensions_landmask,
                          datatype=numpy.dtype(numpy.bool),
                          units=''
                          )
            field.attach_storage(self.get_field_handler(fieldname))
            field.dimensions = dimensions_landmask
            logging.debug('verif dimensions.values() %s',field.dimensions.values())
        else:
            return super(MWFERSLandmaskNCFile, self).read_field(fieldname)
        return field

    def read_values(self, fieldname, slices=None):
        """
        """
        if fieldname == 'landmask':
            native = 'quality_flag'
        elif fieldname == 'lon':
            native = 'longitude'
        elif fieldname == 'lat':
            native = 'latitude'
        else:
            native = fieldname
        res = super(MWFERSLandmaskNCFile, self).read_values(
                native,
                slices=slices
                )
#         logging.debug('ers mask slice %s',slices)
        if fieldname == 'landmask':
            res[res==0] = -1 #sea
            res[res==2] = 0 #land
        return res

#     def get_bbox(self):
#         """Return the bounding box of the feature, as a tuple
#         (lonmin, latmin, lonmax, latmax)
#         """
#         lon =self.read_values('lon')
#         lat = self.read_values('lat')
#         lonmax = numpy.amax(lon)
#         lonmin = numpy.amin(lon)
#         latmax = numpy.amax(lat)
#         latmin = numpy.amin(lat)
#         return(lonmin,latmin,lonmax,latmax)
#         return (self.get_handler().west_longitude,
#                 self.get_handler().south_latitude,
#                 self.get_handler().east_longitude,
#                 self.get_handler().north_latitude
#                 )

