# -*- coding: utf-8 -*-
"""
.. mod:: cerbere.mapper.navolandmaskncfile

Mapper class for NAVO Land mask in netCDF

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import collections

import netCDF4
import numpy
from datetime import datetime
import logging
from .ncfile import NCFile
from ..datamodel.field import Field
from ..datamodel.variable import Variable


class NAVOLandmaskNCFile(NCFile):
    """Mapper class for the NAVO Land mask in netCDF

    The NAVO mask is a vectorized 2-dimensional array.
    in netcdf file sea is -1
    and land is 0
    """

    def get_geolocation_field(self, fieldname):
        return  fieldname

#     def get_matching_dimname(self, geodimname):
#         return geodimname
    def get_matching_dimname(self, geodimname):
        if geodimname == 'y':
            res = 'lat'
        elif geodimname == 'x':
            res = 'lon'
        else:
            res = geodimname
        return res

    def get_standard_dimname(self, geodimname):
        return geodimname

    def get_dimsize(self, dimname):
        if dimname == 'time':
            return 1
        else:
            super(NAVOLandmaskNCFile, self).get_dimsize(dimname)
    def get_fieldnames(self):
        return ['lon','lat','landmask','time']
    def get_dimensions(self,fieldname):
        if fieldname == 'landmask':
            native = 'dst'
        else:
            native = fieldname
        return super(NAVOLandmaskNCFile, self).get_dimensions(native
                                                              )
    def read_times(self, slices=None):
        """Read time values of a file"""
        times = netCDF4.num2date(datetime.datetime(1, 1, 1, 0, 0, 0))
        return numpy.ma.array([times])

    def read_field_attributes(self, fieldname):
        """
        """
        if fieldname == 'time':
            return {}
        elif fieldname == 'landmask':
            return super(NAVOLandmaskNCFile, self).read_field_attributes('dst')
        else:
            return super(NAVOLandmaskNCFile, self).read_field_attributes(fieldname)

    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        # special implementation case for time field which is not
        # available as a variable a land mask
        if fieldname == 'time':
            # create a field for time
            variable = Variable(shortname=fieldname,
                                description='time',
                                authority=self.get_naming_authority(),
                                standardname='time'
                                )
            field = Field(variable,
                          dimensions=collections.OrderedDict([('time', 1)]),
                          datatype=numpy.dtype(numpy.int64),
                          units='seconds since 0001-01-01 00:00:00'
                          )
            field.attach_storage(self.get_field_handler(fieldname))
            field.dimensions=collections.OrderedDict([('time', 1)])
            
        elif fieldname == 'landmask':
            variable = Variable(shortname=fieldname,
                            description='landmask',
                            authority=self.get_naming_authority(),
                            standardname='landmask'
                            )
            dimensions_landmask = self.get_full_dimensions('dst')
            field = Field(variable,
                      dimensions=dimensions_landmask,
                      datatype=numpy.dtype(numpy.bool),
                      units=''
                      )
            field.attach_storage(self.get_field_handler(fieldname))
            field.dimensions = dimensions_landmask
        else:
            return super(NAVOLandmaskNCFile, self).read_field(fieldname)
        return field

    def read_values(self, fieldname, slices=None):
        """
        """
        if fieldname == 'time':
            return numpy.ma.array([0])
        elif fieldname == 'landmask':
            tmp = super(NAVOLandmaskNCFile, self).read_values('dst',slices=slices)
            return tmp
        else:
            return super(NAVOLandmaskNCFile, self).read_values(
                fieldname,
                slices=slices
                )

    def get_start_time(self):
        """Return start of temporal coverage"""
        return datetime(1, 1, 1, 0, 0, 0)

    def get_end_time(self):
        """Return end of temporal coverage"""
        return datetime(1, 1, 1, 0, 0, 0)

    def get_bbox(self):
        """Return the bounding box of the feature, as a tuple
        (lonmin, latmin, lonmax, latmax)
        """
        lon =self.read_values('lon')
        lat = self.read_values('lat')
        lonmax = numpy.amax(lon)
        lonmin = numpy.amin(lon)
        latmax = numpy.amax(lat)
        latmin = numpy.amin(lat)
        return(lonmin,latmin,lonmax,latmax)
#         return (self.get_handler().west_longitude,
#                 self.get_handler().south_latitude,
#                 self.get_handler().east_longitude,
#                 self.get_handler().north_latitude
#                 )
