# -*- coding: utf-8 -*-
"""
cerbere.mapper.aquariusl3podaacncfile
============================

Mapper class for AQUARIUS PODAAC L3 CAPv3.0 netcdf files 

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Antoine Grouazel <agrouaze@ifremer.fr>
.. codeauthor:: Antoine Grouazel <agrouaze@ifremer.fr>
"""
import datetime
import os
import calendar
import logging
from dateutil import parser

import numpy
import netCDF4
from collections import OrderedDict
from ..datamodel.field import Field
from ..datamodel.variable import Variable

from .. import READ_ONLY
from .ncfile import NCFile

import re
VIRTUALFIELD_DESCR = {
                      'time': '15th of the data acquisition month at 00:00:00'
                      }
VIRTUALFIELD_STDNAME = {
                      'time': 'time'
                      }
VIRTUALFIELD_UNITS = {
                'time': 'seconds since 1981-01-01 00:00:00'
                }
class AquariusL3PODAACNCFile(NCFile):
    def __init__(self, url=None, mode=READ_ONLY,
                 ncformat='NETCDF4', **kwargs):
        super(AquariusL3PODAACNCFile, self).__init__(url=url,
                                              mode=mode,
                                              ncformat=ncformat,
                                              **kwargs)
        
#     def get_dimensions(self, fieldname=None):
#         if fieldname != 'lat':
#             res = super(AquariusL3PODAACNCFile, self).get_dimensions(fieldname=fieldname)
#         else:
#             res = ('lat')
#         return res
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
        dims = self._handler.dimensions.keys()
        if 'idlon' in dims:
            dim_matching = {
                    'time': 'time',
                    'x': 'idlon',
                    'y': 'idlat',
                    }
        else:
            dim_matching = {
                    'time': 'time',
                    'x': 'lon',
                    'y': 'lat',
                    }
        return dim_matching[geodimname]

    def get_geolocation_field(self, fieldname):
        return fieldname
    
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

        # cell/row are reversed
        matching = {
            'time': 'time',
            'lon':'x',
            'lat':'y',
            'idlat':'y',
            'idlon':'x'
            }
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None

    def get_dimensions(self, fieldname=None):
        """
        """
        return super(AquariusL3PODAACNCFile, self).get_dimensions(fieldname)

    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature

        Excludes the geolocation fields in the returned list
        '''
        geophyvars = super(AquariusL3PODAACNCFile, self).get_fieldnames()
        geophyvars.append('time')
        geophyvars.remove('Equirectangular')
        return geophyvars

    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        if fieldname in ['time']:
            variable = Variable(
                    shortname=fieldname,
                    description=VIRTUALFIELD_DESCR[fieldname],
                    authority=self.get_naming_authority(),
                    standardname=VIRTUALFIELD_STDNAME[fieldname]
                    )
            field = Field(
                    variable,
                    OrderedDict([('time', 1)]),
                    datatype=numpy.dtype(numpy.float32),
                    units=VIRTUALFIELD_UNITS[fieldname]
                    )
            field.attach_storage(self.get_field_handler(fieldname))
        else:
            field = super(AquariusL3PODAACNCFile, self).read_field(fieldname)
        return field


    def read_values(self, fieldname, slices=None):
        """Read the data of a field.

        Args:
            fieldname (str): name of the field which to read the data from

            slices (list of slice, optional): list of slices for the field if
                subsetting is requested. A slice must then be provided for each
                field dimension. The slices are relative to the opened view
                (see :func:open) if a view was set when opening the file.

        Return:
            MaskedArray: array of data read. Array type is the same as the
                storage type.
        """
        if fieldname == 'time':
            tmp = os.path.basename(self._url)[3:9]
            datefile = datetime.datetime.strptime(tmp,'%Y%m')
            values = netCDF4.date2num(datefile,VIRTUALFIELD_UNITS[fieldname])
        else:
            values = super(AquariusL3PODAACNCFile, self).read_values(fieldname,
                                                               slices)
            values.mask = False
        return values
    def get_start_time(self):
        yearmonth_file = self.read_values('time')
        yearmonth_dt_file = netCDF4.num2date(yearmonth_file,VIRTUALFIELD_UNITS['time'])
        year = yearmonth_dt_file.year
        month = yearmonth_dt_file.month
        st = datetime.datetime(year,month,1)
        return st
    def get_end_time(self):
        yearmonth_file = self.read_values('time')
        yearmonth_dt_file = netCDF4.num2date(yearmonth_file,VIRTUALFIELD_UNITS['time'])
        year = yearmonth_dt_file.year
        month = yearmonth_dt_file.month
        b,e = calendar.monthrange(year,month)
        st = datetime.datetime(year,month,e)
        return st
    def get_cycle_number(self):
        return None


