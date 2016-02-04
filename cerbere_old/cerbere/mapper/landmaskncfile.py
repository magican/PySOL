# -*- coding: utf-8 -*-
"""
cerbere.mapper.landmaskncfile
============================

Mapper class for landmask 10deg and 2deg netcdf files 

:copyright: Copyright 2015 Ifremer / Cersat.
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
                      'time': 'dummy date',
                      'landmask':'sea/land mask',
                      }
VIRTUALFIELD_STDNAME = {
                      'time': 'time',
                      'landmask':'mask',
                      }
VIRTUALFIELD_UNITS = {
                      
                'time': 'days since 1981-01-01 00:00:00'
                }
class LandMaskNCFile(NCFile):
    def __init__(self, url=None, mode=READ_ONLY,
                 ncformat='NETCDF4', **kwargs):
        super(LandMaskNCFile, self).__init__(url=url,
                                              mode=mode,
                                              ncformat=ncformat,
                                              is_reference=True,
                                              **kwargs)
        return
    
    def get_dimsize(self, dimname):
        """returns the size of a dimension

        If a view was set when opening the file, the size of the view subset
        is returned.

        Args:
            dimname (string): name of the dimension for which the size is
                inquired. Can be provided as the native or standard dimension
                name.

        Return:
            int: size of the dimension. 
        """
        if dimname == 'time':
            res = 1
        else:
            res = super(LandMaskNCFile, self).get_dimsize(dimname)
        return res
    
    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature
 
        Excludes the geolocation fields in the returned list
        '''
        geophyvars = super(LandMaskNCFile, self).get_fieldnames()
        geophyvars.append('landmask')
        geophyvars.append('lon')
        geophyvars.append('lat')
        geophyvars.append('time')
        return geophyvars
    
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
                    'x': 'lon',
                    'y': 'lat',
                    }
        return dim_matching[geodimname]
    
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
        elif fieldname =='landmask':
            tmpfield = super(LandMaskNCFile, self).read_field(VIRTUALFIELD_STDNAME[fieldname])
            dims = tmpfield.dimensions
            variable = Variable(
                    shortname=fieldname,
                    description=VIRTUALFIELD_DESCR[fieldname],
                    authority=self.get_naming_authority(),
                    standardname=VIRTUALFIELD_STDNAME[fieldname]
                    )
            field = Field(
                    variable,
                    dims,
                    datatype=numpy.dtype(numpy.float32),
                    )
            field.attach_storage(self.get_field_handler(fieldname))
        else:
            field = super(LandMaskNCFile, self).read_field(fieldname)
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
            dummy = datetime.datetime(2000,1,1)
#             dummy = datetime.datetime.max
            values = netCDF4.date2num(dummy,VIRTUALFIELD_UNITS[fieldname])
        elif fieldname =='landmask':
            values = super(LandMaskNCFile, self).read_values(VIRTUALFIELD_STDNAME[fieldname],
                                                               slices)
            #apply convention sea=-1 and land=0
            values[values==1] = -1
            values[values==2] = 0
        else:
            values = super(LandMaskNCFile, self).read_values(fieldname,
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

