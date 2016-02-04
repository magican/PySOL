# -*- coding: utf-8 -*-
"""
cerbere.mapper.smosncfile
============================

Mapper class for SMOS netcdf files L3 

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Antoine Grouazel <agrouaze@ifremer.fr>
.. codeauthor:: JAntoine Grouazel <agrouaze@ifremer.fr>
"""
import datetime
import os
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

VIRTUALFIELD_DESCR = {'sea_surface_salinity': 'Practical Sea Surface Salinity',
                      'sea_surface_temperature': 'ECMWF sea surface temperature'
                      }
VIRTUALFIELD_STDNAME = {'sea_surface_salinity': 'sea_surface_salinity',
                      'sea_surface_temperature': 'sea_surface_temperature'
                      }
VIRTUALFIELD_UNITS = {
                'sea_surface_salinity': 'P.S.S.',
                'sea_surface_temperature': 'Kelvins'
                }
class SMOSNCFile(NCFile):

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
#                     'cell': 'nj',
#                     'row': 'ni',
                'x': 'longitude',
                'y': 'latitude',
                'date_start':'date_start',
                'date_stop':'date_stop',
                }
        return dim_matching[geodimname]


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
#                 'ni': 'row',
#                 'nj': 'cell',
            'longitude':'x',
            'latitude':'y',
            'date_start':'date_start',
            'date_stop':'date_stop',
            }
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None

    def get_dimensions(self, fieldname=None):
        """
        """
        return super(SMOSNCFile, self).get_dimensions(fieldname)

    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature

        Excludes the geolocation fields in the returned list
        '''
        geophyvars = super(SMOSNCFile, self).get_fieldnames()
        geophyvars.remove('sss')
        geophyvars.remove('sst')
        geophyvars.append('sea_surface_salinity')
        geophyvars.append('sea_surface_temperature')
        return geophyvars

    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        if fieldname in ['sea_surface_salinity','sea_surface_temperature']:
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
        else:
            field = super(SMOSNCFile, self).read_field(fieldname)
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
        if fieldname == 'sea_surface_salinity':
            field_request = 'sss'
        elif fieldname == 'sea_surface_temperature':
            field_request = 'sst'
        else:
            field_request = fieldname
        values = super(SMOSNCFile, self).read_values(field_request,
                                                               slices)
        return values

    def get_cycle_number(self):
        return None


    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage"""
        tmp = self.read_field('date_start')
        print dir(tmp)
        res = netCDF4.num2date(tmp.get_values(), tmp.units, calendar='standard')
        return res

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage"""
        tmp = self.read_field('date_stop')
        res = netCDF4.num2date(tmp.get_values(), tmp.units, calendar='standard')
        return res
