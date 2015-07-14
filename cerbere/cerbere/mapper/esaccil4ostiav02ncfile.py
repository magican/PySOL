# -*- coding: utf-8 -*-
"""
cerbere.mapper.esaccil4ostiav02ncfile
============================

Mapper class for ESA CCI L4 OSTIA v02.0 netcdf files 

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Antoine Grouazel <agrouaze@ifremer.fr>
.. codeauthor:: Antoine Grouazel <agrouaze@ifremer.fr>
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

VIRTUALFIELD_DESCR = {
                      'sea_surface_temperature': 'analysed sea surface temperature'
                      }
VIRTUALFIELD_STDNAME = {
                      'sea_surface_temperature': 'sea_surface_temperature'
                      }
VIRTUALFIELD_UNITS = {
                'sea_surface_temperature': 'Kelvin'
                }
class ESACCIL4OstiaV02NCFile(NCFile):

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
                'bnds':'bnds',
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
            'lon':'x',
            'lat':'y',
            'bnds':'bnds',
            }
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None

    def get_dimensions(self, fieldname=None):
        """
        """
        return super(ESACCIL4OstiaV02NCFile, self).get_dimensions(fieldname)

    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature

        Excludes the geolocation fields in the returned list
        '''
        geophyvars = super(ESACCIL4OstiaV02NCFile, self).get_fieldnames()
        geophyvars.remove('analysed_sst')
        geophyvars.append('sea_surface_temperature')
        return geophyvars

    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        if fieldname in ['sea_surface_temperature']:
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
            field = super(ESACCIL4OstiaV02NCFile, self).read_field(fieldname)
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
        if fieldname == 'sea_surface_temperature':
            field_request = 'analysed_sst'
        else:
            field_request = fieldname
        values = super(ESACCIL4OstiaV02NCFile, self).read_values(field_request,
                                                               slices)
        return values

    def get_cycle_number(self):
        return None


