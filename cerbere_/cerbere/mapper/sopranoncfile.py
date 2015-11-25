# -*- coding: utf-8 -*-
"""
cerbere.mapper.sopranoncfile
============================

Mapper class for SOPRANO netcdf files by CLS

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import copy
import collections

import netCDF4
import numpy
from datetime import datetime, timedelta

from .ncfile import NCFile
from ..datamodel.field import Field
from ..datamodel.variable import Variable


class SopranoNCFile(NCFile):
    """Mapper class for the Soprano files from CLS.

    The Soprano files are in netCDF format but don't have a time variable
    (time is in a global attribute). This class is inherited from the
    :class:`NCFile` mapper but extended to mimic the presence of a `time`
    variable.
    """

    def get_geolocation_field(self, fieldname):
        matching = {'time': 'time',
                    'lon': 'longitude',
                    'lat': 'latitude'}
        if fieldname in matching:
            return matching[fieldname]
        else:
            return None

    def get_matching_dimname(self, geodimname):
        matching = {'time': 'time',
                    'row': 'az_size',
                    'cell': 'ra_size'}
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None

    def get_standard_dimname(self, geodimname):
        matching = {'time': 'time',
                    'az_size': 'row',
                    'ra_size': 'cell'}
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None

    def get_dimsize(self, dimname):
        if dimname == 'time':
            return 1
        else:
            super(SopranoNCFile, self).get_dimsize(dimname)

    def read_times(self, slices=None):
        """Read time values of a file"""
        times = netCDF4.num2date(
            datetime.strptime(
                    self.get_handler().SOURCE_START_DATE.split('.')[0],
                    '%Y%m%d%H%M%S'
                    )
            )
        return numpy.ma.array([times])

    def read_field_attributes(self, fieldname):
        """
        """
        native_fieldname = self.get_geolocation_field(fieldname)
        if native_fieldname is None:
            native_fieldname = fieldname
        if native_fieldname == 'time':
            return {}
        else:
            attrs = self.get_handler().variables[native_fieldname].__dict__
            return copy.copy(attrs)

    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        # special implementation case for time field which is not
        # available as a variable in Soprano files
        if fieldname != 'time':
            return NCFile.read_field(self, fieldname)
        else:
            # create a field for time
            variable = Variable(
                    shortname=fieldname,
                    description='acquisition time of image',
                    authority=self.get_naming_authority(),
                    standardname='time'
                    )
            field = Field(
                    variable,
                    collections.OrderedDict([('time', 1)]),
                    datatype=numpy.dtype(numpy.int64),
                    units='seconds since 1981-01-01 00:00:00'
                    )
            field.attach_storage(self.get_field_handler(fieldname))
        return field

    def read_values(self, fieldname, slices=None):
        """
        """
        if fieldname == 'time':
            if self._handler is None:
                self._handler = self.get_handler()
            timeval = self.get_start_time()
            return numpy.ma.array([netCDF4.date2num(
                                    timeval,
                                    'seconds since 1981-01-01 00:00:00')
                                  ])
        else:
            return super(SopranoNCFile, self).read_values(fieldname,
                                                          slices=slices)

    def get_start_time(self):
        """Return start of temporal coverage"""
        start = datetime.strptime(
            self.get_handler().SOURCE_START_DATE.split('.')[0],
            '%Y%m%d%H%M%S'
            )
        return start

    def get_end_time(self):
        """Return end of temporal coverage"""
        start = datetime.strptime(
            self.get_handler().SOURCE_START_DATE.split('.')[0],
            '%Y%m%d%H%M%S'
            )
        delta = timedelta(
                    seconds=float(self.get_handler().SOURCE_ACQ_DURATION)
                    )
        return start + delta

    def get_bbox(self):
        """
        return the bounding box of the feature, as a tuple
        (lonmin, latmin, lonmax, latmax)
        """
        return (self.get_handler().PROCESSING_LIMIT_WEST,
                self.get_handler().PROCESSING_LIMIT_SOUTH,
                self.get_handler().PROCESSING_LIMIT_EAST,
                self.get_handler().PROCESSING_LIMIT_NORTH
                )