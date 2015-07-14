# -*- coding: utf-8 -*-
"""
cerbere.mapper.clsncfile
========================

Mapper class for netcdf files by CLS/AVISO

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import copy
import collections
from dateutil import parser

import netCDF4
import numpy
from datetime import datetime, timedelta

from .ncfile import NCFile
from cerbere import DEFAULT_TIME_UNITS
from ..datamodel.field import Field
from ..datamodel.variable import Variable


class CLSNCFile(NCFile):
    """Mapper class for the netcdf files from CLS/AVISO.

    The CLS files are in netCDF format but don't have a time variable
    (time is in a global attribute). They are not really CF compliant
    neither.

    This class is inherited from the :class:`NCFile` mapper but extended
    to mimic the presence of a `time` variable.
    """

    def get_geolocation_field(self, fieldname):
        matching = {'time': 'time',
                    'lon': 'NbLongitudes',
                    'lat': 'NbLatitudes'}
        if fieldname in matching:
            return matching[fieldname]
        else:
            return None

    def get_matching_dimname(self, geodimname):
        """
        Return the equivalent name in the native format for a standard
        dimension.

        Args
        ====
        dimname (string)
            standard dimension name

        Returns
        =======
        string,the native name for the dimension. Return `dimname` if the
        input dimension has no standard name.

        See
        ===
        see :func:`get_standard_dimname` for the reverse operation
        """
        matching = {'time': 'time',
                    'x': 'NbLongitudes',
                    'y': 'NbLatitudes'}
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None

    def get_standard_dimname(self, geodimname):
        """
        Returns the equivalent standard dimension name for a
        dimension in the native format.

        Args
        ====
        dimname (string)
            native dimension name

        Return
        ======
        string, the (translated) standard name for the dimension. Return
        `dimname` if the input dimension has no equivalent standard name.

        See
        ===
        see :func:`get_matching_dimname` for the reverse operation
        """
        matching = {'time': 'time',
                    'NbLongitudes': 'x',
                    'NbLatitudes': 'y'}
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None

    def get_dimsize(self, dimname):
        if dimname == 'time':
            return 1
        else:
            super(CLSNCFile, self).get_dimsize(dimname)

    def read_times(self, slices=None):
        """Read time values of a file"""
        times = netCDF4.num2date(
            parser.parse(self.get_handler().variables['Grid_0001'].date),
            DEFAULT_TIME_UNITS
            )
        return numpy.ma.array([times])

    def read_field_attributes(self, fieldname):
        """
        """
        if fieldname == 'time':
            return {}
        else:
            return super(CLSNCFile, self).read_field_attributes(fieldname)

    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        # special implementation case for time field which is not
        # available as a variable in CLS files
        if fieldname != 'time':
            return NCFile.read_field(self, fieldname)
        else:
            # create a field for time
            variable = Variable(
                    shortname=fieldname,
                    description='time',
                    authority=self.get_naming_authority(),
                    standardname='time'
                    )
            field = Field(
                    variable,
                    collections.OrderedDict([('time', 1)]),
                    datatype=numpy.dtype(numpy.int64),
                    units=DEFAULT_TIME_UNITS
                    )
            field.attach_storage(self.get_field_handler(fieldname))
        return field

    def get_fieldnames(self):
        """
        Returns the list of geophysical fields stored for the feature

        Excludes the geolocation fields in the returned list
        """
        fieldnames = super(CLSNCFile, self).get_fieldnames()
        fieldnames.remove('LatLon')
        fieldnames.remove('LatLonMin')
        fieldnames.remove('LatLonStep')
        return fieldnames

    def read_values(self, fieldname, slices=None):
        """
        """
        if fieldname == 'time':
            if self._handler is None:
                self._handler = self.get_handler()
            times = self.read_times()
            return times
        else:
            data = super(CLSNCFile, self).read_values(fieldname,
                                                      slices=slices)
            if fieldname == 'lon':
                return numpy.ma.concatenate([data[data.size/2:],
                                            data[0:data.size/2]]
                                           )
            elif 'Grid_' in fieldname:
                data = data.transpose()
                print "SHAPE", data.shape, fieldname
                xsize = data.shape[1]
                data = numpy.ma.concatenate([data[:, xsize/2:],
                                            data[:, 0:xsize/2]], axis=1
                                           )
                print "SHAPE", data.shape
            return data

    def get_start_time(self):
        """Return start of temporal coverage"""
        return parser.parse(self.get_handler().variables['Grid_0001'].date)

    def get_end_time(self):
        """Return end of temporal coverage"""
        return parser.parse(self.get_handler().variables['Grid_0001'].date)
