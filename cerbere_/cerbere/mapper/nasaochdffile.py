# encoding: utf-8
"""
cerbere.mapper.nasaochdffile
============================

Mapper class for the ocean color HDF files from NASA

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
import datetime
from collections import OrderedDict

from numpy import dtype, float32, resize
from netCDF4 import num2date

from .. import READ_ONLY, DEFAULT_TIME_UNITS
from .hdffile import HDFFile
from ..datamodel.field import Field
from ..datamodel.variable import Variable


class NASAOCHDFFile(HDFFile):
    """"mapper class for ocean color HDF files from NASA"""
    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        """
        """
        HDFFile.__init__(self, url=url, mode=mode, **kwargs)
        return

    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature
        '''
        fields = self.get_handler().datasets().keys()
        # remove here time/space information to keep only geophysical fields
        for field in ['year', 'day', 'msec', 'latitude', 'longitude']:
            if field in fields:
                fields.remove(field)
        return fields

    def get_dimensions(self, fieldname=None):
        """
        Return the standard dimension names of a file or a field in the file

        :keyword fieldname: the field from which to get the dimension names.
            For a geolocation field, use the cerbere standard name 
            (time, lat, lon), though native field name will work too.
        :type fieldname: str

        :return: the standard dimensions of the field or file.
        :rtype: tuple of strings
        """
        dims = list(super(NASAOCHDFFile, self).get_dimensions(fieldname))
        # correct a bug of pyhdf when the same dimension is used twice in a
        # field dimension
        if fieldname == 'sen_mat':
            dims[dims.index('vec')] = 'vec_1'
            dims.append('vec_2')
        return tuple(dims)

    def get_geolocation_field(self, fieldname):
        # time here is a virtual variable (composed from msec and year)
        matching = {'lat': 'latitude',
                    'lon': 'longitude',
                    'time': 'time'
                    }
        if fieldname in matching:
            return matching[fieldname]
        else:
            return None

    def get_matching_dimname(self, dimname):
        """
        Return the equivalent name in the native format for a standard
        dimension.

        This is a translation of the standard names to native ones. It is used
        for internal purpose only and should not be called directly.

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        :arg dimname: standard dimension name
        :type dimname: str

        :rtype: str
        :return: return the native name for the dimension.
            return `dimname` if the input dimension has no standard name.

        *See*

        see :func:`get_standard_dimname` for the reverse operation
        """
        matching = {
                'row': 'Number of Scan Lines',
                'cell': 'Number of Pixel Control Points',
                'time': 'time',
                'vec_1': 'vec',
                'vec_2': 'vec'
                }
        if dimname in matching:
            return matching[dimname]
        else:
            return dimname

    def get_standard_dimname(self, dimname):
        """
        Returns the equivalent standard dimension name for a
        dimension in the native format.

        This is a translation of the native names to standard ones. It is used
        for internal purpose and should not be called directly.

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        :arg dimname: native dimension name
        :type dimname: str

        :rtype: str
        :return: return the (translated) standard name for the dimension.
            return `dimname` if the input dimension has no standard name.

        *See*

        see :func:`get_matching_dimname` for the reverse operation
        """
        matching = {
                'Number of Scan Lines': 'row',
                'Pixels per Scan Line': 'cell',
                'Number of Pixel Control Points': 'cell',
                'time': 'time'
                }
        if dimname in matching:
            return matching[dimname]
        else:
            return dimname

    def read_field(self, fieldname):
        """
        """
        native_fieldname = self.get_geolocation_field(fieldname)
        if native_fieldname is None:
            native_fieldname = fieldname
        if native_fieldname == 'time':
            # time is a virtual variable comining year and msec, and giving
            # time value at each swath pixel, according to swath data model
            dims = OrderedDict([
                ('row', self.get_dimsize('Number of Scan Lines')),
                ('cell', self.get_dimsize('Pixels per Scan Line')),
                ])
            var = Variable(
                            shortname='time',
                            description='time',
                            authority=self.get_naming_authority(),
                            standardname='time'
                            )
            field = Field(
                            var,
                            OrderedDict(dims),
                            datatype=dtype(float32),
                            units=DEFAULT_TIME_UNITS
                            )
            field.attach_storage(self.get_field_handler(fieldname))
        else:
            field = super(NASAOCHDFFile, self).read_field(fieldname)
            attributes = self.read_field_attributes(fieldname)
            if 'standard_name' in attributes:
                field.variable.standardname = attributes['standard_name']
            if 'long_name' in attributes:
                field.variable.description = attributes['long_name']
            if 'units' in attributes:
                field.units = attributes['units']
        return field

    def read_values(self, fieldname, slices=None):
        """
        """
        if fieldname == 'time':
            # time is a virtual variable comining year and msec, and giving
            # time value at each swath pixel, according to swath data model
            nbrows, nbcells \
                = self.get_handler().select('latitude').dimensions().values()
            year = HDFFile.read_values(self, 'year')
            msec = HDFFile.read_values(self, 'msec')
            day = HDFFile.read_values(self, 'day')
            yearmin = year.min()
            yearmax = year.max()
            if yearmin == yearmax - 1:
                diffyearmin = (datetime.datetime(yearmin, 1, 1)\
                       - datetime.datetime(1981, 1, 1)).total_seconds()
                diffyearmax = (datetime.datetime(yearmax, 1, 1)\
                       - datetime.datetime(1981, 1, 1)).total_seconds()
                year = (year - yearmin) * diffyearmax + (yearmax - year) * diffyearmin
            elif yearmin == yearmax:
                diffyearmin = (datetime.datetime(yearmin, 1, 1)\
                       - datetime.datetime(1981, 1, 1)).total_seconds()
                year[:] = diffyearmin
            else:
                raise Exception("Data span three years?")
            values = year + day * 86400. + msec / 1000.
            values = resize(values, (nbcells, nbrows,)).transpose()
            if slices is None:
                return values
            else:
                return values[slices]
        else:
            values = super(NASAOCHDFFile, self).read_values(fieldname, slices)
            return values

    def get_start_time(self):
        year = self.read_global_attribute('Start Year')
        day = self.read_global_attribute('Start Day')
        msec = self.read_global_attribute('Start Millisec')
        date = (datetime.datetime(year, 1, 1)\
                       - datetime.datetime(1981, 1, 1)).total_seconds()\
                    + day * 86400. + msec / 1000.
        return num2date(date, DEFAULT_TIME_UNITS)

    def get_end_time(self):
        year = self.read_global_attribute('End Year')
        day = self.read_global_attribute('End Day')
        msec = self.read_global_attribute('End Millisec')
        date = (datetime.datetime(year, 1, 1)\
                       - datetime.datetime(1981, 1, 1)).total_seconds()\
                    + day * 86400. + msec / 1000.
        return num2date(date, DEFAULT_TIME_UNITS)
