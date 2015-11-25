# encoding: utf-8
"""
cerbere.mapper.qscathdffile
=================================

Mapper class for the QuikSCAT files from PODAAC

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import logging
import array
import datetime
from collections import OrderedDict

import netCDF4
import numpy

from .. import READ_ONLY, DEFAULT_TIME_UNITS
from ..datamodel.field import Field
from ..datamodel.variable import Variable
from .hdffile import HDFFile


class QSCATHDFFile(HDFFile):
    """"mapper class for PODAAC QuikSCAT HDF files"""
    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        """
        """
        HDFFile.__init__(self, url=url, mode=mode, **kwargs)
        return

    def open(self, datamodel_geolocation_dims=None):
        """
        Open a QuikSCAT HDF file or series of QuikSCAT HDF files and return the
        handler.
        """
        res = super(QSCATHDFFile, self).open(datamodel_geolocation_dims)
        if self._mode == READ_ONLY:
            if not self._url is None:
                latvalues = HDFFile.read_values(self, 'lat')
                latvalues = numpy.ma.masked_equal(latvalues, 0)
                latshape = latvalues.shape
                # QuikSCAT lat/lon contains 0 filled rows that must be
                # trimmed
                #self.__trim_nan_rows(latvalues)
                # QuikSCAT lat/lon contains unfilled columns that must be
                # trimmed
                #self.__trim_nan_cols(latvalues)
        return res

    def get_geolocation_field(self, fieldname):
        """
        return the equivalent field name in the file format for a standard
        geolocation field (lat, lon, time).

        Used for internal purpose and should not be called directly.

        :param fieldname: name of the standard geolocation field (lat, lon
            or time)
        :type fieldname: str

        :rtype: str or None
        :return: name of the corresponding field in the native file format.
            Returns None if no matching is found
        """
        matching = {'lat': 'wvc_lat',
                    'lon': 'wvc_lon',
                    'time': 'wvc_row_time'
                    }
        if fieldname in matching:
            return matching[fieldname]
        else:
            return None

    def get_matching_dimname(self, geodimname):
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
                'row': 'Wind_Vector_Cell_Row',
                'cell': 'Wind_Vector_Cell',
                'time': 'wvc_row_time'
                }
        if geodimname in matching:
            return matching[geodimname]
        else:
            return geodimname

    def get_standard_dimname(self, geodimname):
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
                'Wind_Vector_Cell_Row': 'row',
                'Wind_Vector_Cell': 'cell',
                'wvc_row_time': 'time'
                }
        if geodimname in matching:
            return matching[geodimname]
        else:
            return geodimname

    def __trim_nan_rows(self, data):
        '''
        trim the nan filled rows from the file content.
        Introduced to remove unfilled rows from QuiSCAT L2B files.
        This method is called once.

        The file content accessed by any method (such as read_values, etc...)
        will then be the trimmed arrays of values for any variable of the file.

        Warning
        -------
        This method can only be called on 2 dimensional variables (lat,lon)
        '''
        nbrows = data.shape[0]
        k = nbrows - 1
        while data.mask[k, :].all() == True:
            k = k - 1
        self.max_valid_row = k
        k = 0
        while data.mask[k, :].all() == True:
            k = k + 1
        self.min_valid_row = k
        return

    def __trim_nan_cols(self, data):
        '''
        trim the nan filled columns from the file content.
        Introduced to remove unfilled columns from QuiSCAT L2B files.
        This method is called once.

        The file content accessed by any method (such as read_values, etc...)
        will then be the trimmed arrays of values for any variable of the file.

        Warning
        -------
        This method can only be called on 2 dimensional variables (lat,lon)
        '''
        nbcols = data.shape[-1]
        k = nbcols - 1
        while data.mask[:, k].all() == True:
            k = k - 1
        self.max_valid_col = k
        k = 0
        while data.mask[:, k].all() == True:
            k = k + 1
        self.min_valid_col = k
        return

    def get_file_offset(self, x, axis):
        '''
        returns the true offset in a file of an indice in a given axis
        (dimension indice).
        This is provided when one want to locate a value in file, since
        the arrays of values returned by get_values may be trimmed of invalid
        rows in the actual file content.
        '''
        if axis == 0:
            return x + self.min_valid_col
        if axis == 1:
            return x + self.min_valid_row

    def get_valid_rows(self):
        """
        """
        return (self.min_valid_row, self.max_valid_row)

    def get_valid_cols(self):
        """
        """
        return (self.min_valid_col, self.max_valid_col)

    def read_field(self, fieldname):
        """
        """
        native_fieldname = self.get_geolocation_field(fieldname)
        if native_fieldname is None:
            native_fieldname = fieldname
        if native_fieldname == 'wvc_row_time':
            dims = OrderedDict([
                ('row', self.get_dimsize('Wind_Vector_Cell_Row')),
                ('cell', self.get_dimsize('Wind_Vector_Cell')),
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
                            datatype=numpy.dtype(numpy.float32),
                            units=DEFAULT_TIME_UNITS
                            )
            field.attach_storage(self.get_field_handler(fieldname))
            return field
        else:
            return HDFFile.read_field(self, fieldname)

#     def __fill_missing_values(self, data):
#         '''
#         Interpolate missing values
#         '''
#         slice_col = slice(self.min_valid_col, self.max_valid_col + 1, 1)
#         slice_row = slice(self.min_valid_row, self.max_valid_row + 1, 1)
#         valid_data = data[slice_row, slice_col].filled(numpy.NaN)
#         bad_data_col = numpy.where(numpy.isnan(valid_data).any(0) == True)
#         indices = numpy.arange(valid_data.shape[0])
#         for i in bad_data_col[0]:
#             nan = numpy.isnan(valid_data[:, i]).ravel()
#             #print i, nan
#             if len(nan) > 0:
#                 not_nan = ~nan
#                 finterp = InterpolatedUnivariateSpline(indices[not_nan],
#                                                        valid_data[:, i][not_nan],
#                                                        k=1)
#                 #print finterp
#                 tmp = finterp(indices[nan])
#                 if not isinstance(tmp, numpy.ndarray):
#                     tmp = numpy.array([tmp])
#                 valid_data[indices[nan], i] = tmp[:]
#         valid_data = numpy.ma.masked_where(numpy.isnan(valid_data), valid_data)
#         data[slice_row, slice_col] = valid_data
#         return data

    def read_values(self, fieldname, slices=None):
        """
        """
        std_fieldname = self.get_standard_dimname(fieldname)
        if std_fieldname == 'time':
            vd = self._vdata.attach('wvc_row_time')
            values = [datetime.datetime.strptime(
                                        array.array('B', v[0]).tostring(),
                                        '%Y-%jT%H:%M:%S.%f'
                                        ) for v in vd[:]]
            nbrows, nbcells =\
                self.get_handler().select('wvc_lat').dimensions().values()
            values = numpy.array([netCDF4.date2num(values,
                                                   DEFAULT_TIME_UNITS)])
            values = numpy.resize(values, (nbcells, nbrows,)).transpose()
            if slices:
                values = values[slices]
        else:
            native_fieldname = self.get_matching_dimname(fieldname)
            values = HDFFile.read_values(self, native_fieldname, slices)
        if std_fieldname == 'lon':
            ind = numpy.ma.where(values >= 180.)
            if len(ind[0]) > 0:
                values[ind] = values[ind] - 360.
        if std_fieldname == 'lon' or std_fieldname == 'lat':
            if std_fieldname == 'lon':
                values_bis = HDFFile.read_values(self, 'lat', slices)
            elif std_fieldname == 'lat':
                values_bis = HDFFile.read_values(self, 'lon', slices)
            # mask null values
            values = numpy.ma.masked_where(
                            (values == 0.) & (values_bis == 0.),
                            values)
            #values = self.__fill_missing_values(values)
            if std_fieldname == 'lat':
                ind = numpy.ma.where(values > 90.)
                values[ind] = 90.
                ind = numpy.ma.where(values < -90.)
                values[ind] = -90.
            elif std_fieldname == 'lon':
                ind = numpy.ma.where(values > 180.)
                if len(ind[0]) > 0:
                    values[ind] = values[ind] - 360.
                ind = numpy.ma.where(values < -180.)
                if len(ind[0]) > 0:
                    values[ind] = values[ind] + 360.
        return values

    def read_global_attributes(self):
        globattr = {}
        for attr, val in self.get_handler().attributes().items():
            tmp = val.split('\n')
            if int(tmp[1]) == 1:
                globattr[attr] = tmp[2]
            else:
                globattr[attr] = []
                for item in tmp[2:]:
                    globattr[attr].append(item)
        return globattr

    def get_start_time(self):
        """Return start of temporal coverage"""
        startdate = str(self.read_global_attribute('RangeBeginningDate'))
        starttime = str(self.read_global_attribute('RangeBeginningTime'))
        start = datetime.datetime.strptime(startdate + ' ' + starttime,
                                           '%Y-%j %H:%M:%S.%f'
                                           )
        return start

    def get_end_time(self):
        """Return end of temporal coverage"""
        enddate = str(self.read_global_attribute('RangeEndingDate'))
        endtime = str(self.read_global_attribute('RangeEndingTime'))
        end = datetime.datetime.strptime(enddate + ' ' + endtime,
                                           '%Y-%j %H:%M:%S.%f'
                                           )
        return end

    def get_collection_id(self):
        """return the identifier of the product collection"""
        return str(self.read_global_attribute('ShortName'))

    def get_product_version(self):
        """return the product version"""
        return str(self.read_global_attribute('build_id'))

    def get_orbit_number(self):
        """return the orbit number"""
        return str(self.read_global_attribute('rev_number'))

    def get_cycle_number(self):
        """return the cycle number"""
        return None

    def get_bbox(self):
        lats = self.read_values('lat')
        lons = self.read_values('lon')
        return (lons.min(), lats.min(), lons.max(), lats.max())
