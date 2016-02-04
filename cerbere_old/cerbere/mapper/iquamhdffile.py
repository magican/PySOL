# encoding: utf-8
"""
cerbere.mapper.iquamhdffile
=================================

Mapper class for the iQuam datasets from NOAA / NESDIS / STAR

:copyright: Copyright 2013 Pelamis Scientific Software Ltd.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: David Poulter <david.poulter@pelamis.co.uk>
.. codeauthor:: David Poulter <david.poulter@pelamis.co.uk>
"""

import array
import datetime
from collections import OrderedDict
from itertools import izip

import netCDF4
import numpy
from dateutil.relativedelta import relativedelta

from .. import READ_ONLY, DEFAULT_TIME_UNITS
from ..datamodel.field import Field
from ..datamodel.variable import Variable
from .hdffile import HDFFile

import time

class IQUAMHDFFile(HDFFile):
    """"Cerbere mapper class for NOAA NESDIS STAR iQuam files"""
    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        """Return and instance of IQUAMHDFFile for a givel file.
        
        Accepts:
        
        url
           The location of an input file, which is to be mapped.
           
        Arguments:
        
        mode
           The opening mode for the mapper, defaults to 'r' for 
           read only. Set to 'w' to create a new file with the name
           given in url.
        """
        self.__datetime_cache = {}
        self.__start_time = None
        self.__end_time = None
        super(IQUAMHDFFile, self).__init__(url=url, mode=mode, **kwargs)
        return

    def open(self, datamodel_geolocation_dims=None):
        """
        Open a QuikSCAT HDF file or series of QuikSCAT HDF files and return 
        the handler.
        """
        return super(IQUAMHDFFile, self).open(datamodel_geolocation_dims)

    def get_geolocation_field(self, fieldname):
        # time is a virtual field
        matching = {
            'lat': 'Latitude',
            'lon': 'Longitude',
            'time': 'time',
            'station': 'ID',
            }
        if fieldname in matching:
            return matching[fieldname]
        else:
            return None
            
    def get_dimensions(self, fieldname=None, full=False):
        if fieldname == 'time':
            return ('time',)
        else:
            return super(IQUAMHDFFile, self).get_dimensions(fieldname=fieldname, full=full)

    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature
        '''
        ignore = ['Year', 'Day', 'Month', 'Hour', 'Minute', 'Latitude', 'Longitude']
        fields = [_ for _ in self.get_handler().datasets().keys() if _ not in ignore]
        # remove here time/space information to keep only geophysical fields
        for field in ['time', 'lat', 'lon']:
            if field in fields:
                fields.remove(self.get_geolocation_field(field))
                
        return fields

    def get_matching_dimname(self, geodimname):
        """
        Return the equivalent name in the current data mapper for a standard
        feature dimension, or None if not found.
        """
        # all dimensions are fake in iQuam HDF files, so the dimension name must be
        # assessed from the variable itself.
        matching = {
            'lat': self.get_dimensions(self.get_geolocation_field('lat')),
            'lon': self.get_dimensions(self.get_geolocation_field('lon')),
            'time': self.get_dimensions(self.get_geolocation_field('time')),
            'station': self.get_dimensions(self.get_geolocation_field('station'))
            }
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None
        # return None

    def get_standard_dimname(self, geodimname):
        """Returns the equivalent standard dimension for a geolocation
        dimension in the current data mapper.
        """
        matching = {
            'station': 'fakeDim8',
            'time': 'time',
            'lat': 'lat',
            'lon': 'lon'
            }
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None
        # return None

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
        if native_fieldname == 'time':
            dims = OrderedDict([
                ('time', None),
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
        elif native_fieldname == 'ID':
            dims = OrderedDict([
                ('time', None),
                ('ID_Width', 8)
                ])
            var = Variable(
                            shortname='Platform_ID',
                            description='Platform_ID',
                            authority=self.get_naming_authority(),
                            standardname='Platform_ID'
                            )
            field = Field(
                            var,
                            OrderedDict(dims),
                            datatype=numpy.dtype(numpy.uint8),
                            units=DEFAULT_TIME_UNITS
                            )
            field.attach_storage(self.get_field_handler(fieldname))
            return field
                    
        else:
            dims = OrderedDict([
                ('time', None),
                ])
            oldfield = HDFFile.read_field(self, fieldname)
            field = Field(
                            oldfield.variable,
                            OrderedDict(dims),
                            datatype=oldfield.datatype,
                            units=oldfield.units
                            )
            field.attach_storage(self.get_field_handler(fieldname))
            return field
        # else:
        #     return HDFFile.read_field(self, fieldname)

    def date2num(self, datetime_values, time_units=DEFAULT_TIME_UNITS):
        """Lazy evaluation of netCDF4.date2num."""
        result = []
        for datetime_value in datetime_values:
            if datetime_value in self.__datetime_cache:
                result.append(self.__datetime_cache[datetime_value])
            else:
                value = netCDF4.date2num(datetime_value, DEFAULT_TIME_UNITS)
                self.__datetime_cache[datetime_value] = value
                result.append(value)
        return result
            
    def read_values(self, fieldname, slices=None, indices=None):
        """
        """
        if fieldname == 'time':
            if not slices:
                slices = Ellipsis
            time_fields = ['Year', 'Month', 'Day', 'Hour', 'Minute']
            times = izip(*[self._handler.select(_).get()[slices].tolist() for _ in time_fields])
            datetimes = [datetime.datetime(*_) for _ in times]
            
            #Only calculate the unique datetimes!
            unique_datetimes = numpy.array(list(set(datetimes)))
            unique_numbers = netCDF4.date2num(unique_datetimes, DEFAULT_TIME_UNITS)
            mapper = dict(zip(unique_datetimes, unique_numbers))
            values = numpy.array([mapper.__getitem__(_) for _ in datetimes])
            
        # elif fieldname == 'ID':
        #     values = self._handler.select('ID')[:,:][slices].view('S8')
        # elif fieldname == "corrected_ID":
        #     return self.get_stations
        elif fieldname == 'lon':
            values = HDFFile.read_values(self, fieldname, slices=slices, indices=indices)
            values[values > 180.] -= 360.
        else:
            values = HDFFile.read_values(self, fieldname, slices=slices, indices=indices)

        return values

    def write_field(self, field):
        values = field.get_values()
        if values is None:
            return
        if field.name == 'time':
            if isinstance(values[0], datetime.datetime):
                values = netCDF4.date2num(values, field.units)
        elif field.name == 'Platform_ID':
            values = values.view('unit8')
        logging.debug("Writing field %s", field.name)
        ncvar = self.get_handler().variables[field.name]
        fielddims = ncvar.dimensions
        if ncvar.shape == (0,) \
                 or values.shape == () \
                 or ncvar.shape == values.shape:
            ncvar[:] = values
            if field.qc_levels is not None:
                self.get_handler().variables[field.name + '_qc_level'][:]\
                    = field.qc_levels.values[:]
            if field.qc_details is not None:
                self.get_handler().variables[field.name + '_qc_details'].set_auto_maskandscale(True)
                self.get_handler().variables[field.name + '_qc_details'][:]\
                    = field.qc_details.values[:]
        elif fielddims[0] == 'time' \
                and len(self.get_handler().dimensions['time']) <= 1:
            # case where there is dummy time dimension
            # (for instance for a simple grid or a swath)
            slc = [0]
            for dim in fielddims[1:]:
                slc.append((slice(None, None)))
            ncvar[tuple(slc)] = values
        else:
            raise Exception('Dimensions of fields and netCDF variable do not match')
        return

    # def read_global_attributes(self):
    #     globattr = {}
    #     for attr, val in self.get_handler().attributes().items():
    #         print attr, val
    #         tmp = val.split('\n')
    #         if int(tmp[1]) == 1:
    #             globattr[attr] = tmp[2]
    #         else:
    #             globattr[attr] = []
    #             for item in tmp[2:]:
    #                 globattr[attr].append(item)
    #     return globattr

    def get_start_time(self):
        """Return start of temporal coverage"""
        start_time = str(self.read_global_attribute('START_TIME'))
        if '24:' in start_time:
            result = datetime.datetime.strptime(start_time.replace('24:', '00:'), "%Y-%m-%d %H:%M") 
            return result + relativedelta(days=1)
        return datetime.datetime.strptime(start_time, "%Y-%m-%d %H:%M")
    
    def __get_times_helper(self):
        values = self.read_field('time').get_values(indices={'time':[0, -1]})
        self.__start_time, self.__end_time = netCDF4.num2date(values, DEFAULT_TIME_UNITS)

    def get_end_time(self):
        """Return end of temporal coverage"""
        end_time = str(self.read_global_attribute('END_TIME'))
        if '24:' in end_time:
            result = datetime.datetime.strptime(end_time.replace('24:', '00:'), "%Y-%m-%d %H:%M") 
            return result + relativedelta(days=1)
        return datetime.datetime.strptime(end_time, "%Y-%m-%d %H:%M")    

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

    def get_stations(self, field = 'ID', slices=None, indices=None, cache=True, **kwargs):
        """Return a list with the platform names for each observation. """
        values = self.read_field(field).get_values(indices=indices, slices=slices)
        values = values.view('S8').reshape(values.shape[0]).tolist()        
        return ','.join(values).replace(' ', '').split(',')

    def get_bbox(self):
        lats = self.read_values('lat')
        lons = self.read_values('lon')
        return (lons.min(), lats.min(), lons.max(), lats.max())
