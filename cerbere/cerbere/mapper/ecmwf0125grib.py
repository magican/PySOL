# -*- coding: utf-8 -*-
"""
cerbere.mapper.gribfile
=======================

Mapper classs for grib format

:copyright: Copyright 2014 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Antoine Grouazel <antoine.grouazel@ifremer.fr>
.. codeauthor:: Antoine Grouazel <antoine.grouazel@ifremer.fr>
"""

import os
import logging
import datetime
import copy
import collections
from collections import OrderedDict

import numpy
import pygrib
from netCDF4 import date2num

from .. import READ_ONLY, WRITE_NEW, DEFAULT_TIME_UNITS
from .gribfile import GribFile
from ..datamodel.field import Field
from ..datamodel.variable import Variable

VIRTUALFIELD_DESCR = {'wind_speed': 'wind speed',
                      'wind_direction': 'wind direction (meteorological convention)'
                      }
VIRTUALFIELD_STDNAME = {'wind_speed': 'wind_speed',
                      'wind_direction': 'wind_to_direction'
                      }
VIRTUALFIELD_UNITS = {
                'wind_speed': 'm s-1',
                'wind_direction': 'degrees'
                }
class ECMWF0125Grib(GribFile):
    """Mapper class to read ECMWF 0125 grib files"""

    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        """
        Initialize a Grib file mapper
        """
        super(ECMWF0125Grib, self).__init__(url=url, mode=mode, **kwargs)
        self._messages = collections.OrderedDict()
        self._geodim_model2mapper= collections.OrderedDict()
        self._geodim_mapper2model = collections.OrderedDict()
        return

        
    def read_values(self, fieldname, slices=None):
        """
        """
        if fieldname == 'lat':
            firstfield = self._messages[self._messages.keys()[0]][0]
            if firstfield.gridType != 'regular_ll':
                values,junklon = self._messages['10v'][0].latlons()
            else:
                values = firstfield.distinctLatitudes
        elif fieldname == 'lon':
            firstfield = self._messages[self._messages.keys()[0]][0]
            if firstfield.gridType != 'regular_ll':
                junklat,values = self._messages['10v'][0].latlons()
                
            else:
                values = firstfield.distinctLongitudes
#             values = numpy.roll(values, values.shape[0] / 2, axis=0)
            over180 = values[:,values.shape[1] / 2 :]
            under180 = values[:,0:values.shape[1] / 2]
            values = numpy.concatenate((over180,under180), axis=1)
            ind = numpy.ma.where(values >= 180.)
            values[ind] = values[ind] - 360.
        elif fieldname == 'time':
            firstfield = self._messages[self._messages.keys()[0]][0]
#             values = numpy.array([date2num(firstfield.analDate,DEFAULT_TIME_UNITS)])
            ddate = str(firstfield.dataDate)
            dtime = str(firstfield.dataTime)
            completdate = datetime.datetime.strptime(ddate+dtime,'%Y%m%d%H%M')
            print ddate,type(ddate),dtime,completdate
            values = numpy.array([date2num(completdate,DEFAULT_TIME_UNITS)])
            return values
        elif fieldname == "wind_speed":
            u10 = self.read_values('10u')
            v10 = self.read_values('10v')
            u10=u10[slices]
            v10=v10[slices]
            u = numpy.sqrt(u10 * u10 + v10 * v10)
            shaap=numpy.shape(u10)
            u = numpy.reshape(u,(shaap[1],shaap[2]))
            return u
        elif fieldname == "wind_direction":
            u10 = self.read_values('10u')
            v10 = self.read_values('10v')
            u10=u10[slices]
            v10=v10[slices]
            dir = uv2dir(u10, v10)
            shaap=numpy.shape(u10)
            dir = numpy.reshape(dir,(shaap[1],shaap[2]))
            return dir
        else:
            fieldmsg = self._messages[fieldname][0]
            if len(self._messages[fieldname]) > 1:
                # levels to be implemented
                raise NotImplementedError
            values = fieldmsg.values
            values = numpy.roll(values, values.shape[1] / 2, axis=1)
        if slices:
            return values[slices]
        else:
            return values

    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        # Variable
        if not fieldname in ['time', 'lat', 'lon']:
            attrs = self.read_field_attributes(fieldname)
        else:
            attrs = {}
        if 'name' in attrs:
            descr = attrs['name']
        else:
            descr = None
        var = Variable(shortname=fieldname, description=descr)
        dims = collections.OrderedDict()
        if not fieldname in ['time', 'lat', 'lon']:
            if attrs['gridType'] == 'regular_ll':
                dims['time'] = 1
                dims['y'] = attrs['Nj']
                dims['x'] = attrs['Ni']
            elif attrs['gridType'] == 'reduced_gg':
                dims['time'] = 1
                dims['y'] = attrs['Nj']
                dims['x'] = attrs['Nj']*2 #works fine for N640 transformed in regular grid
            else:
                raise NotImplementedError
            if fieldname in ['wind_speed', 'wind_direction']:
                # create a virtual field
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
                        units=VIRTUALFIELD_UNITS[fieldname]
                        )
                field.attach_storage(self.get_field_handler(fieldname))
            else:   #generic case
                
                datatype = self._messages[fieldname][0]['values'].dtype
                field = Field(var, dims, datatype=datatype)
                field.attach_storage(self.get_field_handler(fieldname))
                # MetaData
                field.attributes = {}
                field.units = None
                if 'units' in attrs:
                    field.units = attrs['units']
                field.valid_min = None
                field.valid_max = None
                if 'minimum' in attrs:
                    field.valid_min = attrs['minimum']
                if 'maximum' in attrs:
                    field.valid_min = attrs['maximum']
                return field
        elif fieldname == 'time':
            dims['time'] = 1
            units='standard gregorian date'
            field = Field(var, dims, datatype=numpy.dtype(numpy.float64))
            field.attach_storage(self.get_field_handler(fieldname))
            field.units=units
        else:
            firstfield = self._messages[self._messages.keys()[0]][0]
            logging.debug('firstfield.gridType: %s',firstfield.gridType)
#             if firstfield.gridType == 'regular_ll':
            #get ride of the difference of grid type since pygrib abstract this directly
            if fieldname == 'lat':
                dims['y'] = firstfield.Nj
                units = 'degrees_north'
                var.description = 'latitude'
            elif fieldname == 'lon':
                dims['x'] = firstfield.Ni
                units = 'degrees_east'
                var.description = 'longitude'
            field = Field(
                        var,
                        dims,
                        datatype=numpy.dtype(numpy.float64),
                        units=units)

            field.attach_storage(self.get_field_handler(fieldname))
#             else:
#                 raise NotImplementedError
        return field

    def get_fieldnames(self):
        """
        """
        fieldnames = self._messages.keys()
        fieldnames.extend(["wind_speed", "wind_direction","lon","lat"])
        return fieldnames


    def read_field_attributes(self, fieldname):
        """
        """
        if fieldname == 'wind_speed' or fieldname == 'wind_direction':
            fieldname = '10u'
        if fieldname in self._messages:
            attrs = {}
            for k in self._messages[fieldname][0].keys():
                if not k in [
                        'values',
                        'distinctLongitudes',
                        'distinctLatitudes',
                        'latitudes',
                        'longitudes',
                        'latLonValues',
                        'codedValues',
                        'analDate',
                        'validDate'
                        ]:
                    
                    attrs[k] = self._messages[fieldname][0][k]
        else:
            logging.error('there is no %s message in this grib file ',fieldname)
            raise NotImplementedError
        return attrs

    def get_bbox(self):
        """returns the bounding box of the feature

        Return:
        tuple (lonmin, latmin, lonmax, latmax)
        """
        return (-180., -90., 180., 90.)
    def get_spatial_resolution_in_deg(self):
        """Returns the average spatial resolution in degrees"""
        return None
    
    def get_standard_dimname(self,geodimname):
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
                'time': 'time',
                'lon': 'x',
                'longitude': 'x',
                'lat': 'y',
                'latitude': 'y',
                'mes': 'station',
                'station': 'station',
                'ni': 'cell',
                'cell': 'cell',
                'ra_size': 'cell',
                'nj': 'row',
                'row': 'row',
                'az_size': 'row',
                'depth': 'z',
                'height': 'z'
                }
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None
