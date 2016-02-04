#!/usr/bin/env python
# -*- coding: utf-8 -*-


import logging
import datetime
from collections import OrderedDict

import netCDF4
import numpy

from .abstractfeature import AbstractFeature
from ..mapper import abstractmapper
from .variable import Variable
from .field import Field
from .. import CF_AUTHORITY, DEFAULT_TIME_UNITS

__all__ = ['PointCollection']


class PointCollection(AbstractFeature):
    """
    Implements a set of randomly sample points
    """
    def __init__(self, identifier=None, title=None, description=None,
                 source=None, metadata=None,
                fields=None,
                longitudes=None, latitudes=None, depths=None, times=None):
        super(PointCollection, self).__init__(
                        identifier=identifier,
                        title=title,
                        description=description,
                        source=source,
                        metadata=metadata,
                        fields=fields,
                        longitudes=longitudes,
                        latitudes=latitudes,
                        depths=depths,
                        times=times)
        return

    def get_datetimes(self, slices=None, indices=None, cache=True, **kwargs):
        '''
        return the time values as datetime objects (instead of internal
         numerical representation)

        Kwargs
        slices: array of slice (example: [slice(0,10,1),slice(0,10,1)] )
            array of the slices to extract from the feature. slices and
             indices can not be both specified.

        indices: array of integers
            array indices of a single datetime value to extract from the swath

        cache: bool
            if True, cache the result in memory when reading the data from a
             file. Otherwise the data are read again each time the method is
              called.

        Beware:
        -------
        Conversion to datetime objects can be very slow
        '''
        field = self._geolocation_fields['time']
        values = self._get_field_values(field, slices, indices, cache, **kwargs)
        
        if not isinstance(values, datetime.datetime)\
                or not isinstance(values[0], datetime.datetime):

            #Enhancment, speeds up iQuam by ~80%
            unique_numbers = numpy.unique(values)
            if isinstance(unique_numbers, numpy.ma.core.MaskedArray):
                unique_datetimes = netCDF4.num2date(unique_numbers.data, DEFAULT_TIME_UNITS)
            else:
                unique_datetimes = netCDF4.num2date(unique_numbers, DEFAULT_TIME_UNITS)
            mapper = dict(zip(unique_numbers, unique_datetimes))
            values = [mapper.__getitem__(_) for _ in values]
        return values
        
    def get_stations(self, field = 'Platform_ID', slices=None, indices=None, cache=True, format=None, **kwargs):
        """Return a list with the platform names for each observation. """
        field = self.get_field(field)
        values = self._get_field_values(field, slices, indices, cache, **kwargs)
        # iQuam station names are stored as ASCII (uint8) values in an N*8 array.
        # This code converts those values to a list of strings
        if format == 'iQuam':
            values = values.view('S8').reshape(values.shape[0]).tolist()        
            values = ','.join(values).replace(' ', '').split(',')
        return values

    def extract_subset(self, boundaries=None, slices=None, indices=None):
        """
        extract a subset from the trajectory. The created subset is a new object without any reference to the source.

        Args
        ----
        boundaries (tuple)
            area of the subset to extract, defined as llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat

        slices (dictionary of slices)
            indices for /time/ dimension of the subset to extract from the source data

        """
        if boundaries and slices:
            raise Exception("Boundaries and slices can not be both provided.")
        if boundaries:
            # get corresponding slices
            raise Exception("Not yet implemented")
        if not indices:
            subtraj = PointCollection(
                latitudes=self.extract_field('lat', slices=slices),
                longitudes=self.extract_field('lon', slices=slices),
                times = self.extract_field('time', slices=slices), 
                metadata = self.metadata)
            for field in self.get_fieldnames():
                subtraj.add_field( self.extract_field(field, slices=slices) )
        else:
            subtraj = PointCollection(
                latitudes=self.extract_field('lat', indices=indices),
                longitudes=self.extract_field('lon', indices=indices),
                times = self.extract_field('time', indices=indices), 
                metadata = self.metadata)
            for field in self.get_fieldnames():
                subtraj.add_field( self.extract_field(field, indices=indices) )           
        return subtraj

    def get_geolocation_dimnames(self):
        """
        Returns the geolocation dimension names defining a point collection
        """
        return ['time']

    def get_geolocation_dimsizes(self):
        """
        Returns the geolocation dimension sizes of the grid object
        """
        return OrderedDict([('time',
                             len(self.get_times()))])

    def get_geolocation_field_dimnames(self, dimname):
        """
        Returns the dimensions of a geolocation field
        """
        return ('time',)

    def save(self, output, attrs={}):
        """
        Save the point collection to a storage (file,...)

        Args:
            output (:class:`~cerbere.mapper.abstractmapper.AbstractMapper`): storage
                object where to save the feature data.

            attrs (dict): the global metadata (attributes) of the feature, as
                a dictionary where keys are the attributes names.
                See STANDARD_ATTRIBUTES_VALUES in abstractmapper class to see
                a list of standard attributes
        """
        if output.is_writable():
            # creating dimensions
            output.create_dim('time', self.get_times().size)
            dims = ['time']
            for v in self._fields.keys():
                if not v in self._geolocation_fields:
                    for d in self._fields[v].dimensions:
                        if not d in dims:
                            output.create_dim(d, self._fields[v].get_dimsize(d))
                            dims.append(d)
            # creating metadata
            globalattr = attrs
            
            self._bbox = self.get_bbox()
            globalattr['geospatial_lat_min'] = self._bbox[1]
            globalattr['geospatial_lon_min'] = self._bbox[0]
            globalattr['geospatial_lat_max'] = self._bbox[3]
            globalattr['geospatial_lon_max'] = self._bbox[2]
            globalattr['time_coverage_start'] = self.get_datetimes()[0].strftime('%Y-%m-%dT%H:%M:%SZ')
            globalattr['time_coverage_end'] = self.get_datetimes()[-1].strftime('%Y-%m-%dT%H:%M:%SZ')
            globalattr['cdm_data_type'] = 'point'
            output.write_global_attributes(globalattr)
            
            # if self._bbox == None:
            #     self._bbox = self.get_bbox()
            # globalattr['geospatial_lat_min'] = self._bbox[1]
            # globalattr['geospatial_lon_min'] = self._bbox[0]
            # globalattr['geospatial_lat_max'] = self._bbox[3]
            # globalattr['geospatial_lon_max'] = self._bbox[2]
            # if not 'time_coverage_start' in attrs:
            #     if type(self.get_datetimes()) == type(datetime.datetime.now()):
            #         globalattr['time_coverage_start']\
            #             = self.get_datetimes().strftime('%Y-%m-%dT%H:%M:%SZ')
            #     else:
            #         globalattr['time_coverage_start']\
            #             = self.get_datetimes()[0].strftime('%Y-%m-%dT%H:%M:%SZ')
            # if not 'time_coverage_stop' in attrs:
            #     if type(self.get_datetimes()) == type(datetime.datetime.now()):
            #         globalattr['time_coverage_stop']\
            #             = self.get_datetimes().strftime('%Y-%m-%dT%H:%M:%SZ')
            #     else:
            #         globalattr['time_coverage_stop']\
            #             = self.get_datetimes()[-1].strftime('%Y-%m-%dT%H:%M:%SZ')
            # globalattr['cdm_feature_type'] = 'trajectory'
            # output.write_global_attributes(globalattr)
            # creating records
            for geof in self._geolocation_fields:
                if not self._geolocation_fields['time']:
                    raise Exception('No time information defined')
                if not self._geolocation_fields[geof]:
                    logging.debug('Missing geolocation variable : %s', geof)
                else:
                    output.create_field(self._geolocation_fields[geof])
            for v in self._fields.keys():
                if not v in self._geolocation_fields:
                    output.create_field(self._fields[v])
        # saving records
        for geof in self._geolocation_fields:
            if self._geolocation_fields[geof]:
                output.write_field(self._geolocation_fields[geof])
        for dataf in self._fields.keys():
            output.write_field(self._fields[dataf])
        output.close()
        return

