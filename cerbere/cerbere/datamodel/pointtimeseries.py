# -*- coding: utf-8 -*-
"""
cerbere.datamodel.pointtimeseries
=================================

Model class for the time series at a fixed point feature

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import numpy
import collections

from .variable import Variable
from .field import Field
from .abstractfeature import AbstractFeature
from .. import CF_AUTHORITY

__all__ = ['PointTimeSeries']

class PointTimeSeries(AbstractFeature):
    """
    Model class for the point time series
    """
    def __init__(self, identifier=None, title=None,
                description=None, source=None, metadata=None,
                fields=None,
                longitude=None, latitude=None, depth=None, times=None):
        bbox = None
        latitudes = None
        if latitude is not None:
            latitudes = numpy.array([latitude])
        longitudes = None
        if longitude is not None:
            longitudes = numpy.array([longitude])
        depths = None
        if depth is not None:
            depths = numpy.array([depth])
        super(PointTimeSeries, self).__init__(
                    identifier=identifier,
                    title=title,
                    description=description,
                    source=source,
                    metadata=metadata,
                    fields=fields,
                    longitudes=longitudes,
                    latitudes=latitudes,
                    depths=depths,
                    times=times,
                    bbox=bbox
                    )
        return

    def get_geolocation_dimnames(self):
        """
        Returns the geolocation dimension names defining a grid
        """
        return ['time', 'station']

    def get_geolocation_field_dimnames(self, dimname):
        """
        Returns the dimensions of a geolocation field
        """
        if dimname == 'time':
            return ('time',)
        else:
            return ('station',)

    def get_geolocation_dimsizes(self):
        """
        Returns the geolocation dimension sizes of the grid object
        """
        return collections.OrderedDict([('time',
                                         len(self.get_times())),
                                        ('station', 1)])

    def save(self, output, attrs=None):
        """
        Save the time series to a storage (file,...)

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
            output.create_dim('time', None)
            output.create_dim('station', 1)
            dims = ['time', 'station']
            for v in self._fields.keys():
                if not v in self._geolocation_fields:
                    for d in self._fields[v].dimensions:
                        if not d in dims:
                            output.create_dim(d, self._fields[v].get_dimsize(d))
                            dims.append(d)
            # creating metadata
            if attrs:
                globalattr = attrs
            else:
                globalattr = {}
            if self._bbox == None:
                self._bbox = (self.get_lat().min(),
                              self.get_lon().min(),
                              self.get_lat().max(),
                              self.get_lon().max())
            globalattr['geospatial_lat_min'] = self._bbox[0]
            globalattr['geospatial_lon_min'] = self._bbox[1]
            globalattr['geospatial_lat_max'] = self._bbox[2]
            globalattr['geospatial_lon_max'] = self._bbox[3]
            time0 = self.get_datetimes(slices={'time': 0})[0]
            time1 = self.get_datetimes(slices={'time': -1})[0]
            globalattr['time_coverage_start'] = \
                    min(time0, time1).strftime('%Y%m%dT%H%M%S')
            globalattr['time_coverage_end'] = \
                    max(time0, time1).strftime('%Y%m%dT%H%M%S')
            globalattr['cdm_data_type'] = 'station'
            output.write_global_attributes(globalattr)
            # creating records
            for geof in self._geolocation_fields:
                if self._geolocation_fields['time'] is None:
                    raise Exception('No time information defined')
                output.create_field(self._geolocation_fields[geof])
            for dataf in self._fields.keys():
                if not dataf in self._geolocation_fields:
                    output.create_field(self._fields[dataf])
        else:
            raise Exception("Mapper object is not writable")
        # saving records
        for geof in self._geolocation_fields:
            output.write_field(self._geolocation_fields[geof])
        for dataf in self._fields.keys():
            output.write_field(self._fields[dataf])
        output.close()
        return

# The area under this line insluces changes made whilst developing the 
# in situ data server. It is not clear if they are helpful yet, but they 
# probably will be. They are left below for reference, and to be utilised 
# at a later date.
# 
# 
# # -*- coding: utf-8 -*-
# """
# cerbere.datamodel.pointtimeseries
# =================================
# 
# Model class for the time series at a fixed point feature
# 
# :copyright: Copyright 2013 Ifremer / Cersat.
# :license: Released under GPL v3 license, see :ref:`license`.
# 
# .. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
# .. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
# """
# import datetime
# import collections
# 
# import numpy
# import netCDF4
# 
# from .variable import Variable
# from .field import Field
# from .abstractfeature import AbstractFeature
# from .. import CF_AUTHORITY, DEFAULT_TIME_UNITS
# 
# __all__ = ['PointTimeSeries']
# 
# class PointTimeSeries(AbstractFeature):
#     """
#     Model class for the point time series
#     """
#     def __init__(self, identifier=None, title=None,
#                 description=None, source=None, metadata=None,
#                 fields=None,
#                 longitude=None, latitude=None, depth=None, times=None):
#         bbox = None
#         latitudes = None
#         if latitude is not None:
#             latitudes = numpy.array([latitude])
#         longitudes = None
#         if longitude is not None:
#             longitudes = numpy.array([longitude])
#         depths = None
#         if depth is not None:
#             depths = numpy.array([depth])
#         super(PointTimeSeries, self).__init__(
#                     identifier=identifier,
#                     title=title,
#                     description=description,
#                     source=source,
#                     metadata=metadata,
#                     fields=fields,
#                     longitudes=longitudes,
#                     latitudes=latitudes,
#                     depths=depths,
#                     times=times,
#                     bbox=bbox
#                     )
#         return
# 
#     def extract_subset(self, boundaries=None, slices=None):
#         """
#         extract a subset from the trajectory. The created subset is a new object without any reference to the source.
# 
#         Args
#         ----
#         boundaries (tuple)
#             area of the subset to extract, defined as llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat
# 
#         slices (dictionary of slices)
#             indices for /time/ dimension of the subset to extract from the source data
# 
#         """
#         if boundaries and slices:
#             raise Exception("Boundaries and slices can not be both provided.")
#         if boundaries:
#             # get corresponding slices
#             raise Exception("Not yet implemented")
#         # elif slices:
#             # if len(slices) != 1:
#             #     raise Exception ("only slices for the time dimension must be provided")
#         subtraj = PointTimeSeries(
#             latitude=self.extract_field('lat', slices=slices),
#             longitude=self.extract_field('lon', slices=slices),
#             times = self.extract_field('time', slices=slices), 
#             metadata = self.metadata)
#         for field in self.get_fieldnames():
#             subtraj.add_field( self.extract_field(field, slices=slices) )
#         return subtraj
# 
#     def get_datetimes(self, slices=None, indices=None, cache=True, **kwargs):
#         '''
#         return the time values as datetime objects (instead of internal
#          numerical representation)
# 
#         Kwargs
#         slices: array of slice (example: [slice(0,10,1),slice(0,10,1)] )
#             array of the slices to extract from the feature. slices and
#              indices can not be both specified.
# 
#         indices: array of integers
#             array indices of a single datetime value to extract from the swath
# 
#         cache: bool
#             if True, cache the result in memory when reading the data from a
#              file. Otherwise the data are read again each time the method is
#               called.
# 
#         Beware:
#         -------
#         Conversion to datetime objects can be very slow
#         '''
#         field = self._geolocation_fields['time']
#         values = self._get_field_values(field, slices, indices, cache, **kwargs)
# 
#         if not isinstance(values, datetime.datetime)\
#                 or not isinstance(values[0], datetime.datetime):
# 
#             #Enhancment, speeds up iQuam by ~80%
#             unique_numbers = numpy.unique(values)
#             unique_datetimes = netCDF4.num2date(unique_numbers, DEFAULT_TIME_UNITS)
#             mapper = dict(zip(unique_numbers, unique_datetimes))
#             values = numpy.vectorize(mapper.get)(values)
#         return values
# 
#     def get_geolocation_dimnames(self):
#         """
#         Returns the geolocation dimension names defining a grid
#         """
#         return ['time', 'station']
# 
#     def get_geolocation_field_dimnames(self, dimname):
#         """
#         Returns the dimensions of a geolocation field
#         """
#         if dimname == 'time':
#             return ('time',)
#         else:
#             return ('station',)
# 
#     def get_geolocation_dimsizes(self):
#         """
#         Returns the geolocation dimension sizes of the grid object
#         """
#         return collections.OrderedDict([('time',
#                                          len(self.get_times())),
#                                         ('station', 1)])
# 
#     def save(self, output, attrs=None):
#         """
#         Save the grid to a storage (file,...)
# 
#         Args:
#             output (instance of a `mapper` class)
#                 storage object where to save the feature data
# 
#         Kwargs:
#             attrs:dictionary
#                 specifies the global metadata (attributes) 
#                 of the feature. 
#                 See STANDARD_ATTRIBUTES_VALUES in abstractmapper class to see
#                 a list of standard attributes
#         """
#         if output.is_writable():
#             # creating dimensions
#             output.create_dim('time', None)
#             output.create_dim('station', 1)
#             dims = ['time', 'station']
#             for v in self._fields.keys():
#                 if not v in self._geolocation_fields:
#                     for d in self._fields[v].dimensions:
#                         if not d in dims:
#                             output.create_dim(d, self._fields[v].get_dimsize(d))
#                             dims.append(d)
#             # creating metadata
#             if attrs:
#                 globalattr = attrs
#             else:
#                 globalattr = {}
#             if self._bbox == None:
#                 self._bbox = (self.get_lat().min(),
#                               self.get_lon().min(),
#                               self.get_lat().max(),
#                               self.get_lon().max())
#             globalattr['geospatial_lat_min'] = self._bbox[0]
#             globalattr['geospatial_lon_min'] = self._bbox[1]
#             globalattr['geospatial_lat_max'] = self._bbox[2]
#             globalattr['geospatial_lon_max'] = self._bbox[3]
#             time0 = self.get_datetimes(slices={'time': 0})[0]
#             time1 = self.get_datetimes(slices={'time': -1})[0]
#             globalattr['time_coverage_start'] = \
#                     min(time0, time1).strftime('%Y%m%dT%H%M%S')
#             globalattr['time_coverage_stop'] = \
#                     max(time0, time1).strftime('%Y%m%dT%H%M%S')
#             globalattr['cdm_feature_type'] = 'station'
#             # output.write_global_attributes(globalattr)
#             # creating records
#             for geof in self._geolocation_fields:
#                 if self._geolocation_fields['time'] is None:
#                     raise Exception('No time information defined')
#                 output.create_field(self._geolocation_fields[geof])
#             for dataf in self._fields.keys():
#                 if not dataf in self._geolocation_fields:
#                     output.create_field(self._fields[dataf])
#         else:
#             raise Exception("Mapper object is not writable")
#         # saving records
#         for geof in self._geolocation_fields:
#             output.write_field(self._geolocation_fields[geof])
#         for dataf in self._fields.keys():
#             output.write_field(self._fields[dataf])
#         output.close()
#         return
