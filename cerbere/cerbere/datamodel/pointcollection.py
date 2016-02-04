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
    Implements a set of randomly sampled points
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

    def extract_subset(self, boundaries=None, slices=None, indices=None):
        """Extract a subset from the trajectory. The created subset is a new
        object without any reference to the source.

        Args:
            boundaries (tuple, optional): area of the subset to extract,
                defined as llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat.
            slices (dict): dictionary of slices for `time` dimension of the
                subset to extract from the source data
        """
        if boundaries and slices:
            raise ValueError("Boundaries and slices can not be both provided.")
        if boundaries:
            # get corresponding slices
            raise NotImplementedError()
        if not indices:
            subcoll = PointCollection(
                latitudes=self.extract_field('lat', slices=slices),
                longitudes=self.extract_field('lon', slices=slices),
                times=self.extract_field('time', slices=slices),
                metadata=self.metadata)
            for field in self.get_fieldnames():
                subcoll.add_field(self.extract_field(field, slices=slices))
        else:
            subcoll = PointCollection(
                latitudes=self.extract_field('lat', indices=indices),
                longitudes=self.extract_field('lon', indices=indices),
                times=self.extract_field('time', indices=indices),
                metadata=self.metadata)
            for field in self.get_fieldnames():
                subcoll.add_field(self.extract_field(field, indices=indices))
        return subcoll

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

    def get_geolocation_field_dimnames(self, fieldname):
        """
        Returns the dimensions of a geolocation field
        """
        return ('time',)

    def save(self, output, attrs={}):
        """
        Save the point collection to a storage mapper (file,...)

        Args:
            output (:class:`~cerbere.mapper.abstractmapper.AbstractMapper`):
                storage object where to save the feature data.

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
                if v not in self._geolocation_fields:
                    for d in self._fields[v].dimensions:
                        if d not in dims:
                            output.create_dim(d,
                                              self._fields[v].get_dimsize(d))
                            dims.append(d)

            # creating metadata
            if attrs:
                globalattr = attrs
            else:
                globalattr = {}
            lonmin, latmin, lonmax, latmax = self.get_bbox()
            globalattr['geospatial_lat_min'] = latmin
            globalattr['geospatial_lon_min'] = lonmin
            globalattr['geospatial_lat_max'] = latmax
            globalattr['geospatial_lon_max'] = lonmax
            tmptime = self.get_start_time()
            if tmptime is not None:
                globalattr['time_coverage_start'] = tmptime
            tmptime = self.get_end_time()
            if tmptime is not None:
                globalattr['time_coverage_end'] = tmptime
            globalattr['cdm_data_type'] = 'point'
            if 'title' not in globalattr and self.title is not None:
                globalattr['title'] = self.title
            if 'summary' not in globalattr and self.description is not None:
                globalattr['summary'] = self.description
            if 'id' not in globalattr and self.identifier is not None:
                globalattr['id'] = self.identifier
            output.write_global_attributes(globalattr)

            # create records
            for geof in self._geolocation_fields:
                if not self._geolocation_fields['time']:
                    raise Exception('No time information defined')
                if not self._geolocation_fields[geof]:
                    logging.debug('Missing geolocation variable : %s', geof)
                else:
                    output.create_field(self._geolocation_fields[geof])
            for v in self._fields.keys():
                if v not in self._geolocation_fields:
                    output.create_field(self._fields[v])
        # saving records
        for geof in self._geolocation_fields:
            if self._geolocation_fields[geof]:
                output.write_field(self._geolocation_fields[geof])
        for dataf in self._fields.keys():
            output.write_field(self._fields[dataf])
        output.close()
