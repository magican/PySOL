# -*- coding: utf-8 -*-
"""
cerbere.datamodel.grid
======================

Model class for the grid feature

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import logging
import datetime
import collections

import numpy
from netCDF4 import num2date

from .variable import Variable
from .field import Field

from ..mapper.abstractmapper import AbstractMapper
from .abstractfeature import AbstractFeature

__all__ = ['DEFAULT_PROJECTION', 'Projection', 'Grid']

# Plate Carree projection
DEFAULT_PROJECTION = '+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +units=m'


class Projection(object):
    """Definition of a grid projection"""
    def __init__(self,
                 proj4_definition=DEFAULT_PROJECTION,
                 identifier='regular'
                 ):
        self.definition = proj4_definition
        self.identifier = identifier
        return

    def is_cylindrical(self):
        """return True if the projection is cylindrical"""
        return self.identifier == 'regular'


class Grid(AbstractFeature):
    """
    Model class for the grid feature, ie a two-dimensional array on fixed
    projection, resolution and boundaries
    """
    def __init__(self, identifier=None, title=None,
                description=None, source=None, metadata=None, \
                fields=None, \
                longitudes=None, latitudes=None, depths=None, times=None,
                bbox=None,
                spatial_resolution=None,
                projection=Projection(),
                ):
        """

        Kwargs:
        longitudes (tuple or numpy array or Field):
            Grid longitudes can be initialized :
              * with a tuple (min, max, step) : only for a regular grid. The
               longitudes are automatically generated.
              * a numpy array of float values
              * a Field (for instance copied from another grid object)

        latitudes (tuple or numpy array or Field):
            Grid latitudes can be initialized :
              * with a tuple (min, max, step) : only for a regular grid. The
               latitudes are automatically generated.
              * a numpy array of float values
              * a Field (for instance copied from another grid object)
        """
        # projection parameters
        self.projection = projection
        if not bbox:
            if latitudes is not None and longitudes is not None\
                    and type(latitudes) is tuple\
                    and type(longitudes) is tuple:
                minlon, maxlon, step = longitudes
                minlat, maxlat, step = latitudes
                bbox = (minlon, minlat, maxlon, maxlat)
        if longitudes is not None and type(longitudes) is tuple:
            # generate coordinates
            minlon, maxlon, step = longitudes
            lons = numpy.arange(minlon + step / 2., maxlon, step)
        else:
            lons = longitudes
        if latitudes is not None and type(latitudes) is tuple:
            minlat, maxlat, step = latitudes
            lats = numpy.arange(minlat + step / 2., maxlat, step)
        else:
            lats = latitudes
        super(Grid, self).__init__(
                    identifier=identifier,
                    title=title,
                    description=description,
                    source=source,
                    metadata=metadata,
                    fields=fields,
                    longitudes=lons,
                    latitudes=lats,
                    depths=depths,
                    times=times,
                    bbox=bbox,
                    spatial_resolution=spatial_resolution
                    )
        return

    def is_unique_grid_time(self):
        """
        Return True if a unique time is associated with the grid
        (like in L4 products), False if there is time value per pixel
        (like in L3)
        """
        dims = self.get_geolocation_field('time').dimensions.values()
        return (len(dims) == 1)

    def _get_time_dimensions(self, values):
        """
        return the name and size of each dimension of the `time` geolocation
        field of the feature class.

        Internal function only used when the :class:`Field` instance does not
        exist yet. `values` are used to discriminate between different cases
        for which the number of dimensions for `time` may vary (such as
        :class:`Grid` which can have time defined as a single value or a 2D
        array.

        :param values: time values, as a array of floats/integers or a Field
        object
        :type values: :class:`numpy.ma.MaskedArray`

        :return: dimensions of time field for the feature class
        :rtype: OrderedDict<dimension name, dimension size>
        """
        if len(values.shape) == 1:
            return collections.OrderedDict([('time', 1)])
        elif len(values.shape) == 2:
            nj, ni = values.shape
            return collections.OrderedDict([
                            ('y', nj),
                            ('x', ni)
                            ])
        else:
            raise Exception("Bad shape for time values")

    def get_geolocation_dimnames(self):
        """
        Returns the geolocation dimension names defining a grid
        """
        return ['y', 'x']

    def get_geolocation_field_dimnames(self, fieldname):
        """
        Returns the dimension names of a geolocation field
        """
        if fieldname == 'depth' or fieldname == 'height':
            return ('z',)
        elif fieldname == 'time':
            if self.is_unique_grid_time():
                return ('time',)
            else:
                return ('y', 'x',)
        else:
            if self.projection.is_cylindrical():
                return {'lat': 'y', 'lon': 'x'}[fieldname]
            else:
                return ('y', 'x',)

    def get_geolocation_dimsizes(self):
        """
        Returns the geolocation dimension sizes of the grid object
        """
        return collections.OrderedDict([
            ('y', self.get_geolocation_field('lat').get_dimsize('y')),
            ('x', self.get_geolocation_field('lon').get_dimsize('x'))
            ])

    def save(self, output=None, attrs=None):
        """
        Save the grid to a storage (file,...)

        Args:
            output (:class:`~cerbere.mapper.abstractmapper.AbstractMapper`): storage
                object where to save the feature data.

            attrs (dict): the global metadata (attributes) of the feature, as
                a dictionary where keys are the attributes names.
                See STANDARD_ATTRIBUTES_VALUES in abstractmapper class to see
                a list of standard attributes
        """
        if not output:
            mapper = self.get_mapper()
        else:
            mapper = output
        if mapper.is_writable():
            # creating dimensions
            if self.is_unique_grid_time():
                mapper.create_dim('time', None)
            else:
                mapper.create_dim('time', 1)
            if self.projection.is_cylindrical():
                mapper.create_dim(
                        'lat',
                        self.get_geolocation_field('lat').get_dimsize('y')
                        )
                mapper.create_dim(
                        'lon',
                        self.get_geolocation_field('lon').get_dimsize('x')
                        )
                dim_translation = {'y': 'lat', 'x': 'lon'}
            else:
                mapper.create_dim(
                        'y',
                        self.get_geolocation_field('lat').get_dimsize('y')
                        )
                mapper.create_dim(
                        'x',
                        self.get_geolocation_field('lon').get_dimsize('x')
                        )
                dim_translation = None
            # create additional dimensions
            dims = ['y', 'x', 'lat', 'lon', 'time']
            for v in self._fields.keys():
                if v not in self._geolocation_fields:
                    for d in self._fields[v].dimensions:
                        if d not in dims:
                            output.create_dim(d, self._fields[v].get_dimsize(d))
                            dims.append(d)
            # creating metadata
            if attrs:
                globalattr = attrs
            else:
                globalattr = {}
            if 'title' not in globalattr and self.title is not None:
                globalattr['title'] = self.title
            if 'summary' not in globalattr and self.description is not None:
                globalattr['summary'] = self.description
            if 'id' not in globalattr and self.identifier is not None:
                globalattr['id'] = self.identifier
            lonmin, latmin, lonmax, latmax = self.get_bbox()
            globalattr['geospatial_lat_min'] = latmin
            globalattr['geospatial_lon_min'] = lonmin
            globalattr['geospatial_lat_max'] = latmax
            globalattr['geospatial_lon_max'] = lonmax
            if self.projection.identifier == 'regular':
                if (len(self.get_lat()) > 1) and (len(self.get_lon()) > 1):
                    globalattr['geospatial_lat_resolution'] \
                        = self.get_lat()[1] - self.get_lat()[0]
                    globalattr['geospatial_lon_resolution'] \
                        = self.get_lon()[1] - self.get_lon()[0]
            tmptime = self.get_start_time()
            if tmptime is not None:
                globalattr['time_coverage_start'] = tmptime
            tmptime = self.get_end_time()
            if tmptime is not None:
                globalattr['time_coverage_end'] = tmptime
            globalattr['cdm_data_type'] = 'grid'
            mapper.write_global_attributes(globalattr)
            # creating records
            for geof in self._geolocation_fields:
                if self._geolocation_fields['time'] is None:
                    raise Exception('No time information defined')
                if self._geolocation_fields[geof] is None:
                    logging.warning('Missing geolocation variable : %s', geof)
                else:
                    mapper.create_field(self._geolocation_fields[geof],
                                        dim_translation,
                                        feature='Grid'
                                        )
            for dataf in self._fields.keys():
                if dataf not in self._geolocation_fields:
                    mapper.create_field(self._fields[dataf],
                                        dim_translation,
                                        feature='Grid')
        # saving records
        for geof in self._geolocation_fields:
            field = self._geolocation_fields[geof]
            if field and not field.is_saved():
                mapper.write_field(self._geolocation_fields[geof])
        for dataf in self._fields.keys():
            if not self._fields[dataf].is_saved():
                mapper.write_field(self._fields[dataf])
        mapper.sync()
        return

    def get_spatial_resolution(self):
        """Return the spatial resolution of the feature, in degrees"""
        if self.spatial_resolution is None:
            mapper = self.get_mapper()
            if mapper:
                return mapper.get_spatial_resolution_in_deg()
        else:
            return self.spatial_resolution

    def extract_subset(
            self, boundaries=None, slices=None, fields=None, padding=False):
        """Extract a subset feature from the grid.

        The created subset is a new Grid object without any reference to
        the source.

        Args:
            boundaries (tuple): area of the subset to extract, defined as
                llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat.

            slices (dict): indices for /time/ dimension of the subset to
                extract from the source data. If None, the complete feature
                is extracted.

            fields (list): list of field names to extract. If None, all fields
                are extracted.

            padding (bool): Passed to extract_field method to ensure padding
                with _FillValues for points outside of the bounds of the this
                feature (used only in conjuncture with slices.
        """
        if boundaries and slices:
            raise Exception("Boundaries and slices can not be both provided.")
        if boundaries:
            # get corresponding slices
            llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat = boundaries
            llx1, lly1 = self.latlon2xy(llcrnrlat, llcrnrlon)
            urx2, ury2 = self.latlon2xy(urcrnrlat, urcrnrlon)
            slices = {'y': slice(min(lly1, ury2), max(lly1, ury2), 1),
                      'x': slice(min(llx1, urx2), max(llx1, urx2), 1)}
        if self.is_unique_grid_time():
            timefield = self.extract_field('time')
        else:
            timefield = self.extract_field('time',
                                           slices=slices,
                                           padding=padding)
        subgrid = Grid(
                        latitudes=self.extract_field('lat',
                                                     slices=slices,
                                                     padding=padding),
                        longitudes=self.extract_field('lon',
                                                      slices=slices,
                                                      padding=padding),
                        times=timefield,
                        projection=self.projection,
                        metadata=self.metadata,
#                        bbox=self.get_bbox()
                        )
        if fields is None:
            fields = self.get_fieldnames()
        elif not type(fields) is list:
            raise Exception("fields must be a list")
        for field in fields:
            subgrid.add_field(self.extract_field(field, slices=slices, padding=padding))
        return subgrid

    def extract_spatialsection(self, lat1, lon1, lat2, lon2, fieldnames=None):
        """
        """
        # TBD
        pass

    def latlon2slice(self, lat, lon):
        """Returns the slice corresponding to the provided lat/lon locations

        Args:
            lat (float) : latitude
            lon (float): longitude

        Returns:
            slice
        """
        lats = self.get_lat()
        lons = self.get_lon()
        if self.projection.is_cylindrical():
#             y = lats[numpy.abs(lats - lat).argmin()]
#             x = lons[numpy.abs(lons - lon).argmin()]
            x = numpy.abs(lons - lon).argmin()
            y = numpy.abs(lats - lat).argmin()
#             logging.debug('nearest grid lon %s lat %s',lons[x],lats[y])
            ydim_name_in_attached_storage = self.get_mapper().get_matching_dimname('y')
            xdim_name_in_attached_storage = self.get_mapper().get_matching_dimname('x')
#             logging.debug('latlon2slice dimension name found for y : %s',ydim_name_in_attached_storage)
            return collections.OrderedDict([(ydim_name_in_attached_storage, slice(y,y+1,1)), (xdim_name_in_attached_storage, slice(x,x+1,1))])
        else:
            raise NotImplementedError
