# -*- coding: utf-8 -*-
"""
cerbere.datamodel.gridtimeseries
================================

Model class for the time series of grid feature

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import logging
import copy
import collections
import gc
import pdb
import numpy
from netCDF4 import num2date, date2num

from .abstractfeature import AbstractFeature
from .grid import Projection, Grid

__all__ = ['GridTimeSeries']


class GridIterator:
    """Iterator for looping over a the grids of the time series"""
    def __init__(self, gridtimeseries):
        self.gridtimeseries = gridtimeseries
        self.index = 0

    def __iter__(self):
        return self

    def next(self):
        if self.index == len(self.gridtimeseries.get_times()):
            raise StopIteration
        self.index = self.index + 1
        return self.gridtimeseries.extract_grid(self.index-1)


class GridTimeSeries(AbstractFeature):
    """Class implementing a time series of grids"""

    def __init__(self,
                 identifier=None,
                 title=None,
                 description=None, source=None, metadata=None,
                 fields=None,
                 longitudes=None, latitudes=None, depths=None, times=None,
                 bbox=None,
                 projection=Projection()):
        """Constructor for GridTimeSeries

        Args:
            longitudes (tuple or numpy array or Field): Grid longitudes can be
                initialized:
                  * with a tuple (min, max, step) : only for a regular grid.
                    The longitudes are automatically generated.
                  * a numpy array of float values
                  * a Field (for instance copied from another grid object)

            latitudes (tuple or numpy array or Field): Grid latitudes can be
                initialized:
                  * with a tuple (min, max, step) : only for a regular grid.
                    The latitudes are automatically generated.
                  * a numpy array of float values
                  * a Field (for instance copied from another grid object)
        """
        if bbox:
            bbox = bbox
        else:
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
        super(GridTimeSeries, self).__init__(
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
                    bbox=bbox
                    )
        # projection parameters
        self.projection = projection
        return

    def __iter__(self):
        """Returns the iterator"""
        return GridIterator(self)

    def get_geolocation_dimnames(self):
        """
        Returns the geolocation dimension names defining a gridTimeSeries
        """
        return ['time', 'y', 'x']

    def get_geolocation_dimsizes(self):
        """
        Returns the geolocation dimension sizes of the gridTimeSeries object
        """
        #special trick to handle unregular grid
        lat = self.get_lat()
        lon = self.get_lon()
        if len(lat.shape)==1:
            return collections.OrderedDict([('time', len(self.get_times())),
                                            ('y', len(lat)),
                                            ('x', len(lon))
                                            ])
        elif len(lat.shape)==2:
            return collections.OrderedDict([('time', len(self.get_times())),
                                            ('y', lat.shape[0]),
                                            ('x', lon.shape[1])
                                            ])

    def get_geolocation_field_dimnames(self, fieldname):
        """
        Returns the dimensions of a geolocation field
        """
        if fieldname == 'depth':
            return ('depth',)
        elif fieldname == 'time':
            return ('time',)
        elif fieldname == 'lat':
                return ('y', 'x',)
        elif fieldname == 'lon':
                return ('x', 'y',)

    def time2t(self, time):
        """Get the indice of the closest value in the list of time steps"""
        numtime = date2num(time, self.get_time_units())
        return min(numpy.abs(self.get_times() - numtime).argmin(),
                   len(self.get_times()) - 1
                   )

    def add_grid(self, grid, step=None):
        """Add or update a time step to a grid time series.

        If the time step already exists, the corresponding values will be
        replaced.

        Args:
            step (int): specify directly the index whithin the time series
                where to insert the grid (replacing the current values). If
                None, the time of the grid is used to search for the
                corresponding index within the time series.
        """
        time = num2date(grid.get_times()[0], grid.get_time_units())
        logging.debug("Source grid time : %s", time)
        epsilon = 1
        if step is None:
            step = 0
            steptime = num2date(self.get_times()[step],
                                self.get_time_units()
                                )
            while abs(time - steptime).total_seconds() > epsilon:
                step += 1
        if step >= len(self.get_times()):
            raise Exception('Time not found')
        else:
            logging.debug('Index in time series : %d', step)
            # CASE 1 : time step already exists => update
            # check if variable (and corresponding record) already exists
            for fieldname in grid.get_fieldnames():
                logging.debug( 'fieldname %s',fieldname)
                data = grid.get_values(fieldname)
                if self.has_field(fieldname):
                    logging.debug( 'gridtimeserie field already exist')
                    logging.debug( "ADDING %s %s", fieldname, step)
                    field = self.get_field(fieldname)
                    logging.debug('field %s', field)
                    field._values[step, :, :] = data
                else:
                    logging.debug( '[gridtimeserie] field %s doesnt exist yet in ts',fieldname)
                    field_to_add = grid.get_field(fieldname)
                    self.add_field(field_to_add)
        return

    def extract_grid(self, index=None, time=None):
        """Extract a grid from a grid time series

        The grid to extract is defined either by its time step index or a time
        value. Both arguments are exclusive

        Args:
            index (int): index of the time step to extract from the grid time
                series
            time (:class:`datetime`): time of the grid step to extract from the
                grid time series

        Return:
            :class:`Grid`: extracted grid.
        """
        if time is not None and index is not None:
            raise Exception("index and time arguments are exclusive")
        elif index is None and time is None:
            raise Exception("No time step specified")
        if time is not None:
            # get closest time step in series to time
            times = self.get_datetimes()
            closest_time = min(times,
                               key=lambda val: abs(time - val))
            index = numpy.where(times == closest_time)[0]
        g = Grid()
        new_mapper = self.get_mapper().__class__(
            url=self.get_url(),
            **self.get_mapper().args)
        g.load(new_mapper, view={'time': slice(index, index + 1)})
        return g
    
    def is_unique_grid_time(self):
        """
        Return True if a unique time is associated with the grid
        (like in L4 products), False if there is time value per pixel
        (like in L3)
        """
        dims = self.get_geolocation_field('time').dimensions.values()[0]
        return (dims == 1)
    
    def save(self, output=None, attrs=None):
        """
        Save the gridtimeserie to a storage (file,...)

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
                mapper.create_dim('time', self.get_geolocation_field('time').get_dimsize('time'))
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
            if not 'title' in globalattr and self.title is not None:
                globalattr['title'] = self.title
            if not 'summary' in globalattr and self.description is not None:
                globalattr['summary'] = self.description
            if not 'id' in globalattr and self.identifier is not None:
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
            if not 'time_coverage_start' in globalattr:
                times = self.get_times()
                mintime = num2date(times.min(),
                                           self.get_time_units())
                globalattr['time_coverage_start'] \
                        = mintime.strftime('%Y%m%dT%H%M%S')
            if not 'time_coverage_end' in globalattr:
                times = self.get_times()
                maxtime = num2date(times.max(),
                                           self.get_time_units())
                globalattr['time_coverage_end'] \
                        = maxtime.strftime('%Y%m%dT%H%M%S')
            globalattr['cdm_data_type'] = 'gridtimeserie'
            mapper.write_global_attributes(globalattr)
            # creating records
            for geof in self._geolocation_fields:
                if self._geolocation_fields['time'] is None:
                    raise Exception('No time information defined')
                if self._geolocation_fields[geof] is None:
                    logging.warning('Missing geolocation variable : %s', geof)
                else:
                    mapper.create_field(self._geolocation_fields[geof],
                                        dim_translation
                                        )
            for dataf in self._fields.keys():
                if dataf not in self._geolocation_fields:
                    mapper.create_field(self._fields[dataf],
                                        dim_translation,
                                        feature='GridTimeSeries')
        # saving records
        for geof in self._geolocation_fields:
            field = self._geolocation_fields[geof]
            
            if field and not field.is_saved():
                mapper.write_field(self._geolocation_fields[geof])
        for dataf in self._fields.keys():
            if not self._fields[dataf].is_saved():
                mapper.write_field(self._fields[dataf])
        return
