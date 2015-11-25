# encoding: utf-8
"""
.. module::cerbere.datamodel.abstractfeature

Abstract class for all data feature objects

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
import logging
import datetime
import copy
import collections
from abc import ABCMeta, abstractmethod

import netCDF4
import numpy

from .variable import Variable
from .field import Field
from ..mapper.abstractmapper import AbstractMapper
from .. import DEFAULT_TIME_UNITS, CF_AUTHORITY

__all__ = ['AbstractFeature']

# Default convention for metadata and data model


class AbstractFeature(object):
    """
    Abstract class for all feature objects

    Args:
        fields (dict): a dictionary of :class:`Field` objects to be added
            to the feature content, where the key is a field identifier

        bbox (tuple): describes the limites of the covered geographical
            area as a tuple (lonmin,latmin,lonmax,latmax)

        times: must be either an array of datetime objects or a tuple
            (values, units) where values is an array of numbers and units
            is expressed in CF convention (ex: 'seconds since
            1981-01-01 00:00:00')
    """
    __metaclass__ = ABCMeta

    def __init__(self, identifier=None, title=None, description=None,
                 source=None, metadata=None, fields=None,
                 longitudes=None, latitudes=None, depths=None, times=None,
                 bbox=None,
                 spatial_resolution=None):
        """
        """
        object.__init__(self)
        if fields is not None:
            assert (type(fields) is dict), "fields must be a dictionary"
        self._datastorage = None
        # geolocation fields (time,lat,lon)
        self._geolocation_fields = {'time': None,
                                    'lat': None,
                                    'lon': None,
                                    'z': None
                                    }
        self.spatial_resolution = spatial_resolution
        if longitudes is not None:
            self.set_lon(longitudes)
        if latitudes is not None:
            self.set_lat(latitudes)
        if depths is not None:
            self.set_z(depths)
        self._bbox = None
        if bbox:
            self._bbox = bbox
        elif longitudes is not None and latitudes is not None:
            self._bbox = self.get_bbox()

        self._wkt_bbox = None
        if times is not None:
            if type(times) is list:
                times = numpy.array(times)
            if isinstance(times, numpy.ndarray) and\
                    isinstance(times[0], datetime.datetime):
                val = numpy.array(netCDF4.date2num(times,
                                                units=DEFAULT_TIME_UNITS))
                self.set_times(val,
                               units=DEFAULT_TIME_UNITS)
            elif isinstance(times, datetime.datetime):
                self.set_times(netCDF4.date2num([times],
                                                units=DEFAULT_TIME_UNITS),
                               units=DEFAULT_TIME_UNITS)
            elif isinstance(times, Field):
                self._geolocation_fields['time'] = times
            elif type(times) is tuple:
                if len(times) != 2:
                    raise Exception("If you provide times as numbers, you must"
                                    "also provide the units as a tuple "
                                    "(values(array), units(string))")
                self.set_times(times[0], units=times[1])
            else:
                raise Exception("Unknown type for times: ", type(times))
        self.identifier = identifier
        self.title = title
        self.description = description
        self.observation_source = source
        self.metadata = metadata
        if not metadata:
            self.metadata = {}
        self._fields = fields
        if fields is None:
            self._fields = {}
        return

    def __str__(self):
        result = 'FEATURE : %s\n' % self.__class__.__name__
        result = result + 'DIMENSIONS :\n'
        for dim in self.get_geolocation_dimnames():
            result = result + '   .\t%s\n' % dim
        result = result + 'GEOLOCATION :\n'
        for geof, val in self._geolocation_fields.items():
            if val:
                result = result + '   .\t%s\n' % geof
                result = result + '   \t\t+ dimensions :\t%s\n'\
                    % str(val.dimensions)
        result = result + 'DATA :\n'
        for field in self.get_fieldnames():
            result = result + '   .\t%s\n' % field
            result = result + '   \t\t+ dimensions :\t%s\n' \
                % str(self._fields[field].dimensions)
        return result

    @classmethod
    def get_model_name(cls):
        """Return the name of the datamodel"""
        return cls.__name__

    def get_bbox(self):
        '''
        returns the bounding box, e.g. the south/north and east/west
        limits of the feature.

        Return:
            tuple: the bounding box, always expressed as a tuple (lon min,
                lat_min, lon_max, lat_max)
        '''
        if self._bbox is None:
            lats = self.get_lat()
            lons = self.get_lon()
            self._bbox = (lons.min(), lats.min(), lons.max(), lats.max())
        return self._bbox

    def get_wkt_bbox(self):
        """Return the bounding box in WKT format."""
        #Generate a list of sensible names for the string format
        if self._wkt_bbox is None:
            cardinal_names = ['lon_min', 'lat_min', 'lon_max', 'lat_max']

            #Create empty WKT string
            polygon = 'POLYGON (({lon_min} {lat_min}, {lon_max} '\
                '{lat_min}, {lon_max} {lat_max}, {lon_min} {lat_max}, '\
                '{lon_min} {lat_min}))'

            #Generate a dictionary with these names and the bbox info
            bbox = dict(zip(cardinal_names, self.get_bbox()))

            #Return the formated string as WKT.
            self._wkt_bbox = polygon.format(**bbox)
        return self._wkt_bbox

    @abstractmethod
    def get_geolocation_dimnames(self):
        """Returns the geolocation dimensions defining the feature

        Returns:
            list<string>: the list of geolocation dimensions, using their
                standard name.

        Note:
            Abstract function to be overriden.
        """
        raise NotImplementedError

    @abstractmethod
    def get_geolocation_field_dimnames(self, fieldname):
        """Returns the dimensions of a geolocation field.

        Note:
            Abstract function to be overriden.
        """
        raise NotImplementedError

    @abstractmethod
    def get_geolocation_dimsizes(self):
        """Returns the geolocation dimensions defining the feature.

        Note:
            Abstract function to be overriden.
        """
        raise NotImplementedError

    def get_geolocation_dims(self, fieldname):
        """Return the dimensions of a geolocation field

        Return:
            dict: an ordered dictionary where key/values are the dimension
                names and sizes
        """
        dims = collections.OrderedDict()
        for dim in list(self.get_geolocation_field_dimnames(fieldname)):
            dims[dim] = self.get_geolocation_dimsizes()[dim]
        return dims

    def check_storage_matching(self):
        """Check that the data storage fits the feature.

        Verify the existence of required geolocation fields and dimensions.
        """
        dim_validity = self._datastorage.check_geolocation_dimensions(
            self.get_geolocation_dimnames()
            )
        if not dim_validity:
            logging.warning("The mapper's geolocation dimensions don't match "
                            "all the expected ones:")
            for dim in self.get_geolocation_dimnames():
                if not self._datastorage.check_geolocation_dimensions([dim]):
                    logging.warning("  => Missing geolocation dimension : %s",
                                    dim)
        geolocation_validity = self._datastorage.check_geolocation_fields(
                                        ['time', 'lat', 'lon']
                                        )
        if not geolocation_validity:
            logging.warning("The mapper's geolocation fields don't match all "
                            "the expected ones:")
            for f in ['time', 'lat', 'lon']:
                if not self._datastorage.check_geolocation_fields([f]):
                    logging.warning("  => Missing geolocation field : %s", f)
        return (dim_validity and geolocation_validity)

    def _get_field_values(self, field, slices=None, indices=None, cache=False,
                          **kwargs):
        """Return the data stored in a field object.

        Subsets can be returned by using either slices or indices arguments.
        slices and indices are exclusive and can not be both specified. Both
        must be provided as dictionaries where keys are the names of the
        dimensions to subset. Values are slice objects in the case of `slice`
        argument and numbers in the case of `indices` argument. Only the
        subsetted dimensions need to be provided (the full range s assumed for
        the other dimensions).

        Args:
            field (:class:Field): field object from which to read the data

            slices (dictionary): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            indices (list): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            cache (bool): if cache is True, the data read from file are kept in
                memory. The full field data array is kept in cache, slices or
                indices are ignored by caching (though the result returned by
                this function will be a subset if slices or indices are
                provided. Use the `view` concept in mappers instead if you want
                frequent accesses to a subset of a file.
                Default is False.

        Returns:
            numpy.ma.array: the requested data.
        """
        if not isinstance(field, Field):
            raise Exception("field must be a Field class instance.")
        return field.get_values(slices, indices, cache, **kwargs)

    def get_values(self, fieldname, slices=None, indices=None, cache=False,
                   **kwargs):
        """Return the data of a field, given its name.

        Subsets can be returned by using either slices or indices arguments.
        slices and indices are exclusive and can not be both specified. Both
        must be provided as dictionaries where keys are the names of the
        dimensions to subset. Values are slice objects in the case of `slice`
        argument and numbers in the case of `indices` argument. Only the
        subsetted dimensions need to be provided (the full range s assumed for
        the other dimensions).

        Args:
            fieldname (string): shortname of the field.

            slices (dictionary): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            indices (list): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            cache (bool): if cache is True, the data read from file are kept in
                memory. The full field data array is kept in cache, slices or
                indices are ignored by caching (though the result returned by
                this function will be a subset if slices or indices are
                provided. Use the `view` concept in mappers instead if you want
                frequent accesses to a subset of a file.
                Default is False.

        Returns:
            numpy.ma.array: the requested data.
        """
        field = self._fields[fieldname]
        return self._get_field_values(field, slices, indices, cache, **kwargs)

    def get_times(self, slices=None, indices=None, cache=True, **kwargs):
        """Return the times of a feature.

        The time values are returned as numbers. Use the
        :func:`get_time_units` function to convert these values to date or
        time objects.

        Subsets can be returned by using either slices or indices arguments.
        slices and indices are exclusive and can not be both specified. Both
        must be provided as dictionaries where keys are the names of the
        dimensions to subset. Values are slice objects in the case of `slice`
        argument and numbers in the case of `indices` argument. Only the
        subsetted dimensions need to be provided (the full range s assumed for
        the other dimensions).

        Args:
            slices (dictionary): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            indices (list): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            cache (bool): if cache is True, the data read from file are kept in
                memory. The full field data array is kept in cache, slices or
                indices are ignored by caching (though the result returned by
                this function will be a subset if slices or indices are
                provided. Use the `view` concept in mappers instead if you want
                frequent accesses to a subset of a file.
                Default is False.

        Returns:
            numpy.ma.array: the requested data.
        """
        field = self._geolocation_fields['time']
        return self._get_field_values(field, slices, indices, cache, **kwargs)

    def get_time_units(self):
        """Return the time units as a CF convention compliant string"""
        return self._geolocation_fields['time'].units

    def get_datetimes(self, slices=None, indices=None, cache=False, **kwargs):
        """Return the time values of a feature as datetime objects.

        Subsets can be returned by using either slices or indices arguments.
        slices and indices are exclusive and can not be both specified. Both
        must be provided as dictionaries where keys are the names of the
        dimensions to subset. Values are slice objects in the case of `slice`
        argument and numbers in the case of `indices` argument. Only the
        subsetted dimensions need to be provided (the full range s assumed for
        the other dimensions).

        Args:
            slices (dictionary): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            indices (list): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            cache (bool): if cache is True, the data read from file are kept in
                memory. The full field data array is kept in cache, slices or
                indices are ignored by caching (though the result returned by
                this function will be a subset if slices or indices are
                provided. Use the `view` concept in mappers instead if you want
                frequent accesses to a subset of a file.
                Default is False.

        Returns:
            numpy.ma.array<datetime>: the requested data.

        Warning:
            To be used with caution as conversion to datetime objects can take
            a long time. Use only when accessing single elements for instance.
        """
        field = self._geolocation_fields['time']
        values = self._get_field_values(field, slices, indices, cache, **kwargs)
        if not isinstance(values, datetime.datetime)\
                or not isinstance(values[0], datetime.datetime):
            values = netCDF4.num2date(numpy.array(values),
                                      field.units)
        return values

    def _get_time_dimensions(self, values):
        """Return the name and size of each dimension of the `time` geolocation
        field of the feature class.

        Internal function, only used when the :class:`Field` instance does not
        exist yet. `values` are used to discriminate between different cases
        for which the number of dimensions for `time` may vary (such as
        :class:`Grid` which can have time defined as a single value or a 2D
        array.

        Args:
            values (:class:`numpy.ma.MaskedArray`): time values, as a array of
                floats/integers or a Field.

        Return:
            OrderedDict<dimension name, dimension size>: dimensions of time
                field for the feature class
        """
        if isinstance(values, numpy.ndarray) and \
            isinstance(values[0], datetime.datetime):
            return collections.OrderedDict([('time', 1)])
        else:
            dims = collections.OrderedDict(
                    zip(list(self.get_geolocation_field_dimnames('time')),
                        list(values.shape))
                    )
        return dims

    def set_times(self, values, units=DEFAULT_TIME_UNITS):
        """
        Set the times

        Args:
            values (:class:`numpy.ma.MaskedArray` or :class:`Field`): time
                values, as a array of floats/integers or a Field object

            units (string): the time units as a CF convention compliant string.
                Only used when values is an array of floats/integers.
        """
        if not isinstance(values, Field):
            var = Variable(shortname='time',
                           description='time',
                           authority=CF_AUTHORITY,
                           standardname='time'
                           )
            dims = self._get_time_dimensions(values)
            field = Field(var,
                          dims,
                          values=values
                          )
            field.units = units
            self._geolocation_fields['time'] = field
        else:
            self._geolocation_fields['time'] = values
        return

    def set_datetimes(self, values, units=DEFAULT_TIME_UNITS):
        """Same as set_times method but providing a list of datetime as input.
        """
        return self.set_times(netCDF4.date2num(values, units))

    def get_lon(self, slices=None, indices=None, cache=True, **kwargs):
        """Return the longitude values of a feature.

        Subsets can be returned by using either slices or indices arguments.
        slices and indices are exclusive and can not be both specified. Both
        must be provided as dictionaries where keys are the names of the
        dimensions to subset. Values are slice objects in the case of `slice`
        argument and numbers in the case of `indices` argument. Only the
        subsetted dimensions need to be provided (the full range s assumed for
        the other dimensions).

        Args:
            slices (dictionary): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            indices (list): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            cache (bool): if cache is True, the data read from file are kept in
                memory. The full field data array is kept in cache, slices or
                indices are ignored by caching (though the result returned by
                this function will be a subset if slices or indices are
                provided. Use the `view` concept in mappers instead if you want
                frequent accesses to a subset of a file.
                Default is False.

        Returns:
            numpy.ma.array<float>: the requested longitudes.
        """
        field = self._geolocation_fields['lon']
        return self._get_field_values(field, slices, indices, cache, **kwargs)

    def set_lon(self, values):
        """Set the longitudes

        Args:
            values (:class:`numpy.ma.MaskedArray` or :class:`Field`): longitude
                values, as a array of floats or a Field object
        """
        if not isinstance(values, Field):
            var = Variable(
                        shortname='lon',
                        description='longitude',
                        authority=CF_AUTHORITY,
                        standardname='longitude'
                        )
            dims = collections.OrderedDict(
                    zip(list(self.get_geolocation_field_dimnames('lon')),
                        list(values.shape)
                        )
                    )
            field = Field(
                            var,
                            dims,
                            values=values
                            )
            field.units = 'degrees_east'
            self._geolocation_fields['lon'] = field
        else:
            self._geolocation_fields['lon'] = values
        return

    def get_lat(self, slices=None, indices=None, cache=True, **kwargs):
        """Return the latitude values of a feature.

        Subsets can be returned by using either slices or indices arguments.
        slices and indices are exclusive and can not be both specified. Both
        must be provided as dictionaries where keys are the names of the
        dimensions to subset. Values are slice objects in the case of `slice`
        argument and numbers in the case of `indices` argument. Only the
        subsetted dimensions need to be provided (the full range s assumed for
        the other dimensions).

        Args:
            slices (dictionary): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            indices (list): a dictionary where keys are the dimensions
                to slice and values are slice objects.

            cache (bool): if cache is True, the data read from file are kept in
                memory. The full field data array is kept in cache, slices or
                indices are ignored by caching (though the result returned by
                this function will be a subset if slices or indices are
                provided. Use the `view` concept in mappers instead if you want
                frequent accesses to a subset of a file.
                Default is False.

        Returns:
            numpy.ma.array<float>: the requested latitudes.
        """
        field = self._geolocation_fields['lat']
        return self._get_field_values(field, slices, indices, cache, **kwargs)

    def set_lat(self, values):
        """Set the latitudes

        Args:
            values (:class:`numpy.ndarray` or :class:`Field`): longitude
                values, as a array of floats or a Field object.
        """
        if not isinstance(values, Field):
            var = Variable(
                        shortname='lat', description='latitude',
                        authority=CF_AUTHORITY,
                        standardname='latitude'
                        )
            dims = collections.OrderedDict(
                        zip(list(self.get_geolocation_field_dimnames('lat')),
                            list(values.shape))
                        )
            field = Field(
                            var,
                            dims,
                            values=values
                            )
            field.units = 'degrees_north'
            self._geolocation_fields['lat'] = field
        else:
            self._geolocation_fields['lat'] = values
        return

    def get_z(self, slices=None, indices=None, cache=True, **kwargs):
        '''
        if cache is True, data read from storage are kept in memory
        '''
        field = self._geolocation_fields['z']
        return self._get_field_values(field, slices, indices, cache, **kwargs)

    def set_z(self, values, ztype='depth'):
        """Set the depth(s)

        Args:
            values (:class:`numpy.ma.MaskedArray` or :class:`Field`): depth
                values, as a array of floats or a Field object

            ztype (string): type is depth (positive down) or height (positive
                up)
        """
        if not isinstance(values, Field):
            var = Variable(
                        shortname=ztype,
                        description=ztype,
                        authority=CF_AUTHORITY,
                        standardname=ztype
                        )
            dims = collections.OrderedDict(
                        zip(list(self.get_geolocation_field_dimnames('depth')),
                            list(values.shape))
                        )
            field = Field(
                        var,
                        dims,
                        values=values,
                        attributes={'depth':{'positive':'down'},
                                    'height':{'positive':'up'}}[ztype]
                        )
            field.units = 'm'
            self._geolocation_fields['z'] = field
        else:
            self._geolocation_fields['z'] = values

    def load(self, storage_mapper, readonly=True, view=None):
        """Restore the feature from a file (or any type of storage.)

        It does not read the actual data at this point for memory issues,
        only the metadata and geolocation information necessary to initialize
        the feature model. Physical read access to the data from the storage
        is performed transparently when attempting to get the values of a
        field.

        Args:
            storage_mapper (:class:`AbstractMapper`): an instance of mapper
                class, corresponding to the storage of the feature

            view (dict, optional): a dictionary where keys are dimension names
                and values are slices. A view can be set on a file, meaning
                that only the subset defined by this view will be accessible.
                This view is expressed as any subset (see :func:`get_values`).
                For example:

                .. code-block:: python

                   view = {'time':slice(0,0), 'lat':slice(200,300),
                      'lon':slice(200,300)}
        """
        if self._datastorage:
            raise Exception("A mapper is already attached to the feature")
        if isinstance(storage_mapper, AbstractMapper):
            self._datastorage = storage_mapper
            # open mapper object in read_only if not already opened
            logging.debug("opening %s", storage_mapper.get_url())
            if not self._datastorage.is_opened():
                if readonly:
                    dims = self.get_geolocation_dimnames()
                    self._datastorage.open(
                        datamodel_geolocation_dims=dims,
                        view=view,
                        datamodel=self.__class__.get_model_name()
                        )
                else:
                    raise Exception("Update mode not yet implemented")
            # check if mapper contains correct geolocation/dimension
            # information for the current data model
            #
            # check_storage_matching() is a virtual method to be implemented
            # by each derived data model class
            if not self.check_storage_matching():
                raise Exception("mapper structure does not match the feature")
            # Get the dimensions defining the data model
            # the dimensions (defining the data arrays) depend on the feature
            # pattern :
            # timeseries = (time)
            # trajectory = (time)
            # grid = (x, y) or (lon, lat)
            # gridtimeseries = (time, lon, lat) or (time, x, y)
            # swath = (row, cell)
            # ...
            dims = self.get_geolocation_dimnames()
            logging.debug("dimensions: %s ", dims)
            # get geolocation fields
            for geof in ['time', 'lat', 'lon', 'z']:
                try:
                    rec = storage_mapper.read_field(geof)
                    self._geolocation_fields[geof] = rec
                except:
                    logging.debug('geolocation field %s is not available'
                                  % geof)
            # get the list of data fields
            variables = self._datastorage.get_fieldnames()
            logging.debug("fields : ")
            logging.debug(variables)
            # initialize corresponding Field for each variable field
            for var in variables:
                rec = storage_mapper.read_field(var)
                self.add_field(rec, new=False)
            # initialize bbox
            self._bbox = storage_mapper.get_bbox()
            # attributes
            self.get_metadata()
            self.get_start_time()
            self.get_end_time()
            return True
        else:
            raise Exception('storage_mapper must be a mapper class')

    def add_field(self, field, new=True):
        """Add a field to the feature.

        Args:
            field (Field): the field is described by a :class:`Field`
                object, containing the actual array of data values (or only a
                storage descriptor when the data have not yet been read from
                the files) and the description of the observed quantity.
            new (boolean): True if is a new field (default), False if it is
                already saved on disk. For internal usage only as from user
                perspective you will add new fields only to a feature.
        """
        if field.variable is None:
            raise Exception('No variable associated to the record')
        if field.get_name() in self._fields:
            raise Exception("Field already existing in feature. Can not add :",
                            field.get_name())
        self._fields[field.get_name()] = field
        # field is not yet saved
        if field.handler is not None and new:
            field.handler.set_unsaved()
        return

    def has_field(self, fieldname):
        """Return True if the field exists in the feature"""
        return fieldname in self._fields.keys()

    def get_fieldnames(self):
        """Return the names of the data fields stored in the feature"""
        return self._fields.keys()

    def get_geolocation_fields(self):
        """
        Return the names of the geolocation fields
        """
        return self._geolocation_fields.keys()

    def get_field(self, fieldname):
        '''
        Return a field from the feature
        '''
        if fieldname not in self._fields:
            raise Exception("Field is not exising: ", fieldname)
        return self._fields[fieldname]

    def get_geolocation_field(self, fieldname):
        '''
        Return a geolocation field from the feature
        '''
        return self._geolocation_fields[fieldname]

    def get_metadata(self):
        '''
        returns the global attributes (metadata) of the feature
        '''
        if self._datastorage is None:
            return self.metadata
        else:
            metadata = {}
            for att in self._datastorage.read_global_attributes():
                metadata[att] = self._datastorage.read_global_attribute(att)
            return metadata

    def get_start_time(self):
        """Return start time of feature's temporal coverage.

        Returns:
            datetime: start time of temporal coverage. None if no valid time
                was found.
        """
        if self._datastorage is None:
            metadata = self.get_metadata()
            if 'time_coverage_start' in metadata:
                if isinstance(metadata['time_coverage_start'], datetime):
                    return metadata['time_coverage_start']
                else:
                    logging.warning('time_coverage_start is not a datetime'
                                    ' object')
            times = self.get_times()
            if times.size == 0:
                return None
            mintime = netCDF4.num2date(times.min(),
                                       self.get_time_units())
            return mintime
        else:
            return self._datastorage.get_start_time()

    def get_end_time(self):
        """Return end time of feature's temporal coverage.

        Returns:
            datetime: end time of temporal coverage. None if no valid time was
                found.
        """
        if self._datastorage is None:
            metadata = self.get_metadata()
            if 'time_coverage_end' in metadata:
                if isinstance(metadata['time_coverage_end'], datetime):
                    return metadata['time_coverage_end']
                else:
                    logging.warning('time_coverage_end is not a datetime'
                                    ' object')
            times = self.get_times()
            if times.size == 0:
                return None
            maxtime = netCDF4.num2date(times.max(),
                                       self.get_time_units())
            return maxtime
        else:
            return self._datastorage.get_end_time()

    def extract_subset(self, boundaries=None, slices=None, fields=None):
        """Extract a subset feature from the feature.

        The created subset is a new feature object of the same class without
        any reference to the source.

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
        raise NotImplementedError

    def extract_field(
            self, fieldname, slices=None, indices=None, padding=False):
        """
        Create a copy of an field, limiting to a set of slices or indices, and
        padding out as required.

        Args
        :param fieldname: The name of the field to extract
        :param slices: slices of the geolocation dimensions, in case
        sub-setting is required
        :param indices: indices to pe passed to the extraction method.
        :param padding: True to pad out feature with fill values to the extent
        of the dimensions.
        """
        if fieldname in self._geolocation_fields:
            field = self._geolocation_fields[fieldname]
        else:
            field = self._fields[fieldname]
        if field is None:
            return None
        new_dims = copy.copy(field.dimensions)
        if slices:
            for dim in new_dims.keys():
                if dim in slices:
                    #start,end,step = slices[d]
                    if not slices[dim].start:
                        start = 0
                    else:
                        start = slices[dim].start
                    if not slices[dim].stop:
                        stop = field.dimensions[dim]
                    else:
                        stop = slices[dim].stop
                    new_dims[dim] = stop - start
        elif indices:
            for dim in new_dims:
                if dim in indices:
                    new_dims[dim] = len(indices[dim])
                    
        new_field = Field(
            copy.copy(field.variable),
            new_dims,
            fields=copy.copy(field.components),
            datatype=field.datatype
        )
        new_field.set_metadata(field.get_metadata())
        new_field.set_values(
            field.get_values(slices=slices, indices=indices, padding=padding)
        )
        return new_field

    def get_url(self):
        '''
        Return the URL of the storage resource where feature's data are saved.
        '''
        if self._datastorage is None:
            return None
        else:
            return self._datastorage.get_url()

    def get_mapper(self):
        '''
        Return the handler of the storage resource where feature's data are
         saved.
        '''
        return self._datastorage
