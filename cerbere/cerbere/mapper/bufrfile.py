# -*- coding: utf-8 -*-
"""
.. module::cerbere.mapper.bufrfile

Mapper classs for BUFR format

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
import collections
import datetime

import bufr
import numpy
import netCDF4
import logging
from cerbere.mapper import abstractmapper
from cerbere.datamodel import field
from cerbere.datamodel import variable
import cerbere.mapper.slices

UNITS = {'DEGREE TRUE': 'degree_true',
         'SECOND': 's',
         'M/S': 'm s-1',
         'M': 'm',
         'KG/(M**2)S': 'kg m-2 s',
         'K': 'K',
         'dB': 'dB',
         'NUMERIC': None,
         'DEGREE': 'degree'} 

REFERENCE_TIME = datetime.datetime(1991, 1, 1, 0, 0, 0, 0)


class BUFRFile(abstractmapper.AbstractMapper):
    """Mapper class to read BUFR files.

    It requires the python-bufr package available at
    http://www.pytroll.org

    Requires the BUFR tables from emos package in the `BUFR_TABLES` environment
    variable, e.g.

    .. code block: bash

        export BUFR_TABLES=/opt/lib/emos/bufrtables/
    """

    ATTRIBUTE_ENTRIES = {}
    SKIPPED_ENTRIES = ['year', 'month', 'day', 'hour', 'minute', 'second']

    FILLVALUE = 1.7e38

    def __init__(self, url=None, mode=abstractmapper.READ_ONLY, **kwargs):
        """
        Initialize a BUFR file mapper
        """
        super(BUFRFile, self).__init__(url=url, mode=mode, **kwargs)
        return

    def get_geolocation_field(self, fieldname):
        """Return the equivalent field name in the file format for a standard
        geolocation field (lat, lon, time, z).

        Used for internal purpose and should not be called directly.

        Args:
            fieldname (str): name of the standard geolocation field (lat, lon
                or time)

        Return:
            str: name of the corresponding field in the native file format.
                Returns None if no matching is found
        """
        return fieldname

    def get_matching_dimname(self, dimname):
        """Return the equivalent name in the native format for a standard
        dimension.

        This is a translation of the standard names to native ones. As BUFR
        does not have explicit dimensions, we return the standard dimension
        names, which are :

        * row, cell, time for :class:`~cerbere.datamodel.swath.Swath` or
          :class:`~cerbere.datamodel.image.Image`

        Args:
            dimname (str): standard dimension name.

        Returns:
            str: return the native name for the dimension. Return `dimname` if
                the input dimension has no standard name.

        See Also:
            see :func:`get_standard_dimname` for the reverse operation
        """
        return dimname

    def get_standard_dimname(self, dimname):
        """
        Returns the equivalent standard dimension name for a
        dimension in the native format.

        This is a translation of the native names to standard ones. As BUFR
        does not have explicit dimensions, we return the standard dimension
        names, which are :

        Args:
            dimname (string): native dimension name

        Return:
            str: the (translated) standard name for the dimension. Return
            `dimname` if the input dimension has no standard name.

        See Also:
            see :func:`get_matching_dimname` for the reverse operation
        """
        return dimname

    def open(self, view=None,datamodel_geolocation_dims=None,datamodel=None,):
        """Open the file (or any other type of storage)

        Args:
            view (dict, optional): a dictionary where keys are dimension names
                and values are slices. A view can be set on a file, meaning
                that only the subset defined by this view will be accessible.
                Any later use of slices in :func:`read_values` will be relative
                to the view defined here.

                This view is expressed as any subset (see :func:`read_values`).
                For example:

                .. code-block: python

                    view = {'time':slice(0,0), 'row':slice(200,300),
                        'cell':slice(200,300)}

        """
        self._handler = bufr.BUFRFile(self._url)
        super(BUFRFile, self).open(view=view,
                                   datamodel='Swath')
        if datamodel_geolocation_dims:
            self.datamodel_geolocation_dims = datamodel_geolocation_dims
        if datamodel is not None:
            self._feature_type = datamodel
        # force data buffering to get full file structure and content
        fieldnames = []
        units = []
        self._fields = {}
        data = []
        geolocdata = {}
        # read file content
        first = True
        nbrows = 0
        yearidx = None
        monthidx = None
        dayidx = None
        houridx = None
        minuteidx = None
        secondidx = None
        latidx = None
        lonidx = None
        indexes = {}
        duplicated_field_indice = {}
        for record in self._handler:
            nbrows += 1
            if first:
                nbcells = len(record[1].data)
                for entry in record:
                    # create valid name
                    entryname = (entry.name.lower().replace(' ','_').
                                 strip('_').strip('*'))
                    logging.debug('entry name %s',entryname)
                    # identify geolocation information. We assume the first met
                    # entries with geolocation names are the correct ones. The
                    # next ones with the same names are processed as additional
                    # variables
                    if entryname == 'year' and yearidx is None:
                        yearidx = entry.index
                        logging.debug('entry index year = %s',yearidx)
                    elif entryname == 'month' and monthidx is None:
                        monthidx = entry.index
                    elif entryname == 'day' and dayidx is None:
                        dayidx = entry.index
                    elif entryname == 'hour' and houridx is None:
                        houridx = entry.index
                    elif entryname == 'minute' and minuteidx is None:
                        minuteidx = entry.index
                    elif entryname == 'second' and secondidx is None:
                        secondidx = entry.index
                    elif ('latitude' in entryname and
                            latidx is None):
                        latidx = entry.index
                        geolocdata['lat'] = [entry.data]
                    elif ('longitude' in  entryname and
                            lonidx is None):
                        geolocdata['lon'] = [entry.data]
                        lonidx = entry.index
                    # decides if entry should be a field or an attribute
                    elif (entryname not in self.ATTRIBUTE_ENTRIES and
                            entryname not in self.SKIPPED_ENTRIES):
                        if entryname in duplicated_field_indice:
                            # case where different fields with the same name
                            # we had an incremental number and
                            # change the name of the first occurence
                            duplicated_field_indice[entryname] = 0
                            prev_occurence = fieldnames.index(entryname)
                            fieldnames[prev_occurence] = '_'.join([
                                entryname,
                                '%d' % duplicated_field_indice[entryname]
                                ])
                        if entryname in duplicated_field_indice:
                            duplicated_field_indice[entryname] += 1
                            entryname = '_'.join([
                                entryname,
                                '%d' % duplicated_field_indice[entryname]
                                ])
                        fieldnames.append(entryname)
                        indexes[entryname] = entry.index
                        # get unit in standard form
                        unit = entry.unit.strip()
                        if unit in UNITS:
                            units.append(UNITS[unit])
                        else:
                            units.append(unit)
                        data.append([entry.data])
                    elif entry.name in self.ATTRIBUTE_ENTRIES:
                        self.attributes[entryname] = entry.data[0]
                first = False
                for i in range(nbcells):
                    logging.debug('year %s', int(record[yearidx].data[i]))
                geolocdata['time'] = [numpy.array(
                    [(datetime.datetime(int(record[yearidx].data[i]),
                                        int(record[monthidx].data[i]),
                                        int(record[dayidx].data[i]),
                                        int(record[houridx].data[i]),
                                        int(record[minuteidx].data[i]),
                                        int(record[secondidx].data[i])) -
                        REFERENCE_TIME).total_seconds()
                     for i in range(nbcells)
                     ])]
            else:
                # read arrays of data
                current_reccord_range_length = record[secondidx].data.shape[0]
                if nbcells == current_reccord_range_length:
                    for i, fieldname in enumerate(fieldnames):
                        data[i].append(record[indexes[fieldname]].data)
                    geolocdata['lon'].append(record[lonidx].data)
                    geolocdata['lat'].append(record[latidx].data)

                    geolocdata['time'].append(numpy.array(
                        [(datetime.datetime(int(record[yearidx].data[i]),
                                            int(record[monthidx].data[i]),
                                            int(record[dayidx].data[i]),
                                            int(record[houridx].data[i]),
                                            int(record[minuteidx].data[i]),
                                            int(record[secondidx].data[i]))
                            - REFERENCE_TIME).total_seconds()
                         for i in range(nbcells)
                         ]))
        del self._handler
        self._handler = self
        # get dimensions (take into account the view)
        self._dimensions = collections.OrderedDict([
            ('row', nbrows),
            ('cell', nbcells)
            ])
        newslices = cerbere.mapper.slices.get_absolute_slices(
            self.view,
            slices=None,
            dimnames=self._dimensions.keys(),
            dimsizes=self._dimensions.values())
        # get fields and cache data
        self._fields = {}
        for i, fieldname in enumerate(fieldnames):
            varobj = variable.Variable(fieldname, fieldname.replace('_', ' '))
            newfield = field.Field(
                variable=varobj,
                dimensions=self._dimensions,
                datatype=numpy.float32,
                values=numpy.ma.masked_equal(numpy.vstack(data[i]),
                                             self.FILLVALUE)[tuple(newslices)],
                fillvalue=self.FILLVALUE,
                units=units[i],
                )
            self._fields[fieldname] = newfield
        # geolocation
        self._geofields = {}
        varobj = variable.Variable(shortname='lat', description='latitude')
        newfield = field.Field(
            variable=varobj,
            dimensions=self._dimensions,
            datatype=numpy.float32,
            values=numpy.ma.masked_equal(numpy.vstack(geolocdata['lat']),
                                         self.FILLVALUE)[tuple(newslices)],
            fillvalue=self.FILLVALUE,
            units='degrees_north',
            )
        self._geofields['lat'] = newfield
        varobj = variable.Variable(shortname='lon', description='longitude')
        newfield = field.Field(
            variable=varobj,
            dimensions=self._dimensions,
            datatype=numpy.float32,
            values=numpy.ma.masked_equal(numpy.vstack(geolocdata['lon']),
                                         self.FILLVALUE)[tuple(newslices)],
            fillvalue=self.FILLVALUE,
            units='degrees_east',
            )
        self._geofields['lon'] = newfield
        var = variable.Variable(shortname='time', description='time')
        newfield = field.Field(
            variable=var,
            dimensions=self._dimensions,
            datatype=numpy.int32,
            values=numpy.vstack(geolocdata['time'])[tuple(newslices)],
            units=('seconds since %s'
                   % datetime.datetime.strftime(REFERENCE_TIME,
                                                "%Y-%m-%d %H:%M:%S"))
            )
        self._geofields['time'] = newfield
        return self._handler

    def close(self):
        """Close handler on storage"""
        pass

    def read_values(self, fieldname, slices=None):
        """Read the data of a field.

        Args:
            fieldname (str): name of the field which to read the data from

            slices (list of slice, optional): list of slices for the field if
                subsetting is requested. A slice must then be provided for each
                field dimension. The slices are relative to the opened view
                (see :func:open) if a view was set when opening the file.

        Return:
            MaskedArray: array of data read. Array type is the same as the
                storage type.
        """
        if fieldname in ['lat', 'lon', 'time']:
            field = self._geofields[fieldname]
        else:
            field = self._fields[fieldname]
        fulldims = self.get_full_dimensions(fieldname)
        newslices = cerbere.mapper.slices.get_absolute_slices(
            self.view,
            slices=slices,
            dimnames=fulldims.keys(),
            dimsizes=fulldims.values())
        return field._values[tuple(newslices)]

    def read_field(self, fieldname):
        """
        Return the :class:`cerbere.field.Field` object corresponding to
        the requested fieldname.

        The :class:`cerbere.field.Field` class contains all the metadata
        describing a field (equivalent to a variable in netCDF).

        Args:
            fieldname (str): name of the field

        Returns:
            :class:`cerbere.field.Field`: the corresponding field object
        """
        if fieldname in ['lat', 'lon', 'time']:
            return self._geofields[fieldname]
        else:
            return self._fields[fieldname]

    def get_fieldnames(self):
        """Returns the list of geophysical fields stored for the feature.

        The geolocation field names are excluded from this list.

        Returns:
            list<string>: list of field names
        """
        return self._fields.keys()

    def get_dimensions(self, fieldname=None):
        """Return the dimension's standard names of a file or a field in the
        file.

        Args:
            fieldname (str): the name of the field from which to get the
                dimensions. For a geolocation field, use the cerbere standard
                name (time, lat, lon), though native field name will work too.

        Returns:
            tuple<str>: the standard dimensions of the field or file.
        """
        return tuple(self._dimensions.keys())

    def get_dimsize(self, dimname):
        """Return the size of a dimension.

        Args:
            dimname (str): name of the dimension.

        Returns:
            int: size of the dimension.
        """
        return cerbere.mapper.slices.adjust_dimsize(
            self.view, dimname,
            self._dimensions[dimname])

    def read_field_attributes(self, fieldname):
        """Return the specific attributes of a field.

        Args:
            fieldname (str): name of the field.

        Returns:
            dict<string, string or number or datetime>: a dictionary where keys
                are the attribute names.
        """
        return []

    def read_global_attributes(self):
        """Returns the names of the global attributes.

        Returns:
            list<str>: the list of the attribute names.
        """
        return []

    def read_global_attribute(self, name):
        """Returns the value of a global attribute.

        Args:
            name (str): name of the global attribute.

        Returns:
            str, number or datetime: value of the corresponding attribute.
        """
        raise Exception("This attribute does not exist")

    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage.

        Returns:
            datetime: start time of the data in file.
        """
        return netCDF4.num2date(
            self._geofields['time']._values.min(),
            self._geofields['time'].units)

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage.

        Returns:
            datetime: end time of the data in file.
        """
        return netCDF4.num2date(
            self._geofields['time']._values.max(),
            self._geofields['time'].units)

    def write_field(self, fieldname):
        """Writes the field data on disk.

        Args:
            fieldname (str): name of the field to write.
        """
        raise NotImplementedError

    def create_field(self, field, dim_translation=None):
        """Creates a new field in the mapper.

        Creates the field structure but don't write yet its values array.

        Args:
            field (Field): the field to be created.

        See also:
            :func:`write_field` for writing the values array.
        """
        raise NotImplementedError

    def create_dim(self, dimname, size=None):
        """Add a new dimension.

        Args:
            dimname (str): name of the dimension.
            size (int): size of the dimension (unlimited if None)
        """
        raise NotImplementedError

    def write_global_attributes(self, attrs):
        """Write the global attributes of the file.

        Args:
            attrs (dict<string, string or number or datetime>): a dictionary
                containing the attributes names and values to be written.
        """
        raise NotImplementedError

    def get_bbox(self):
        """Returns the bounding box of the feature, as a tuple.

        Returns:
            tuple: bbox expressed as (lonmin, latmin, lonmax, latmax)
        """
        return (self._geofields['lon']._values.min(),
                self._geofields['lat']._values.min(),
                self._geofields['lon']._values.max(),
                self._geofields['lat']._values.max())

    def read_fillvalue(self, fieldname):
        """Read the fill value of a field.

        Args:
            fieldname (str): name of the field.

        Returns:
            number or char or str: fill value of the field. The type is the
                as the type of the data in the field.
        """
        return self.FILLVALUE
