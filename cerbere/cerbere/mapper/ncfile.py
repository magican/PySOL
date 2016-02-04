# -*- coding: utf-8 -*-
"""
.. module::cerbere.mapper.ncfile

Class to read netCDF files using CF conventions.

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
import os
import logging
import datetime
import copy
import collections
from dateutil import parser
#import pdb

import numpy
import netCDF4

from .. import READ_ONLY, WRITE_NEW, DEFAULT_TIME_UNITS
from .abstractmapper import AbstractMapper, FieldHandler, CorruptFileException
from ..datamodel.field import Field
from ..datamodel.variable import Variable


class NCFile(AbstractMapper):
    """Initialize a CF compliant NetCDF file mapper

    The file is not yet open (read mode) nor created (write mode), this
    is done explicitely by the caller when calling open()

    Args:
        ncformat (string): specifies the required NetCDF format version:
            NETCDF3_CLASSIC, NETCDF4, NETCDF4_CLASSIC,....

        geodim_matching (dict): explicitly provides the matching between
            standard names for the geolocation dimensions (keys)and the
            corresponding native names (values) in the file. If not
            provided, the library will try to guess with respect to known
            conventions.

        geofield_matching (dict): explicitly provides the matching between
            standard names for the geolocation fields (keys)and the
            corresponding native names (values) in the file. If not
            provided, the library will try to guess with respect to known
            conventions.

        is_reference (bool): if set to True, then creates a virtual time
            field if missing otherwise any attempt to load the file would
            fail. Used for instance with climatology files with no time
            value provided. The time will be set to 0-1-1 by default
            (COARDS convention).

        fillvalue (float): additional fillvalue to use to mask some
            invalid data. Some fields don't have a correct fill value set,
            and this value must be explicitly provided.

        center_on_greenwhich (bool, optional): if True, shift a grid so
            that it is centered on meridian 0 instead of meridian 180.
            Makes only sense for global grids in cylindrical projection.
            Default is False.
    """
    CONVENTIONS = 'CF-1.6, Unidata Observation Dataset v1.0'

    STANDARD_ATTRIBUTES_VALUES = collections.OrderedDict([
        ('id', ''),
        ('naming_authority', 'fr.ifremer.cersat'),
        ('Metadata_Conventions', 'Unidata Dataset Discovery v1.0'),
        ('standard_name_vocabulary',
         'NetCDF Climate and Forecast (CF) Metadata Convention'),
        ('institution', 'Institut Francais de Recherche et d\'Exploitation de'
            ' la Mer/Centre de Recherche et d\'Exploitation satellitaire'),
        ('institution_abbreviation', 'ifremer/cersat'),
        ('title', ''),
        ('summary', ''),
        ('cdm_data_type', ''),
        ('keywords', ''),
        ('keywords_vocabulary', 'NASA Global Change Master Directory (GCMD)'
            ' Science Keywords'),
        ('standard_name_vocabulary', 'NetCDF Climate and Forecast (CF) '
            'Metadata Convention'),
        ('scientific_project', ''),
        ('acknowledgement', ''),
        ('license', ''),
        ('format_version', ''),
        ('processing_software', 'Cersat/Cerbere 1.0'),
        ('product_version', ''),
        ('uuid', ''),
        ('processing_level', ''),
        ('history', ''),
        ('publisher_name', 'ifremer/cersat'),
        ('publisher_url', 'http,//cersat.ifremer.fr'),
        ('publisher_email', 'cersat@ifremer.fr'),
        ('creator_name', ''),
        ('creator_url', ''),
        ('creator_email', ''),
        ('references', ''),
        ('metadata_link', ''),
        ('source', ''),
        ('source_version', ''),
        ('platform', ''),
        ('platform_type', ''),
        ('sensor', ''),
        ('sensor_type', ''),
        ('band', ''),
        ('spatial_resolution', ''),
        ('geospatial_lat_min', ''),
        ('geospatial_lat_max', ''),
        ('geospatial_lat_units', 'degrees'),
        ('geospatial_lat_resolution', ''),
        ('geospatial_lon_min', ''),
        ('geospatial_lon_max', ''),
        ('geospatial_lon_units', 'degrees'),
        ('geospatial_lon_resolution', ''),
        ('geospatial_vertical_min', ''),
        ('geospatial_vertical_max', ''),
        ('geospatial_vertical_units', 'meters above mean sea level'),
        ('geospatial_vertical_positive', 'up'),
        ('time_coverage_start', ''),
        ('time_coverage_end', ''),
        ('time_coverage_resolution', '')
        ])

    TIME_FORMAT = '%Y%m%dT%H%M%SZ'

    def __init__(self, url=None,
                 mode=READ_ONLY,
                 ncformat='NETCDF4',
                 geodim_matching=None,
                 geofield_matching=None,
                 is_reference=False,
                 fillvalue=None,
                 center_on_greenwhich=False,
                 **kwargs):
        """
        """
        # save additional arguments before calling parent constructor
        kwargs['ncformat'] = ncformat
        kwargs['geodim_matching'] = geodim_matching
        kwargs['geofield_matching'] = geofield_matching
        kwargs['is_reference'] = is_reference
        kwargs['fillvalue'] = fillvalue
        kwargs['center_on_greenwhich'] = center_on_greenwhich
        super(NCFile, self).__init__(url=url, mode=mode, **kwargs)
        # Type of NetCDF format : NETCDF3_CLASSIC, NETCDF_CLASSIC,....
        self.format = ncformat
        # case of reference files (climatology etc) with no time value
        # attached
        self.is_reference = is_reference
        # additional fillvalue
        self.fillvalue = fillvalue
        # center on meridian 0
        self.center_on_greenwhich = center_on_greenwhich
        # initialize matching tables for geolocation information
        self.geodim_std2native = {
            'time': ['time'],
            'x': ['lon', 'longitude'],
            'y': ['lat', 'latitude'],
            'z': ['depth', 'height'],
            'station': ['mes', 'station'],
            'cell': ['ni', 'cell', 'col', 'ra_size', 'columns', 'NUMCELLS',
                     'across_track'],
            'row': ['nj', 'row', 'az_size', 'rows', 'NUMROWS',
                    'along_track'],
            }
        self.geodim_native2std = {
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
            'col': 'cell',
            'nj': 'row',
            'row': 'row',
            'az_size': 'row',
            'depth': 'z',
            'height': 'z',
            'rows': 'row',
            'columns': 'cell',
            'NUMROWS': 'row',
            'NUMCELLS': 'cell',
            'across_track': 'cell',
            'along_track': 'row'
            }
        if geodim_matching:
            for std, native in geodim_matching.items():
                self.geodim_std2native[std] = [native]
                self.geodim_native2std[native] = std
        self.geofield_std2native = {
            'time': ['time'],
            'lon': ['lon', 'longitude'],
            'lat': ['lat', 'latitude'],
            'z': ['depth', 'height'],
            }
        if geofield_matching:
            for std, native in geofield_matching.items():
                self.geofield_std2native[std] = [native]
        return

    def get_geolocation_field(self, fieldname):
        """Return the equivalent field name in the file format for a standard
        geolocation field (lat, lon, time, z).

        Used for internal purpose and should not be called directly.

        Args:
            fieldname (str): name of the standard geolocation field (lat, lon
                or time).

        Return:
            str: name of the corresponding field in the native file format.
                Returns None if no matching is found.
        """
        keys = self.get_handler().variables.keys()
        # case where file is a reference (possibly no time value defined)
        # => creates a dummy time
        if self.is_reference:
            keys.append('time')
        if fieldname in self.geofield_std2native:
            for choice in self.geofield_std2native[fieldname]:
                if choice in keys:
                    return choice
        return None

    def get_matching_dimname(self, geodimname):
        """Return the equivalent name in the native format for a standard
        dimension.

        This is a translation of the standard names to native ones. It is used
        for internal purpose only and should not be called directly.

        The standard dimension names are:

        * x, y, time for :class:`~cerbere.datamodel.grid.Grid`
        * row, cell, time for :class:`~cerbere.datamodel.swath.Swath` or
          :class:`~cerbere.datamodel.image.Image`

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        Args:
            dimname (str): standard dimension name.

        Returns:
            str: return the native name for the dimension. Return `dimname` if
                the input dimension has no standard name.

        See Also:
            see :func:`get_standard_dimname` for the reverse operation
        """
        dims = self.get_handler().dimensions.keys()
        # case where file is a reference (possibly no time value defined)
        # => creates a dummy time dimension
        if self.is_reference:
            dims.append('time')
        if geodimname in self.geodim_std2native:
            for choice in self.geodim_std2native[geodimname]:
                if choice in dims:
                    return choice
        return None

    def get_standard_dimname(self, geodimname):
        """
        Returns the equivalent standard dimension name for a
        dimension in the native format.

        This is a translation of the native names to standard ones. It is used
        for internal purpose and should not be called directly.

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        Args:
            dimname (string): native dimension name

        Return:
            str: the (translated) standard name for the dimension. Return
            `dimname` if the input dimension has no standard name.

        See Also:
            see :func:`get_matching_dimname` for the reverse operation
        """
        if geodimname in self.geodim_native2std:
            return self.geodim_native2std[geodimname]
        else:
            return geodimname

    #@profile
    def _open_netcdf(self):
        """
        Open a NetCDF file
        """
        url = self._url
        logging.debug("opening in mode : %s", self._mode)
        logging.debug("url : %s", url)
        cd2dir = False
        if len(url) >= 256:
            # netCDF4 fails if filename is longer than 256 char
            cd2dir = True
        try:
            # move to file directory if filename is too long
            if cd2dir:
                localpath = os.path.curdir
                os.chdir(os.path.dirname(url))
                url = os.path.basename(url)
            if not self.is_readonly():
                logging.debug("format : %s", self.format)
                handler = netCDF4.Dataset(url, self._mode, format=self.format)
            else:
                handler = netCDF4.Dataset(url)
            self._handler = handler
            if cd2dir:
                # move back to original rep
                os.chdir(localpath)
            return handler
        except:
            raise CorruptFileException("Could not read file %s", self._url)

    def open(self,
             view=None,
             datamodel=None,
             datamodel_geolocation_dims=None):
        """Open the netCDF file

        Args:
            view (dict, optional): a dictionary where keys are dimension names
                and values are slices. A view can be set on a file, meaning
                that only the subset defined by this view will be accessible.
                This view is expressed as any subset (see :func:`get_values`).
                For example::

                view = {'time':slice(0,0), 'lat':slice(200,300),
                'lon':slice(200,300)}

            datamodel (str): type of feature read or written. Internal argument
                only used by the classes from :mod:`~cerbere.datamodel`
                package. Can be 'Grid', 'Swath', etc...

            datamodel_geolocation_dims (list, optional): list of the name of the
                geolocation dimensions defining the data model to be read in
                the file. Optional argument, only used by the datamodel
                classes, in case the mapper class can store different types of
                data models.

        Returns:
            an handler on the opened file
        """
        if self._url is not None:
            handler = self._open_netcdf()
            self._handler = handler
        else:
            if self._mode == WRITE_NEW:
                raise Exception('Can not open file time series in write mode')
            return self._urlseries
        if datamodel is not None:
            self._feature_type = datamodel
        if datamodel_geolocation_dims:
            self.datamodel_geolocation_dims = datamodel_geolocation_dims
        self.view = view
        if self.center_on_greenwhich:
            # check if file content is indeed shifted
            lonvar = self.get_geolocation_field('lon')
            if len(handler.variables[lonvar].dimensions) > 1:
                raise Exception(
                    "This is not a grid or it is not in cylindrical projection"
                    ". Grid centering can not be applied here")
            firstlat = handler.variables[lonvar][0]
            lastlat = handler.variables[lonvar][-1]
            lastlat = lastlat if lastlat <= 180. else lastlat - 360.
            if firstlat < lastlat:
                raise Exception("Grid does not seem to be de-centered")
        return handler

    def close(self):
        """Close file"""
        self._url = None
        if self._handler is not None:
            self._handler.close()
            self._handler = None

    def sync(self):
        """force physical writing of the content on disk."""
        if self._handler is not None:
            self._handler.sync()

    def get_data_reference(self, time=None, proximity='exact'):
        '''
        return a complete reference pointer to the requested data

        time can be either a single time value or a time interval
        '''
        if self._url is not None:
            raise Exception("Not implemented")
        elif not self._urlseries is None:
            if isinstance(time, tuple):
                pass
            # case where a single time step is requested
            elif isinstance(time, datetime.datetime):
                # case where data span over several files
                if not self._urlseries is None:
                    url = self._urlseries.get_url_for_time(self, time, proximity=proximity)
                    if url is None:
                        return None
                    # warning : opening/closing many times the same file 
                    # can introduce serious performance issues
                    dataset = netCDF4.Dataset(url)
                    idx = netCDF4.date2index(time, dataset.variables['time'], select='exact')
                    reference = FieldHandler(url, 'time', index=idx, status=STORED)
                    return reference
            return
        else:
            raise Exception()

    def _do_remove_time(self, stddims):
        """Check if time dimension must be removed."""
        return ((self._feature_type == 'Grid' and
                 'x' in stddims and 'y' in stddims and 'time' in stddims and
                 'bnds' not in stddims) or
                (self._feature_type == 'Swath' and
                 'time' in stddims and len(stddims) > 1))

    #@profile
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
        assert isinstance(slices, (type(None), list, dict))
        native_fieldname = self.get_geolocation_field(fieldname)
        if native_fieldname is None:
            native_fieldname = fieldname
        if self._handler is None:
            self._handler = self.get_handler()
        # dummy time for reference grids
        if self.is_reference and native_fieldname == 'time'\
                and native_fieldname not in self.get_handler().variables:
            return numpy.ma.array([0.])

        var = self._handler.variables[native_fieldname]
        var.set_auto_maskandscale(True)
        dims = list(self._handler.variables[native_fieldname].dimensions)
        stddims = []
        for dim in dims:
            stddim = self.get_standard_dimname(dim)
            stddims.append(stddim if stddim is not None else dim)
        if slices is None and self.view is None and\
                not self.center_on_greenwhich:
            # no slicing whatsoever
            values = numpy.ma.array(var[:])
        else:
            dimsizes = [len(self._handler.dimensions[dim]) for dim in dims]
            # view slicing
            if self.view is not None:
                viewslices = Field.format_slices(self.view, dims)
                viewslices = self._fill_slices(viewslices, dimsizes)
            # standard slicing
            if slices is not None:
                ncslices = copy.copy(slices)
                if (fieldname == 'time' and
                        self.__is_swath() and
                        len(self.get_full_dimensions('time')) == 1):
                    # case of swath files where the time is one-dimensional
                    # along the satellite track (row dimension), e.g. one time
                    # value per scan line : we reconstruct a 2-dimensionsal
                    # time field (row, cell) to match the cerbere swath model
                    if slices is not None:
                        ncslices = [slices[0]]
                    else:
                        ncslices = slices
                if isinstance(ncslices, dict):
                    if self._do_remove_time(stddims) and 'time' in ncslices:
                        ncslices.pop('time')
                    ncslices = Field.format_slices(ncslices, stddims)
                elif isinstance(ncslices, list):
                    if self._do_remove_time(stddims):
                        ncslices.insert(stddims.index('time'), slice(None))
                if self.view is None:
                    ncslices = self._fill_slices(ncslices, dimsizes)
                else:
                    # slices are relative to the opened view
                    viewdimsizes = [self.get_dimsize(dim) for dim in dims]
                    ncslices = self._fill_slices(ncslices, viewdimsizes)
                    for index, (ncsli, vwsli) in enumerate(zip(ncslices,
                                                               viewslices)):
                        vw_first = vwsli.start
                        nc_first = ncsli.start
                        if ncsli.stop is None:
                            nc_last = 0
                        else:
                            nc_last = ncsli.stop - ncsli.step / abs(ncsli.step)
                        first = vw_first + nc_first * vwsli.step
                        last = vw_first + nc_last * vwsli.step
                        start = first
                        step = vwsli.step * ncsli.step
                        stop = last + step / abs(step)
                        if stop == -1:
                            stop = None
                        ncslices[index] = slice(start, stop, step)
            elif self.view is not None:
                ncslices = viewslices
            else:
                ncslices = [slice(None) for i in range(len(dims))]
                ncslices = self._fill_slices(ncslices, dimsizes)
            if self.center_on_greenwhich:
                # shift offsets for longitudes when centering on meridian 0
                # is requested
                if self._feature_type in ['Grid', 'GridTimeSeries']:
                    xdim = self.get_matching_dimname('x')
                    if xdim in dims:
                        ind = dims.index(xdim)
                        xdimsize = self.get_dimsize(xdim)
                        ## Version with indices list
                        # sli = ncslices[ind]
                        # if sli.stop is not None:
                        #     indices = numpy.arange(sli.start, sli.stop, sli.step)
                        # else:
                        #     indices = numpy.arange(sli.start, -1, sli.step)
                        # indices = numpy.mod(indices + xdimsize / 2, xdimsize)
                        # ncslices[ind] = indices
                        # values = numpy.ma.array(var[ncslices])
                        ## Version with slice
                        start = (ncslices[ind].start + xdimsize / 2) % xdimsize
                        if ncslices[ind].stop is not None:
                            stop = (ncslices[ind].stop + xdimsize / 2) % xdimsize
                        else:
                            stop = -1 + xdimsize / 2
                        step = ncslices[ind].step
                        if step > 0 and stop <= start:
                            left_slice = copy.copy(ncslices)
                            left_slice[ind] = slice(start, xdimsize, step)
                            right_start = step - (xdimsize-1-start) % step - 1
                            if right_start >= stop:
                                values = numpy.ma.array(var[left_slice])
                            else:
                                right_slice = copy.copy(ncslices)
                                right_slice[ind] = slice(right_start, stop, step)
                                values = numpy.ma.concatenate(
                                    [var[left_slice], var[right_slice]],
                                    axis=ind
                                )
                        elif step < 0 and start <= stop:
                            left_slice = copy.copy(ncslices)
                            left_slice[ind] = slice(start, None, step)
                            right_start = xdimsize - (abs(step) - start % abs(step))
                            if right_start <= stop:
                                values = numpy.ma.array(var[left_slice])
                            else:
                                right_slice = copy.copy(ncslices)
                                right_slice[ind] = slice(right_start, stop, step)
                                values = numpy.ma.concatenate(
                                    [var[left_slice], var[right_slice]],
                                    axis=ind
                                )
                        else:
                            ncslices[ind] = slice(start, stop, step)
                            values = numpy.ma.array(var[ncslices])
                    else:
                        values = numpy.ma.array(var[ncslices])
                else:
                    raise Exception("Longitude centering only applies to grid features")
            else:
                values = numpy.ma.array(var[ncslices])
        if fieldname == 'lon' or fieldname == 'longitude':
            ind = numpy.ma.where(values >= 180.)
            values[ind] = values[ind] - 360.
        # apply additional fillvalue
        if self.fillvalue is not None:
            values = numpy.ma.masked_values(values, self.fillvalue, copy=False)
        # remove the time dimension from two-dimensional fields in Grid
        # features (Grid model in cerbere is (x, y))
        if self._do_remove_time(stddims):
            shape = list(values.shape)
            shape.pop(stddims.index('time'))
            values = values.reshape(tuple(shape))
        if (fieldname == 'time' and
                self.__is_swath() and
                len(self.get_full_dimensions('time')) == 1):
            # case of swath files where the time is one-dimensional
            # along the satellite track (row dimension), e.g. one time value
            # per scan line : we reconstruct a 2-dimensionsal time field
            # (row, cell) to match the cerbere swath model
            rows = self.get_dimsize('row')
            cols = self.get_dimsize('cell')
            if slices is None:
                shape = (cols, rows)
            else:
                newslices = self._fill_slices(slices, (rows, cols))
                shape = (newslices[1].stop - newslices[1].start,
                         newslices[0].stop - newslices[0].start)
            time = numpy.ma.resize(values, shape).transpose()
            return time
        return values

    def read_fillvalue(self, fieldname):
        """
        """
        native_fieldname = self.get_geolocation_field(fieldname)
        if native_fieldname is None:
            native_fieldname = fieldname
        var = self.get_handler().variables[native_fieldname]
        if '_FillValue' in var.__dict__:
            return var._FillValue
        elif 'fill_value' in var.__dict__:
            return var.fill_value
        else:
            return None

    @classmethod
    def get_netcdf_fillvalue(cls, datatype):
        if datatype == str:
            return str('')
        elif datatype == numpy.dtype('S1'):
            return ' '
        elif datatype == numpy.dtype(numpy.complex64):
            # take type of real and imag components
            comptype = numpy.dtype('float32')
            strtype = '%s%s' % (comptype.kind, comptype.itemsize)
        elif datatype == numpy.dtype(numpy.complex128):
            # take type of real and imag components
            comptype = numpy.dtype('float64')
            strtype = '%s%s' % (comptype.kind, comptype.itemsize)
        else:
            strtype = '%s%s' % (datatype.kind, datatype.itemsize)
        return netCDF4.default_fillvals[strtype]

    def _create_field(self, name, datatype, dims, fillvalue, comp, field):
        dropAttributes = ['dtype']

        if self.format != 'NETCDF3_CLASSIC':
            try:
                if self.format == 'NETCDF4' and datatype == 'int32':
                    # cast from int32 to int64 or the file appeared corrupted
                    # later (bug in netcdf library?).
                    datatype = 'int64'
                    fillvalue = numpy.int64(fillvalue)
                if self.format == 'NETCDF4' and datatype == str:
                    fillvalue = None

                result = self.get_handler().createVariable(
                    name,
                    datatype=datatype,
                    dimensions=dims,
                    zlib=True,
                    complevel=6,
                    shuffle=True,
                    fletcher32=False,
                    fill_value=fillvalue
                    )
            except:
                logging.error("Could not create field %s", name)
                raise
        else:
            result = self.get_handler().createVariable(
                name,
                datatype=datatype,
                dimensions=dims,
                fill_value=fillvalue
            )

        # add variable attributes
        if comp.variable.description:
            result.long_name = comp.variable.description
        if comp.variable.standardname:
            result.standard_name = comp.variable.standardname
        if comp.units:
            result.units = comp.units
        if comp.valid_min:
            result.valid_min = comp.valid_min
        if comp.valid_max:
            result.valid_max = comp.valid_max
        for att in comp.attributes:
            # Stop transferring some attributes which are redundant or
            # incorrect.
            if att.lower() not in dropAttributes:
                result.setncattr(att, comp.attributes[att])
        # create quality fields if defined
        if field.qc_levels:
            result = self.get_handler().createVariable(
                comp.variable.shortname + '_qc_level',
                'i1',
                dimensions=tuple(field.qc_levels.dimensions.keys()),
                fill_value=-127
            )
            result.flag_values = field.qc_levels.levels
            result.flag_meanings = field.qc_levels.meanings
        if field.qc_details:
            result = self.get_handler().createVariable(
                comp.variable.shortname + '_qc_details',
                'i2',
                dimensions=tuple(field.qc_details.dimensions.keys()),
                fill_value=-32767
            )
            result.flag_masks = field.qc_details.mask
            result.flag_meanings = field.qc_details.meanings

    def _create_complex(self,
                        varName,
                        dType,
                        dims,
                        comp,
                        field,
                        fillValue=None):
        """Internal wrapper for the creation of a complex field."""
        if fillValue is None:
            fillValue = self.get_netcdf_fillvalue(dType)
        names = (varName + '_real', varName + '_imag')
        for name in names:
            self._create_field(name, dType, dims, fillValue, comp, field)

    def create_field(self, field, dim_translation=None, feature=None):
        """
        Create netCDF variable for provided field

        Args:
            field (:class:`~cerbere.datamodel.field.Field`): field to create.

            dim_translation (dict, optional): dictionary, alternative dimension
                 names to be used.

            feature (string, optional): class name of the feature. Used to
                implement specific writing conventions for some features
                (ex: additional dummy time dimension for grids).
        """
        try:
            if self.is_readonly():
                raise Exception("File is read only")
            for comp in field.get_components():
                dims = comp.dimensions.keys()
                if dim_translation:
                    # renaming of default dimensions
                    dims = list(comp.dimensions)
                    for idim, dim in enumerate(dims):
                        if dim in dim_translation:
                            dims[idim] = dim_translation[dim]
                # add a dummy time dimension for grids (CF convention) to all
                # fields except lat/lon (and to time only if it is a 2-d field)
                atemporal_fields = ['lat', 'lon', 'lat_bnds', 'lon_bnds']
                if feature == 'Grid' and\
                        (field.variable.shortname not in atemporal_fields or
                         (field.variable.shortname == 'time' and
                          len(field.dimensions) > 1)):
                    if 'time' not in dims:
                        dims.insert(0, 'time')
                fillvalue = None
                if comp.fillvalue is not None:
                    fillvalue = comp.fillvalue

                if comp.datatype == numpy.dtype(numpy.complex128):
                    # netcdf can not store complex and the field has to be
                    # split into two netcdf variables
                    self._create_complex(comp.variable.shortname,
                                         numpy.dtype(numpy.float64),
                                         dims, comp, field)
                elif comp.datatype == numpy.dtype(numpy.complex64):
                    # netcdf can not store complex and the field has to be
                    # split into two netcdf variables
                    self._create_complex(comp.variable.shortname,
                                         numpy.dtype(numpy.float32),
                                         dims, comp, field)
                else:
                    # 1.  Need to convert datetime object to int32
                    if comp.datatype == datetime.datetime:
                        datatype = numpy.dtype(numpy.int32)
                        fillvalue = self.get_netcdf_fillvalue(datatype)
                    # 2.  NETCDF4 cannot handle int64
                    elif self.format != "NETCDF4" and\
                            comp.datatype == numpy.dtype(numpy.int64):
                        datatype = numpy.dtype(numpy.int32)
                        if fillvalue is None:
                            fillvalue = self.get_netcdf_fillvalue(datatype)
                        else:
                            fillvalue = numpy.int32(fillvalue)
                    else:
                        datatype = comp.datatype
                        if fillvalue is None:
                            fillvalue = self.get_netcdf_fillvalue(datatype)

                    self._create_field(comp.variable.shortname, datatype, dims,
                                       fillvalue, comp, field)
            return True
        except:
            logging.error("Error when creating field %s",
                          field.variable.shortname)
            raise

# 
# 
#     def create_field(self, field, dim_translation=None):
#         """
#         Create netCDF variable for provided field
# 
#         dim_translation: dictionary, alternative dimension names to be used
#         """
#         try:
#             if self.is_readonly():
#                 raise Exception("File is read only")
#             ncvars = []
#             varnames = []
#             vardatatypes = []
#             for comp in field.get_components():
#                 if comp.datatype == numpy.dtype(numpy.complex128):
#                     # netcdf can not store complex and the field has to be 
#                     # splitted in two netcdf variables
#                     ncvars.extend([comp, comp])
#                     varnames.extend([comp.variable.shortname + '_real',
#                                      comp.variable.shortname + '_imag'])
#                     vardatatypes.extend([numpy.dtype(numpy.float64),
#                                          numpy.dtype(numpy.float64)])
#                 elif comp.datatype == numpy.dtype(numpy.complex64):
#                     # netcdf can not store complex and the field has to be 
#                     # splitted in two netcdf variables
#                     ncvars.extend([comp, comp])
#                     varnames.extend([comp.variable.shortname + '_real',
#                                      comp.variable.shortname + '_imag'])
#                     vardatatypes.extend([numpy.dtype(numpy.float32),
#                                          numpy.dtype(numpy.float32)])
#                 else:
#                     ncvars.append(comp)
#                     varnames.append(comp.variable.shortname)
#                     vardatatypes.append(comp.datatype)
#             for i, comp in enumerate(ncvars):
#                 units = None
#                 #print comp.variable.shortname, comp.datatype
#                 if comp.datatype == type(datetime.datetime.now()):
#                     # defaut coding for time values
#                     datatype = numpy.dtype(numpy.int32)
#                     units = comp.units
#                 else:
#                     datatype = comp.datatype
#                 # forbids int64
#                 if self.format != "NETCDF4":
#                     if datatype == numpy.dtype(numpy.int64):
#                         logging.warning("Forced type from int64 to int32")
#                         datatype = numpy.dtype(numpy.int32)
#                 if comp.fillvalue is None:
#                     fillvalue = self.get_netcdf_fillvalue(datatype)
#                 else:
#                     fillvalue = comp.fillvalue
#                 dims = comp.dimensions.keys()
#                 if dim_translation:
#                     # renaming of default dimensions
#                     dims = list(comp.dimensions)
#                     for idim, dim in enumerate(dims):
#                         if dim in dim_translation:
#                             dims[idim] = dim_translation[dim]
#                 logging.debug(
#                             "CREATE : %s %s %s %s",
#                             varnames[i], vardatatypes[i], dims, fillvalue
#                             )
#                 if self.format != 'NETCDF3_CLASSIC':
#                     try:
#                         result = self.get_handler().createVariable(
#                                     varnames[i],
#                                     datatype=vardatatypes[i],
#                                     dimensions=dims,
#                                     zlib=True,
#                                     complevel=6,
#                                     shuffle=True,
#                                     fletcher32=False,
#                                     fill_value=fillvalue
#                                     )
#                     except:
#                         logging.error("Could not create field %s",
#                                       varnames[i])
#                         raise
#                 else:
#                     result = self.get_handler().createVariable(
#                                     varnames[i],
#                                     datatype=vardatatypes[i],
#                                     dimensions=dims,
#                                     fill_value=fillvalue
#                                     )
#                 # add variable attributes
#                 if comp.variable.description:
#                     result.long_name = comp.variable.description
#                 if comp.variable.standardname:
#                     result.standard_name = comp.variable.standardname
#                 if comp.units and units is None:
#                     units = comp.units
#                 if units:
#                     result.units = units
#                 if comp.valid_min:
#                     result.valid_min = comp.valid_min
#                 if comp.valid_max:
#                     result.valid_max = comp.valid_max
#                 for att in comp.attributes:
#                     result.setncattr(att, comp.attributes[att])
#                 # create quality fields if defined
#                 if field.qc_levels:
#                     result = self.get_handler().createVariable(
#                                         comp.variable.shortname + '_qc_level',
#                                         'i1',
#                                         dimensions=tuple(field.qc_levels.dimensions.keys()),
#                                         fill_value=-127
#                                         )
#                     result.flag_values = field.qc_levels.levels
#                     result.flag_meanings = field.qc_levels.meanings
#                 if field.qc_details:
#                     result = self.get_handler().createVariable(
#                                         comp.variable.shortname + '_qc_details',
#                                         'i2',
#                                         dimensions=tuple(field.qc_details.dimensions.keys()),
#                                         fill_value=-32767
#                                         )
#                     result.flag_masks = field.qc_details.mask
#                     result.flag_meanings = field.qc_details.meanings
#             return True
#         except:
#             logging.error("Error when creating field %s",
#                           field.variable.shortname)
#             raise

    def write_field(self, field):
        """Write the data of a field in the netCDF file.

        Args:
            field (:class:`~cerbere.datamodel.field.Field`): the field which
                content has to be written in the file.

        .. Warning::
           The field with all its metadata (attributes) must have been
           explicitely created by the mapper before, by calling the
           :func:`create_field` method of this class.
        """
        try:
            values = field.get_values()
            if values is None:
                raise Exception("No values provided for field ", field.name)
            # forbids int64
            if self.format != "NETCDF4":
                if values.dtype == numpy.dtype(numpy.int64):
                    logging.warning("Forced type of data from int64 to int32")
                    values = values.astype(numpy.int32)
            if field.name == 'time':
                # Need to convert incoming datetime object to a floating point
                # time as netcdf cannot handle time properly.
                if isinstance(values.reshape(-1)[0], datetime.datetime):
                    units = DEFAULT_TIME_UNITS
                    if field.units is not None:
                        units = field.units
                    values = netCDF4.date2num(values, units)
                # Else already in proper units (we hope).

            logging.debug("Writing field %s", field.name)
            if field.datatype == numpy.dtype(numpy.complex64)\
                    or field.datatype == numpy.dtype(numpy.complex128):
                # netcdf can not store complex and the field is
                # splitted in two netcdf variables
                ncvar_real = self.get_handler().variables[field.name + '_real']
                ncvar_imag = self.get_handler().variables[field.name + '_imag']
                if ncvar_real.shape == (0,)\
                        or values.shape == ()\
                        or ncvar_real.shape == values.shape:
                    ncvar_real[:] = values.real
                    ncvar_imag[:] = values.imag
                else:
                    raise IndexError(
                        str(
                            "Your input data shape {} is not compatible with your "
                            "field shape {}."
                        ).format(
                            values.shape, ncvar_real.shape
                        )
                    )
            else:
                ncvar = self.get_handler().variables[field.name]
                if ncvar.shape == (0,)\
                        or ncvar.shape == ()\
                        or values.shape == ()\
                        or ncvar.shape == values.shape\
                        or (ncvar.shape == (1,) and values.size == 1):
                    ncvar[:] = values
                elif (len(ncvar.shape) - 1 == len(values.shape) and
                      (ncvar.shape[0] == 1 or ncvar.shape[0] == 0)):
                    # case of grids
                    ncvar[0, :] = values
                else:
                    raise IndexError(
                        str(
                            "Your input data shape {} is not compatible with your "
                            "field shape {}."
                        ).format(
                            values.shape, ncvar.shape
                        )
                    )
            if field.qc_levels is not None:
                self.get_handler().variables[field.name + '_qc_level'][:]\
                    = field.qc_levels.values[:]
            if field.qc_details is not None:
                fieldname = field.name + '_qc_details'
                self.get_handler().variables[fieldname].set_auto_maskandscale(True)
                self.get_handler().variables[fieldname][:]\
                    = field.qc_details.values[:]
        except:
            logging.error("Could not write field {}".format(field.name))
            raise

    def __is_swath(self):
        """Returns True if the stored feature is a Swath or Image."""
        return ('row' in self.get_full_dimensions() and
                'cell' in self.get_full_dimensions())

    def read_field(self, fieldname):
        """
        Return the :class:`cerbere.field.Field` object corresponding to
        the requested fieldname.

        The :class:`cerbere.field.Field` class contains all the metadata
        describing a field (equivalent to a variable in netCDF).

        Args:
            fieldname (str): name of the field

        Return:
             :class:`cerbere.field.Field`: the corresponding field object
        """
        # Creates a dummy time field for reference files with no
        # explicit time in the file
        native_fieldname = self.get_geolocation_field(fieldname)
        if self.is_reference and native_fieldname == 'time'\
                and native_fieldname not in self.get_handler().variables:
            # create a field for time
            variable = Variable(
                shortname=fieldname,
                description='reference time',
                authority=self.get_naming_authority(),
                standardname='time'
                )
            field = Field(
                variable,
                collections.OrderedDict([('time', 1)]),
                datatype=numpy.dtype(numpy.int64),
                units='seconds since 0001-01-01 00:00:00'
                )
            field.attach_storage(self.get_field_handler(fieldname))
            return field
        # Get NetCDF Variable
        naming_auth = self.get_naming_authority()
        if native_fieldname is None:
            native_fieldname = fieldname
        attrs = self.read_field_attributes(native_fieldname)
        if 'long_name' in attrs:
            descr = attrs['long_name']
        else:
            descr = None
        if 'standard_name' in attrs:
            stdname = attrs['standard_name']
        else:
            stdname = None
        var = Variable(
            shortname=fieldname,
            description=descr,
            authority=naming_auth,
            standardname=stdname
            )
        dims = self.get_full_dimensions(native_fieldname)
        if 'scale_factor' in attrs:
            datatype = numpy.dtype(attrs['scale_factor'])
        elif 'add_offset' in attrs:
            datatype = numpy.dtype(attrs['add_offset'])
        else:
            datatype = self.get_handler().variables[native_fieldname].dtype
        if '_FillValue' in attrs:
            fillvalue = attrs['_FillValue']
        elif 'missing_value' in attrs:
            fillvalue = attrs['missing_value']
        else:
            fillvalue = None
        field = Field(var, dims, datatype=datatype, fillvalue=fillvalue)
        field.attach_storage(self.get_field_handler(fieldname))
        # MetaData
        field.attributes = {}
        field.units = None
        if 'units' in attrs:
            field.units = attrs['units']
        field.valid_min = None
        field.valid_max = None
        if 'valid_min' in attrs and 'valid_max' in attrs:
            # sometimes expressed as a string : need type cast
            try:
                field.valid_min = numpy.array(attrs['valid_min']).astype(
                    field.datatype)
                field.valid_max = numpy.array(attrs['valid_max']).astype(
                    field.datatype)
                if 'scale_factor' in attrs:
                    field.valid_min = field.valid_min * attrs['scale_factor']
                    field.valid_max = field.valid_max * attrs['scale_factor']
                if 'add_offset' in attrs:
                    field.valid_min = field.valid_min + attrs['add_offset']
                    field.valid_max = field.valid_max + attrs['add_offset']
            except:
                field.valid_min = attrs['valid_min']
                field.valid_max = attrs['valid_max']
                logging.error("invalid valid_min or valid_max : %s %s",
                              field.valid_min,
                              field.valid_max)
        for att in attrs:
            if att not in ['units', 'scale_factor', 'add_offset',
                           '_FillValue', 'valid_min', 'valid_max']:
                field.attributes[att] = attrs[att]
        if (fieldname == 'time' and
                self.__is_swath() and
                len(self.get_full_dimensions('time')) == 1):
            # case of swath files where the time is one-dimensional
            # along the satellite track (row dimension), e.g. one time value
            # per scan line : we reconstruct a 2-dimensionsal time field
            # (row, cell) to match the cerbere swath model
            rows = self.get_dimsize('row')
            cols = self.get_dimsize('cell')
            field.dimensions = collections.OrderedDict([('row', rows),
                                                        ('cell', cols)])
        return field

    def create_dim(self, dimname, size=None):
        """Add a new dimension.

        Args:
            dimname (str): name of the dimension.
            size (int): size of the dimension (unlimited if None)
        """
        if not self.is_writable():
            raise RuntimeError('Dimension can only be defined in write mode')
        return self.get_handler().createDimension(dimname,
                                                  size)

    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature

        Excludes the geolocation fields in the returned list
        '''
        geophyvars = self.get_handler().variables.keys()
        for geofield in ['lat', 'lon', 'time', 'z']:
            ncgeovar = self.get_geolocation_field(geofield)
            if ncgeovar in geophyvars:
                geophyvars.remove(ncgeovar)
        return geophyvars

    def get_dimensions(self, fieldname=None):
        """
        return the (standard) dimension names of a file or field.

        Standardization of the dimension name is applied when possible.
        """
        if not fieldname:
            dims = self.get_handler().dimensions
        else:
            if fieldname not in self.get_handler().variables:
                raise Exception("Field % not existing in NetCDF file",
                                fieldname)
            dims = self.get_handler().variables[fieldname].dimensions
        stddims = []
        for dim in dims:
            stdd = self.get_standard_dimname(dim)
            if stdd is None:
                stddims.append(dim)
            else:
                stddims.append(stdd)
        # remove the time dimension from two-dimensional fields in Grid
        # features (Grid model in cerbere is (x, y))
        if self._do_remove_time(stddims):
            stddims.remove('time')
        return tuple(stddims)

    def get_full_dimensions(self, fieldname=None):
        """return the name and size of a field's dimensions as an
        ordered dictionary

        Args:
            fieldname (string, optional): name of the specific field
                (netCDF variable) for which the dimensions are requested. If
                not provided, all file dimensions are returned.

        Return:
            OrderedDict: an ordered dict where keys are the dimension
                names and values are their size
        """
        dims = collections.OrderedDict()
        for dim in self.get_dimensions(fieldname):
            standarddim = self.get_standard_dimname(dim)
            dimsize = self.get_dimsize(dim)
            if standarddim:
                dims[standarddim] = dimsize
            else:
                dims[dim] = dimsize
        return dims

    def get_dimsize(self, dimname):
        """returns the size of a dimension

        If a view was set when opening the file, the size of the view subset
        is returned.

        Args:
            dimname (string): name of the dimension for which the size is
                inquired. Can be provided as the native or standard dimension
                name.

        Return:
            int: size of the dimension. 
        """
        dim = self.get_matching_dimname(dimname)
        if dim is None:
            dim = dimname
        dimsize = len(self.get_handler().dimensions[dim])
        if self.view is not None and dim in self.view:
            viewslice = self._fill_slices([self.view[dim]], [dimsize])[0]
            start, stop, step = viewslice.start, viewslice.stop, viewslice.step
            if stop is None:
                stop = -1
            dimsize = 1 + (abs(stop - start) - 1) / abs(step)
        return dimsize

    def read_field_attributes(self, fieldname):
        """
        return the specific storage attributes of a field
        (such as _FillValue, scale_factor, add_offset,...)

        :param fieldname: name of the field
        :type fieldname: str

        :return: a dictionary where keys are the attribute names.
        :rtype: Dict
        """
        return copy.copy(self.get_handler().variables[fieldname].__dict__)

    def read_global_attributes(self):
        """Returns the names of the global attributes.

        Returns:
            list<str>: the list of the attribute names.
        """
        return self.get_handler().ncattrs()

    def write_global_attributes(self, attrs):
        """Write the global attributes of the file.

        Includes a set of predefined standard attributes for CF compliant
        files. Default attributes will be filled automatically with default
        values unless you override them in `attrs` argument. You can also
        add additional attributes in `attrs` dictionary.

        Args:
            attrs (dict<string, string or number or datetime>): a dictionary
                containing the attributes names and values to be written.
        """
        attributes = copy.copy(attrs)
        dataset = self.get_handler()
        if self._mode == WRITE_NEW:
            dataset.Conventions = self.CONVENTIONS
            dataset.netcdf_version_id = netCDF4.getlibversion()
            dataset.date_created\
                = datetime.datetime.now().strftime(self.TIME_FORMAT)
            dataset.date_modified\
                = datetime.datetime.now().strftime(self.TIME_FORMAT)
            # convert datetime objects to string
            for att in attributes:
                if isinstance(attributes[att], datetime.datetime):
                    attributes[att]\
                        = attributes[att].strftime(self.TIME_FORMAT)
            # we try to keep some orderinf in the attributes
            for att in self.STANDARD_ATTRIBUTES_VALUES:
                if att in attributes:
                    dataset.setncattr(att, attributes[att])
                else:
                    dataset.setncattr(att,
                                      self.STANDARD_ATTRIBUTES_VALUES[att])
            # additional attributes
            keys = attributes.keys()
            keys.sort()
            for att in keys:
                if att not in self.STANDARD_ATTRIBUTES_VALUES:
                    dataset.setncattr(att, attributes[att])
        else:
            dataset.date_modified\
                = datetime.datetime.now().strftime('%Y%m%dT%H%M%SZ')
            for att in attributes:
                dataset.setncattr(att, attributes[att])
        return

    def read_global_attribute(self, name):
        if name in self.get_handler().ncattrs():
            return self.get_handler().getncattr(name)
        else:
            return None

    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage.

        Returns:
            datetime: start time of the data in file.
        """
        handler = self.get_handler()
        attrs = handler.ncattrs()
        if 'start_date' in attrs:
            attrdate = handler.getncattr('start_date').replace(' UTC', '')
            # remove Zulu time indication
            attrdate = attrdate.strip('Z')
            if 'start_time' in attrs:
                attrtime = handler.getncattr('start_time')
                # remove Zulu time indication
                attrtime = attrtime.strip('Z')
                # remove UTC time indication
                attrtime = attrtime.replace(' UTC', '')
                # combine time and date
                attrdate = attrdate + 'T' + attrtime
            dt = parser.parse(attrdate)
            return dt
        elif 'time_coverage_start' in attrs:
            try:
                attrdate = handler.getncattr('time_coverage_start')
                attrdate = attrdate.strip('Z')
                return parser.parse(attrdate)
            except ValueError:
                logging.error('Unexpected metadata time format.')
                # Unexpected metadata time format.
                return None
        else:
            return None

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage.

        Returns:
            datetime: end time of the data in file.
        """
        handler = self.get_handler()
        attrs = self.get_handler().ncattrs()
        if 'stop_date' in attrs:
            attrdate = handler.getncattr('stop_date').replace(' UTC', '')
            # remove Zulu time indication
            attrdate = attrdate.strip('Z')
            if 'stop_time' in attrs:
                attrtime = handler.getncattr('stop_time')
                # remove Zulu time indication
                attrtime = attrtime.strip('Z')
                # remove UTC time indication
                attrtime = attrtime.replace(' UTC', '')
                # combine time and date
                attrdate = attrdate + 'T' + attrtime
            dt = parser.parse(attrdate)
            return dt
        elif 'time_coverage_end' in attrs:
            try:
                attrdate = handler.getncattr('time_coverage_end')
                attrdate = attrdate.strip('Z')
                return parser.parse(attrdate)
            except ValueError:
                # Unexpected metadata time format.
                logging.error('Unexpected metadata time format.')
                return None
        elif 'time_coverage_stop' in attrs:
            try:
                attrdate = handler.getncattr('time_coverage_stop')
                attrdate = attrdate.strip('Z')
                return parser.parse(attrdate)
            except ValueError:
                # Unexpected metadata time format.
                logging.error('Unexpected metadata time format.')
                return None
        else:
            return None

    def get_spatial_resolution_in_deg(self):
        """Returns the average spatial resolution in degrees"""
        attrs = self.get_handler().ncattrs()
        latres = filter(lambda x: x in ['geospatial_lat_resolution',
                                        'latitude_resolution'],
                        attrs)
        lonres = filter(lambda x: x in ['geospatial_lon_resolution',
                                        'longitude_resolution'],
                        attrs)
        if len(latres) == 0:
            return None
        latres = self.get_handler().getncattr(latres[0])
        lonres = self.get_handler().getncattr(lonres[0])
        if latres == lonres:
            return float(latres)
        else:
            raise Exception("Conflict: don't know which resolution to choose")

    def get_product_version(self):
        """return the product version"""
        attrs = self.get_handler().ncattrs()
        if 'product_version' in attrs:
            return self.get_handler().getncattr('product_version')
        else:
            return None

    def get_bbox(self):
        '''
        return the bounding box of the feature, as a tuple
        (lonmin, latmin, lonmax, latmax)
        '''
        attrs = self.get_handler().ncattrs()
        choices = {'latmin': ['southernmost_latitude','geospatial_lat_min','south_latitude'],
                   'latmax': ['northernmost_latitude','geospatial_lat_max','north_latitude'],
                   'lonmin': ['westernmost_longitude','geospatial_lon_min','west_longitude'],
                   'lonmax': ['easternmost_longitude','geospatial_lon_max','east_longitude']}
        bbox = []
        for limit in ['lonmin', 'latmin', 'lonmax', 'latmax']:
            for att in choices[limit]:
                fatt = next((x for x in attrs if x == att), None)
                if fatt != None:
                    bbox.append(self.get_handler().getncattr(fatt))
        if len(bbox) != 4:
            return None
        # case ECMWF converted from Grib!
        if bbox[0] == 360. and bbox[2] == 0.:
            bbox[0] = 0.
            bbox[2] = 360.
        return bbox

    def get_cycle_number(self):
        if 'cycle_number' in self.get_handler().ncattrs():
            return self.get_handler().getncattr('cycle_number')
        else:
            return None

    def get_orbit_number(self):
        if 'pass_number' in self.get_handler().ncattrs():
            return self.get_handler().getncattr('pass_number')
        else:
            return None
