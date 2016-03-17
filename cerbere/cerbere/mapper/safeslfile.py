# -*- coding: utf-8 -*-
"""
.. module::cerbere.mapper.safeslfile

Class to read Sentinel-3 SLSTR files (except L2P).

:copyright: Copyright 2015 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
import os
import glob
import datetime
from numpy import ma, dtype, int64
from collections import OrderedDict

from netCDF4 import num2date, date2num

from .. import READ_ONLY
from .abstractmapper import AbstractMapper
from cerbere.datamodel.variable import Variable
from cerbere.datamodel.field import Field
from cerbere.mapper.ncfile import NCFile
from cerbere.mapper.ghrsstncfile import GHRSSTNCFile
import cerbere.mapper.slices

METEO = 'met_tx.nc'

PREFIX = {'n': 'Nadir', 'o': 'Oblique', 'x': ''}

PIXSYNC = {'in': 'PIXSYNC_i', 'io': 'PIXSYNC_i', 'an': 'PIXSYNC_a'}


class SAFESLFile(AbstractMapper):
    """Abstract class for SAFE SLSTR files (except L2P).
    """
    def __init__(self, url=None, sltype=None, mode=READ_ONLY, **kwargs):
        if mode != READ_ONLY:
            raise Exception("This mapper can only be used in read_only mode.")
        if sltype not in ['i', 'c', 'a', 'b']:
            raise Exception("Unknown SLSTR product type")
        super(SAFESLFile, self).__init__(url=url, mode=mode, **kwargs)
        self.__sltype = sltype
        self.__data_handlers = []
        self.__oblique_fields = []
        # coordinate files
        geodetic = "geodetic_%sn.nc" % sltype
        times = "time_%sn.nc" % sltype
        coordinate_files = [geodetic, times]
        # detect the data files and instanciate mappers for each one
        datafiles = []
        for fname in glob.glob(os.path.join(url, "*_%s[n,o].nc" % sltype)):
            if os.path.basename(fname) not in coordinate_files:
                datafiles.append(os.path.basename(fname))
        for f in datafiles:
            fname = os.path.join(url, f)
            self.__data_handlers.append(NCFile(url=fname, mode=mode, **kwargs))
        # instantiate mappers for each coordinate file
        self.__geod_handler = NCFile(os.path.join(url, geodetic),
                                     mode=mode, **kwargs)
        self.__time_handler = NCFile(os.path.join(url, times),
                                     mode=mode, **kwargs)
        self.__fieldlocator = {}
        self.__geofieldlocator = {}
        self.__fieldtranslate = {}
        # offset between nadir and oblique swath edge
        self.nadir_to_oblique_offset = None

    def __is_oblique(self, fieldname):
        """Test if a field corresponds to an oblique view subproduct.

        Returns:
            bool: True f a field corresponds to an oblique view
        """
        return (fieldname in self.__oblique_fields)

    def open(self,
             view=None,
             datamodel=None,
             datamodel_geolocation_dims=None):
        """
        Args:
            view (dict, optional): a dictionary where keys are dimension names
                and values are slices. A view can be set on a file, meaning
                that only the subset defined by this view will be accessible.
                This view is expressed as any subset (see :func:`get_values`).
                For example::

                view = {'row':slice(200,250), 'cell':slice(200,300)}

            datamodel (str): type of feature read or written. Internal argument
                only used by the classes from :mod:`~cerbere.datamodel`
                package. Can be 'Grid', 'Swath', etc...

            datamodel_geolocation_dims (list, optional): list of the name of
                the geolocation dimensions defining the data model to be read
                in the file. Optional argument, only used by the datamodel
                classes, in case the mapper class can store different types of
                data models.

        Returns:
            an handler on the opened file
        """
        # open each related file in the SAFE repo
        if view is None:
            rowview = None
        else:
            rowview = {'row': view['row']}
            # modify view for oblique fields which are narrower
            if 'cell' in view:
                obliqueview = view
                obliqueview['cell'] = (
                    view['cell'] - self.nadir_to_oblique_offset)
        for hdlr in self.__data_handlers:
            f = os.path.basename(hdlr.get_url())
            is_oblique = (f[-4] == 'o')
            newview = view
            if is_oblique and view is not None:
                newview = obliqueview
            hdlr.open(newview, datamodel,
                      datamodel_geolocation_dims)
        self.__geod_handler.open(view, datamodel,
                                 datamodel_geolocation_dims)
        self.__time_handler.open(rowview, datamodel,
                                 datamodel_geolocation_dims)
        # build the two-way dictionaries of fields
        # ...for data
        for hdlr in self.__data_handlers:
            f = os.path.basename(hdlr.get_url())
            is_oblique = (f[-4] == 'o')
            for fieldname in hdlr.get_fieldnames():
                ncvar = hdlr.get_handler().variables[fieldname]
                if 'long_name' in (ncvar.ncattrs()):
                    longname = ncvar.long_name
                    newfieldname = self.__get_fieldname(fieldname,
                                                        longname)
                else:
                    newfieldname = fieldname
                self.__fieldtranslate[newfieldname] = fieldname
                self.__fieldlocator[newfieldname] = hdlr
                if is_oblique:
                    self.__oblique_fields.append(newfieldname)
        # ...for geodetic coordinates
        for fieldname in self.__geod_handler.get_fieldnames():
            self.__geofieldlocator[fieldname] = self.__geod_handler
        # define nadir/oblique offset
        nadir_track_offset = (
            self.__geofieldlocator['latitude_orphan_%sn' % self.__sltype]
            .read_global_attribute('track_offset'))
        oblique_track_offset = (
            self.__fieldlocator['latitude_orphan_%so' % self.__sltype]
            .read_global_attribute('track_offset'))
        self.nadir_to_oblique_offset = int(round(
            nadir_track_offset - oblique_track_offset))

    def close(self):
        """Close handler on storage"""
        for hdlr in self.__data_handlers:
            hdlr.close()
        self.__data_handlers = None
        self.__geod_handler.close()
        self.__time_handler.close()
        self.__geod_handler = None
        self.__time_handler = None

    def get_dimsize(self, dimname):
        """Return the size of a dimension.

        Args:
            dimname (str): name of the dimension.

        Returns:
            int: size of the dimension.
        """
        dim = self.get_matching_dimname(dimname)
        return self.__geod_handler.get_dimsize(dim)

    def get_dimensions(self, fieldname=None):
        """Return the dimension names of a file or a field in the
        file. For temporal and spatial dimensions, the cerbere standard names
        are returned.

        Args:
            fieldname (str): the name of the field from which to get the
                dimensions. For a geolocation field, use the cerbere standard
                name (time, lat, lon), though native field name will work too.

        Returns:
            tuple<str>: the standard dimensions of the field or file.
        """
        if fieldname is None:
            dims = self.__geod_handler.get_dimensions()
            if 'orphan_pixels' in dims:
                # same dimension name (but not size) in oblique/nadir views so
                # we have to create two dimensions since we merge the two views
                newdims = []
                for dim in dims:
                    if dim != 'orphan_pixels':
                        newdims.append(dim)
                    else:
                        newdims.extend([dim + '_n', dim + '_o'])
                return newdims
            return dims
        if fieldname in ['time', 'lat', 'lon', 'z']:
            # Should all have the same dimension as lat
            native_fieldname = self.get_geolocation_field('lat')
            dims = self.__geod_handler.get_dimensions(native_fieldname)
        else:
            handler = self.__fieldlocator[fieldname]
            dims = handler.get_dimensions(
                self.__get_native_fieldname(fieldname))
        # convert geolocation dims to standard names
        newdims = []
        for dim in list(dims):
            if self.__is_oblique(fieldname):
                view = 'o'
            else:
                view = 'n'
            newdims.append(self.get_standard_dimname(dim, view))
        return tuple(newdims)

    def get_matching_dimname(self, dimname):
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
        matching = {'time': 'time', 'row': 'rows', 'cell': 'columns',
                    'z': 'elevation'}
        # remove the oblique/nadir suffix
        if dimname[-2:] in ['_n', '_o']:
            dimname = dimname.strip('_n').strip('_o')
        if dimname in matching:
            return matching[dimname]
        return dimname

    def get_standard_dimname(self, dimname, view=None):
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
        matching = {'time': 'time', 'rows': 'row', 'columns': 'cell',
                    'elevation': 'z'}
        if dimname in matching:
            return matching[dimname]
        # for other dimensions, add the oblique/nadir suffix
        if dimname in ['row', 'cell', 'time', 'z']:
            return dimname
        if view is 'n':
            dimname += '_n'
        elif view is 'o':
            dimname += '_o'
        return dimname

    def get_fieldnames(self):
        """Returns the list of geophysical fields stored for the feature.

        The geolocation field names are excluded from this list.

        Returns:
            list<string>: list of field names
        """
        return self.__fieldlocator.keys()

    def __get_native_fieldname(self, fieldname):
        """Returns the native name of a field.

        Args:
            fieldname (str): name of the field.

        Returns:
            str: the native name of the field. The same as input
                if the field is not a geolocation field.
        """
        if fieldname in ['lat', 'lon', 'time', 'z']:
            return self.get_geolocation_field(fieldname)
        if fieldname in self.__fieldtranslate:
            return self.__fieldtranslate[fieldname]
        return fieldname

    def __get_fieldname(self, fieldname, longname):
        """Returns a unique field name built from the long name for ambiguous
        field names.

        Used because some sub files use the same variable names.

        Args:
            fieldname (str): field name to replace with a new name, if not
                unique.
            longname (str): longname from which to build a new unique field
                name
        Returns:
            str: a unique field name among all the files in a SAFE container.
        """
        if fieldname in ["SST", "SST_uncertainty", "exception"]:
            return longname.replace(' ', '_')
        return fieldname

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
        if fieldname == 'time':
            return 'time'
        matching = {'lat': 'latitude', 'lon': 'longitude', 'time': 'time',
                    'z': 'elevation'}[fieldname]
        native_fieldname = matching + '_%sn' % self.__sltype
        return native_fieldname

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
        if fieldname == 'time':
            rows = self.get_dimsize('row')
            cols = self.get_dimsize('cell')
            variable = Variable(
                shortname='time',
                description='time of measurement',
                authority=None,
                standardname=None
                )
            field = Field(
                variable,
                OrderedDict([('row', rows), ('cell', cols)]),
                datatype=dtype(int64)
                )
            field.attach_storage(self.get_field_handler(fieldname))
            field.units = self.__time_handler.get_handler().\
                variables['time_stamp_%s' % self.__sltype[0]].units
            return field
        elif fieldname in ['lat', 'lon', 'z']:
            native_name = self.get_geolocation_field(fieldname)
            geofield = self.__geofieldlocator[native_name].read_field(
                native_name
                )
            geofield.name = fieldname
            return geofield
        else:
            native_name = self.__get_native_fieldname(fieldname)
            field = self.__fieldlocator[fieldname].read_field(native_name)
            field.name = fieldname
            field.attach_storage(self.get_field_handler(fieldname))
            dims = field.dimensions
            renamed_dims = OrderedDict()
            for dim in dims:
                newdim = dim
                if dim not in ['row', 'cell', 'time', 'z']:
                    if self.__is_oblique(fieldname):
                        newdim += '_o'
                    else:
                        newdim += '_n'
                renamed_dims[newdim] = dims[dim]
            field.dimensions = renamed_dims
            # for oblique view, the swath is narrower. It is padded with
            # dummy values to stack it over nadir view fields
            if fieldname in self.__oblique_fields:
                if 'cell' in field.get_dimnames():
                    field.dimensions['cell'] = self.get_dimsize('cell')
            return field

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
        native_name = self.__get_native_fieldname(fieldname)
        rowslice = None
        if slices is not None:
            rowslice = [slices[0]]
        if fieldname == 'time':
            suffix = self.__sltype
            SCANSYNC = self.__time_handler.read_values('SCANSYNC')[0]
            PIXSYNC_i = self.__time_handler.read_values(
                'PIXSYNC_%s' % suffix)[0]
            prefix = PREFIX['n']
            first_scan_i\
                = self.__time_handler.read_values(
                    '%s_First_scan_%s' % (prefix, suffix),
                    slices=rowslice)
            first_min_ts\
                = self.__time_handler.read_values(
                    '%s_Minimal_ts_%s' % (prefix, suffix),
                    slices=rowslice)
            scanfield = 'scan_%sn' % suffix
            pixelfield = 'pixel_%sn' % suffix
            indices_handler = self.__fieldlocator[scanfield]
            scan = indices_handler.read_values(scanfield,
                                               slices=slices)
            pixel = indices_handler.read_values(pixelfield,
                                                slices=slices)
            time = first_min_ts.reshape((-1, 1))\
                + (scan - first_scan_i.reshape((-1, 1)))\
                * SCANSYNC + pixel * PIXSYNC_i
            # mask wrong times (which occur in test data)
            maxdate = date2num(self.get_end_time(),
                               "microseconds since 2000-01-01T00:00:00Z")
            mindate = date2num(self.get_start_time(),
                               "microseconds since 2000-01-01T00:00:00Z")

            time = ma.masked_where(
                ((time < mindate) | (time > maxdate)),
                time,
                copy=False
                )
            return time
        elif fieldname in ['lat', 'lon', 'z']:
            return self.__geofieldlocator[native_name].read_values(native_name,
                                                                   slices)
        elif self.__is_oblique(fieldname):
            # oblique views has not the same cell dimension than nadir
            # we want to stack fields from both views by padding fillvalues
            celldim = None
            try:
                dims = self.__fieldlocator[fieldname].get_dimensions(
                    native_name
                    )
                empty = False
                if 'cell' in dims:
                    dimsizes = [self.get_dimsize(dim) for dim in dims]
                    celldim = list(dims).index('cell')
                    nadir_slices = cerbere.mapper.slices.get_nice_slices(
                        slices,
                        dimsizes)
                    sli = nadir_slices[celldim]
                    nad_start, nad_end, step = sli.start, sli.stop, sli.step
                    obliquewidth = (self.__fieldlocator[fieldname]
                                    .get_dimsize('cell'))
                    obl_start, obl_end = nad_start, nad_end
                    obl_start = max(0,
                                    obl_start - self.nadir_to_oblique_offset)
                    if obl_start > obliquewidth:
                        obl_start = obliquewidth
                    obl_end = max(0,
                                  min(obliquewidth,
                                      obl_end - self.nadir_to_oblique_offset)
                                  )
                    newslices = list(nadir_slices)
                    newslices[celldim] = slice(obl_start, obl_end, step)
                    if obl_start >= obl_end:
                        empty = True
                else:
                    # case of some fields such as orphan pixels which don't
                    # have a cell dimension
                    newslices = list(slices) if slices is not None else None
            except ValueError:
                raise

            # read values in oblique view grid
            if celldim is not None and empty:
                values = ma.masked_all(
                    (),
                    dtype=self.__fieldlocator[fieldname]._handler.variables[native_name].dtype)
            else:
                values = self.__fieldlocator[fieldname].read_values(
                    self.__fieldtranslate[fieldname],
                    newslices)
            # padding for missing oblique values to match nadir view grid
            if celldim is not None:
                shape = cerbere.mapper.slices.get_shape_from_slice(
                    cerbere.mapper.slices.get_nice_slices(slices, dimsizes))
                padded_values = ma.masked_all(
                    tuple(shape),
                    dtype=values.dtype,
                    )
                if values.shape == () or min(list(values.shape)) <= 0:
                    # empty result => return padded values only
                    return padded_values
                padded_slice = []
                for dim in dims:
                    if dim != 'cell':
                        padded_slice.append(slice(None, None, None))
                    else:
                        offset = (self.nadir_to_oblique_offset -
                                  nad_start + obl_start)
                        padded_slice.append(slice(
                            offset,
                            offset + values.shape[celldim]
                        ))
                padded_values[padded_slice] = values
                return padded_values
            else:
                return values
        else:
            return self.__fieldlocator[fieldname].read_values(native_name,
                                                              slices)

    def read_fillvalue(self, fieldname):
        """Read the fill value of a field.

        Args:
            fieldname (str): name of the field.

        Returns:
            number or char or str: fill value of the field. The type is the
                as the type of the data in the field.
        """
        return self.__fieldlocator[fieldname].read_fillvalue(fieldname)

    def read_global_attributes(self):
        """Returns the names of the global attributes.

        Returns:
            list<str>: the list of the attribute names.
        """
        # all files seem to have the same list og global attributes.
        return self.__geod_handler.read_global_attributes()

    def read_global_attribute(self, name):
        """Returns the value of a global attribute.

        Args:
            name (str): name of the global attribute.

        Returns:
            str, number or datetime: value of the corresponding attribute.
        """
        # all files seem to have the same list or global attributes.
        return self.__geod_handler.read_global_attribute(name)

    def read_field_attributes(self, fieldname):
        """Return the specific attributes of a field.

        Args:
            fieldname (str): name of the field.

        Returns:
            dict<string, string or number or datetime>: a dictionary where keys
                are the attribute names.
        """
        return self.__fieldlocator[fieldname].read_field_attributes(fieldname)

    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage.

        Returns:
            datetime: start time of the data in file.
        """
        varname = 'time_stamp_%s' % self.__sltype[0]
        vardate = self.__time_handler.get_handler().variables[varname]
        return num2date(vardate[0], vardate.units)

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage.

        Returns:
            datetime: end time of the data in file.
        """
        # WRONG!!!
        varname = 'time_stamp_%s' % self.__sltype[0]
        vardate = self.__time_handler.get_handler().variables[varname]
        return num2date(vardate[-1], vardate.units)

    def get_bbox(self):
        """Returns the bounding box of the feature, as a tuple.

        Returns:
            tuple: bbox expressed as (lonmin, latmin, lonmax, latmax)
        """
        return None

    def write_global_attributes(self, attrs):
        """Write the global attributes of the file.

        Args:
            attrs (dict<string, string or number or datetime>): a dictionary
                containing the attributes names and values to be written.
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

    def write_field(self, fieldname):
        """Writes the field data on disk.

        Args:
            fieldname (str): name of the field to write.
        """
        raise NotImplementedError


class SAFESLIRFile(SAFESLFile):
    """Mapper class for S-3 SLSTR IR files
    """
    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        super(SAFESLIRFile, self).__init__(
            url=url, sltype='i', mode=mode, **kwargs)


class SAFESL500AFile(SAFESLFile):
    """Mapper class for S-3 SLSTR 500m and A-Stripe files
    """
    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        super(SAFESL500AFile, self).__init__(
            url=url, sltype='a', mode=mode, **kwargs)


class SAFESL500BFile(SAFESLFile):
    """Mapper class for S-3 SLSTR 500m and B-Stripe files
    """
    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        super(SAFESL500BFile, self).__init__(
            url=url, sltype='b', mode=mode, **kwargs)


class SAFESL500TDIFile(SAFESLFile):
    """Mapper class for S-3 SLSTR 500m TDI files in nadir view
    """
    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        super(SAFESL500TDIFile, self).__init__(
            url=url, sltype='c', mode=mode, **kwargs)


class SAFESLWSTFile(GHRSSTNCFile):
    """Mapper class for S-3 SLSTR WST files
    """
    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        if url is not None and 'L2P.nc' not in url:
            fullurl = os.path.join(url, 'L2P.nc')
        super(SAFESLWSTFile, self).__init__(
            url=fullurl, mode=mode, **kwargs)

    def read_values(self, fieldname, slices=None):
        values = GHRSSTNCFile.read_values(self, fieldname, slices=slices)
        if fieldname in ['lat', 'lon']:
            values = ma.masked_values(values, -0.000999)
        return values
