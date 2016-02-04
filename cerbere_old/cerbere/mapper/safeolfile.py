# -*- coding: utf-8 -*-
"""
.. module::cerbere.mapper.safeslfile

Class to read Sentinel-3 OLCI files.

:copyright: Copyright 2015 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
import os
import glob
from numpy import ma, dtype, int64
from collections import OrderedDict

from netCDF4 import num2date

from .. import READ_ONLY
from .abstractmapper import AbstractMapper
from cerbere.datamodel.variable import Variable
from cerbere.datamodel.field import Field
from cerbere.mapper.ncfile import NCFile

DATAFILES = {
    "L1B": ["*_radiance.nc"],
    "L2LAND": ["iwv.nc", "lqsf.nc", "ogvi.nc", "otci.nc", "rc_ogvi.nc"],
    "L2WATER": ["chl_nn.nc", "chl_oc4me.nc", "iop_nn.nc", "*_reflectance.nc",
                "par.nc", "trsp.nc", "tsm_nn.nc", "w_aer.nc", "wqsf.nc"]
    }


class SAFEOLFile(AbstractMapper):
    """Abstract class for SAFE OLCI files.

    This mapper concatenates together the files within a SAFE folder that
    share the same dimensions.

    url: the path to the product SAFE folder

    L1B
    geo_coordinates.nc, instrument_data.nc, *_radiance.nc, time_coordinates.nc

    L2 LAND
    geo_coordinates.nc, instrument_data.nc, iwv.nc, lqsf.nc, ogvi.nc, otci.nc,
    rc_ogvi.nc, time_coordinates.nc

    L2 WATER
    geo_coordinates.nc, instrument_data.nc, chl_nn, chl_oc4me.nc, iop_nn.nc,
    *_reflectance.nc, par.nc, trsp.nc, tsm_nn, w_aer.nc, wqsf.nc,
    time_coordinates.nc
    """
    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        if mode != READ_ONLY:
            raise Exception("This mapper can only be used in read_only mode.")
        super(SAFEOLFile, self).__init__(url=url, mode=mode, **kwargs)
        self.__data_handlers = []
        self.__fields = {}
        # ancillary files
        cartesian = "geo_coordinates.nc"
        instrument = "instrument_data.nc"
        times = "time_coordinates.nc"

        # handlers for ancillary fields
        self.__time_handler = NCFile(os.path.join(url, times),
                                     mode=mode, **kwargs)
        self.__fields[times] = []
        self.__coord_handler = NCFile(os.path.join(url, cartesian),
                                      mode=mode, **kwargs)
        self.__fields[cartesian] = []
        self.__instr_handler = NCFile(os.path.join(url, instrument),
                                      mode=mode, **kwargs)
        self.__fields[instrument] = []
        # get product type
        safefolder = os.path.basename(url)
        if "_OL_1_ERR" in safefolder:
            datafiles = DATAFILES["L1B"]
        elif "OL_2_LRR" in safefolder:
            datafiles = DATAFILES["L2LAND"]
        elif "OL_2_WRR" in safefolder:
            datafiles = DATAFILES["L2WATER"]
        # detect the data files and instanciate mappers for each one
        for f in datafiles:
            if '*' in f:
                fnames = glob.glob(os.path.join(url, f))
                for fname in fnames:
                    self.__data_handlers.append(NCFile(url=fname, mode=mode,
                                                       **kwargs))
                    self.__fields[f] = []
            else:
                fname = os.path.join(url, f)
                self.__data_handlers.append(NCFile(url=fname, mode=mode,
                                                   **kwargs))
                self.__fields[f] = []
        self.__fieldlocator = {}
        self.__geofieldlocator = {}

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

            datamodel_geolocation_dims (list, optional): list of the name of the
                geolocation dimensions defining the data model to be read in
                the file. Optional argument, only used by the datamodel
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
        for hdlr in self.__data_handlers:
            hdlr.open(view, datamodel,
                      datamodel_geolocation_dims)
        self.__coord_handler.open(view, datamodel,
                                  datamodel_geolocation_dims)
        self.__time_handler.open(rowview, datamodel,
                                 datamodel_geolocation_dims)
        self.__instr_handler.open(view, datamodel,
                                  datamodel_geolocation_dims)
        # build the two-way dictionaries of fields
        # ...for data
        for hdlr in self.__data_handlers:
            self.__fields[os.path.basename(hdlr.get_url())]\
                = hdlr.get_fieldnames()
            for fieldname in hdlr.get_fieldnames():
                self.__fieldlocator[fieldname] = hdlr

    def close(self):
        """Close handler on storage"""
        for hdlr in self.__data_handlers:
            hdlr.close()
        self.__data_handlers = None
        self.__coord_handler.close()
        self.__time_handler.close()
        self.__instr_handler.close()
        self.__coord_handler = None
        self.__time_handler = None
        self.__instr_handler = None

    def get_dimsize(self, dimname):
        """Return the size of a dimension.

        Args:
            dimname (str): name of the dimension.

        Returns:
            int: size of the dimension.
        """
        dim = self.get_matching_dimname(dimname)
        return self.__coord_handler.get_dimsize(dim)

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
            return self.__coord_handler.get_dimensions()
        if fieldname in ['time', 'lat', 'lon']:
            # Should all have the same dimension as lat
            native_fieldname = self.get_geolocation_field('lat')
            dims = self.__coord_handler.get_dimensions(native_fieldname)
        else:
            handler = self.__fieldlocator[fieldname]
            dims = handler.get_dimensions(fieldname)
        # convert geolocation dims to standard names
        newdims = []
        for dim in list(dims):
            newdims.append(self.get_standard_dimname(dim))
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
        matching = {'time': 'time', 'row': 'rows', 'cell': 'columns'}
        if dimname in matching:
            return matching[dimname]
        return dimname

    def get_standard_dimname(self, dimname):
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
        matching = {'time': 'time', 'rows': 'row', 'columns': 'cell'}
        if dimname in matching:
            return matching[dimname]
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
        MATCHES = {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'}
        if fieldname in MATCHES:
            return MATCHES[fieldname]
        return None

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
                variables['time_stamp'].units
            return field
        elif fieldname in ['lat', 'lon']:
            native_name = self.get_geolocation_field(fieldname)
            geofield = self.__coord_handler.read_field(
                native_name
                )
            geofield.name = fieldname
            geofield.attach_storage(self.get_field_handler(fieldname))
            return geofield
        else:
            native_name = self.__get_native_fieldname(fieldname)
            return self.__fieldlocator[native_name].read_field(native_name)

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
        if fieldname == 'time':
            if slices is not None:
                tslices = [slices[0]]
            else:
                tslices = slices
            time = self.__time_handler.read_values('time_stamp',
                                                   slices=tslices)
            # reshape as a 2D field
            rows = self.get_dimsize('row')
            cols = self.get_dimsize('cell')
            if slices is None:
                shape = (cols, rows)
            else:
                newslices = self._fill_slices(slices, (rows, cols))
                shape = (newslices[1].stop - newslices[1].start,
                         newslices[0].stop - newslices[0].start)
            time = ma.resize(time, shape).transpose()
            return time
        elif fieldname in ['lat', 'lon']:
            return self.__coord_handler.read_values(native_name,
                                                    slices)
        else:
            return self.__fieldlocator[native_name].read_values(native_name,
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
        return self.__coord_handler.read_global_attributes()

    def read_global_attribute(self, name):
        """Returns the value of a global attribute.

        Args:
            name (str): name of the global attribute.

        Returns:
            str, number or datetime: value of the corresponding attribute.
        """
        # all files seem to have the same list or global attributes.
        return self.__coord_handler.read_global_attribute(name)

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
        varname = 'time_stamp'
        vardate = self.__time_handler.get_handler().variables[varname]
        return num2date(vardate[0], vardate.units)

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage.

        Returns:
            datetime: end time of the data in file.
        """
        # WRONG!!!
        varname = 'time_stamp'
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
