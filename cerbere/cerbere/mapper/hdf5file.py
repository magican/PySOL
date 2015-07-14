# encoding: utf-8
"""
cerbere.mapper.hdffile
======================

Mapper class for HDF files

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""


import os
import logging
import copy
from collections import OrderedDict

import h5py
import numpy

from .. import READ_ONLY
from .abstractmapper import AbstractMapper
from ..datamodel.field import Field
from ..datamodel.variable import Variable
import cerbere.mapper.slices

MODES = {READ_ONLY: 'r'}


class HDF5File(AbstractMapper):
    """Generic mapper class for HDF5 files"""

    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        """
        """
        AbstractMapper.__init__(self, url=url, mode=mode, **kwargs)
        self.fieldnames = None
        return

    def open(self,
             view=None,
             datamodel=None,
             datamodel_geolocation_dims=None):
        """Open the HDF5 file
        
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
        self.view = view
        if self.is_writable():
            raise NotImplementedError
        else:
            if not os.path.exists(self._url) and self._mode == READ_ONLY:
                raise Exception("File %s is not existing" % self._url)

        if (self._url is not None) and (self._mode is not None):
            logging.debug("MODE : %s", self._mode)
            self._handler = h5py.File(self._url, MODES[self._mode])
            return self._handler
        else:
            return None

    def close(self):
        self.fieldnames = None
        self.get_handler().close()
        self._handler = None
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
        raise NotImplementedError

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
        raise NotImplementedError

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
        raise NotImplementedError

    def __get_group_fieldnames(self, item):
        """Iterate recursively all groups to collect HDF5 datasets (=fields)"""
        fields = {}
        if not isinstance(item, h5py.Dataset):
            fields[item.name] = [
                self.__get_group_fieldnames(item[key]) for key in item.keys()
                ]
            return fields
        else:
            name = item.name
            self.fieldnames.append(name)
            return name

    def get_fieldnames(self):
        """Returns the list of geophysical fields stored for the feature.

        The geolocation field names are excluded from this list.

        Returns:
            list<string>: list of field names
        """
        if self.fieldnames is not None:
            return self.fieldnames
        self.fieldnames = []
        self.__get_group_fieldnames(self.get_handler())
        return self.fieldnames

    def read_field_attributes(self, fieldname):
        """Return the specific attributes of a field.

        Args:
            fieldname (str): name of the field.

        Returns:
            dict<string, string or number or datetime>: a dictionary where keys
                are the attribute names.
        """
        native_fieldname = self.get_geolocation_field(fieldname)
        # collect group attributes
        groups = native_fieldname.strip('/').split('/')
        attrs = OrderedDict([])
        node = self.get_handler()
        for group in groups:
            node = node[group]
            for att, val in node.attrs.items():
                attrs[att] = val
        return attrs

    def get_dimsize(self, dimname):
        """Return the size of a dimension.

        Args:
            dimname (str): name of the dimension.

        Returns:
            int: size of the dimension.
        """
        dim = self.get_matching_dimname(dimname)
        if dim is None:
            dim = dimname
        rank = int(dim[-1])
        fieldname = dim[:-1]
        return self.get_handler()[fieldname].shape[rank]

    def _get_native_fieldname(self, fieldname):
        """Get field native name with full group path"""
        if fieldname in ['lat', 'lon', 'time', 'z']:
            native_fieldname = self.get_geolocation_field(fieldname)
        else:
            native_fieldname = fieldname
        return native_fieldname

    def _get_field_dimensions(self, fieldname):
        """get field dimensions"""
        native_fieldname = self._get_native_fieldname(fieldname)
        # get dimension names. No way at this point to guess which
        # ones are standard. This has to be done in inherited class.
        dims = []
        for i, dim in enumerate(
                self.get_handler()[native_fieldname].dims):
            if dim.label == '':
                # create dimension name with group and dim rank
                dims.append(fieldname + '%d' % i)
            else:
                dims.append(dim.label)
        return tuple(dims)

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
        if fieldname is None:
            dims = []
            for item in self.get_fieldnames():
                fielddims = self._get_field_dimensions(item)
                dims.extend(list(fielddims))
            return tuple(dims)
        else:
            return self._get_field_dimensions(fieldname)

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
        # get field native name
        native_fieldname = self._get_native_fieldname(fieldname)
        # attributes
        varattrs = copy.copy(self.read_field_attributes(native_fieldname))
        variable = Variable(shortname=fieldname)
        # dimensions
        dims = self.get_full_dimensions(native_fieldname)
        datatype = self.get_handler().get(native_fieldname).dtype
        field = Field(variable,
                      dimensions=dims,
                      datatype=datatype,
                      attributes=varattrs
                      )
        field.attach_storage(self.get_field_handler(fieldname))
        # MetaData
        field.units = None
        field.valid_min = None
        field.valid_max = None
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
        native_fieldname = self._get_native_fieldname(fieldname)
        field = self.read_field(fieldname)
        if slices is None:
            var = self.get_handler().get(native_fieldname)
            values = var[:]
        else:
            dims = self.get_full_dimensions(fieldname)
            newslices = cerbere.mapper.slices.get_absolute_slices(
                view=self.view,
                slices=slices,
                dimnames=dims.keys(),
                dimsizes=dims.values()
                )
            values = self.get_handler().get(native_fieldname)[tuple(newslices)]
        fillvalue = field.fillvalue
        if fillvalue is not None:
            values = numpy.ma.array(values,
                                    fill_value=fillvalue)
        return values

    def read_global_attributes(self):
        """Returns the names of the global attributes.

        Returns:
            list<str>: the list of the attribute names.
        """
        return self.get_handler().attrs.keys()

    def get_product_version(self):
        """return the product version"""
        raise NotImplementedError

    def get_orbit_number(self):
        """In the case of a satellite orbit file, returns the orbit number.

        Returns:
            int: the orbit number
        """
        raise NotImplementedError

    def get_cycle_number(self):
        """In the case of a satellite orbit file, returns the cycle number.

        Returns:
            int: the cycle number
        """
        raise NotImplementedError

    def read_global_attribute(self, attr):
        """Returns the value of a global attribute.

        Args:
            name (str): name of the global attribute.

        Returns:
            str, number or datetime: value of the corresponding attribute.
        """
        return self.get_handler().attrs[attr]

    def write_field(self, fieldname):
        """
        """
        raise NotImplementedError

    def read_fillvalue(self, fieldname):
        """
        """
        raise NotImplementedError

    def create_field(self, field, dim_translation=None):
        """
        """
        raise NotImplementedError

    def create_dim(self, dimname, size=None):
        """
        """
        raise NotImplementedError

    def write_global_attributes(self, attrs):
        """
        write the storage (file) global attributes
        """
        raise NotImplementedError

    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage.

        Returns:
            datetime: start time of the data in file.
        """
        raise NotImplementedError

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage.

        Returns:
            datetime: end time of the data in file.
        """
        raise NotImplementedError

    def get_bbox(self):
        '''
        returns the bounding box of the feature, as a tuple
         (lonmin, latmin, lonmax, latmax)
        '''
        return None

    def get_spatial_resolution_in_deg(self):
        """Returns the average spatial resolution in degrees"""
        return None
