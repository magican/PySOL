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

import numpy
from pyhdf.SD import SD, SDC
from pyhdf.HDF import HDF, HC
from pyhdf.V import *

from .. import READ_ONLY
from .abstractmapper import AbstractMapper
from ..datamodel.field import Field
from ..datamodel.variable import Variable


MODES = {'r':SDC.READ, 'w':SDC.WRITE}


class HDFFile(AbstractMapper):
    '''
    Generic storage class for HDF files
    '''

    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        """
        """
        AbstractMapper.__init__(self, url=url, mode=mode, **kwargs)
        return

    def open(self,
             view=None,
             datamodel=None,
             datamodel_geolocation_dims=None):
        """Open the HDF file

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
        self.view=view
        if self.is_writable():
            raise NotImplementedError
        else:
            if not os.path.exists(self._url):
                raise Exception("File %s is not existing" % self._url)

        if (self._url is not None) and (self._mode is not None):
            logging.debug("MODE : %s", self._mode)
            self._handler = SD(self._url, MODES[self._mode])
            # case of vgroup containing some information
            if self._mode == 'r':
                # open HDF file
                self._hdffile = HDF(self._url, HC.READ)
                # initialize V interface on HDF file
                self._vdata = self._hdffile.vstart()
            return self._handler
        else:
            return None

    def close(self):
        self._vdata.end()                    # terminate V interface
        self._hdffile.close()
        self._handler = None
        self._vdata = None
        self._hdffile = None
        return

    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature
        '''
        fields = self.get_handler().datasets().keys()
        # remove here time/space information to keep only geophysical fields
        for field in ['time', 'lat', 'lon']:
            if field in fields:
                fields.remove(self.get_geolocation_field(field))
        return fields

    def read_field_attributes(self, fieldname):
        """
        return the specific storage attributes of a variable
        (_FillValue, scale_factor, add_offset)
        """
        native_fieldname = self.get_geolocation_field(fieldname)
        if native_fieldname is None:
            native_fieldname = fieldname
        attrs = self.get_handler().select(native_fieldname).attributes()
        return attrs

    def get_dimsize(self, dimname):
        hdfdim = self.get_matching_dimname(dimname)
        if hdfdim is None:
            hdfdim = dimname
        for fieldname in self.get_handler().datasets():
            dims = self.get_handler().select(fieldname).dimensions()
            for dim in dims:
                if dim == hdfdim:
                    return dims[dim]
        return None

    def get_dimensions(self, fieldname=None):
        """
        Return the standard dimension names of a file or a field in the file

        :keyword fieldname: the field from which to get the dimension names.
            For a geolocation field, use the cerbere standard name
            (time, lat, lon), though native field name will work too.
        :type fieldname: str

        :return: the standard dimensions of the field or file.
        :rtype: tuple of strings
        """
        if fieldname is None:
            raise NotImplementedError
        else:
            native_fieldname = self.get_geolocation_field(fieldname)
            if native_fieldname is None:
                native_fieldname = fieldname
            var = self.get_handler().select(native_fieldname)
            if var is None:
                raise Exception("Variable %s not existing in file"\
                                    % native_fieldname)
            dims = OrderedDict(
                    sorted(var.dimensions(full=True).items(),
                           key=lambda t: t[1][1])
                    )
            dims = [self.get_standard_dimname(dim) for dim in dims]
            return tuple(dims)

    def read_field(self, fieldname):
        namingauth = None
        native_fieldname = self.get_geolocation_field(fieldname)
        if native_fieldname is None:
            native_fieldname = fieldname
        varattrs = copy.copy(self.read_field_attributes(fieldname))
        if 'long_name' in varattrs:
            descr = varattrs['long_name']
        else:
            descr = None
        variable = Variable(
                        shortname=fieldname,
                        description=descr,
                        authority=namingauth,
                        standardname=None
                        )
        dims = self.get_full_dimensions(fieldname)
        TYPE_CONVERT = {'4': numpy.dtype(numpy.int8),
                        '5': numpy.dtype(numpy.float32),
                        '20': numpy.dtype(numpy.int8),
                        '21': numpy.dtype(numpy.uint8),
                        '22': numpy.dtype(numpy.int16),
                        '23': numpy.dtype(numpy.uint16),
                        '24': numpy.dtype(numpy.int32)
                        }
        typestr = self.get_handler().select(native_fieldname).info()[3]
        rec = Field(
                variable,
                dims,
                datatype=TYPE_CONVERT[str(typestr)]
                )
        rec.attach_storage(self.get_field_handler(fieldname))
        # MetaData
        rec.units = None
        if 'units' in varattrs:
            rec.units = varattrs['units']
        rec.valid_min = None
        rec.valid_max = None
        rec.attributes = {}
        if ('valid_min' in varattrs and 'valid_max' in varattrs)\
                 or 'valid_range' in varattrs:
            if 'valid_range' in varattrs:
                rec.valid_min, rec.valid_max = varattrs['valid_range']
            else:
                rec.valid_min = varattrs['valid_min']
                rec.valid_max = varattrs['valid_max']
            if 'scale_factor' in varattrs:
                rec.valid_min = rec.valid_min * varattrs['scale_factor']
                rec.valid_max = rec.valid_max * varattrs['scale_factor']
            if 'add_offset' in varattrs:
                rec.valid_min = rec.valid_min + varattrs['add_offset']
                rec.valid_max = rec.valid_max + varattrs['add_offset']
        for att in varattrs:
            if not att in ['units', 'scale_factor', 'add_offset',
                           '_FillValue', 'valid_min', 'valid_max',
                           'scale_factor_err', 'add_offset_err',
                           'valid_range', 'calibrated_nt', 'SDS_type',
                           'long_name', 'bad_value_scaled',
                           'bad_value_unscaled']:
                rec.attributes[att] = varattrs[att]
        return rec

    def read_values(self, fieldname, slices=None):
        native_fieldname = self.get_geolocation_field(fieldname)
        if native_fieldname is None:
            native_fieldname = fieldname
        var = self.get_handler().select(native_fieldname)
        if slices is None:
            values = var.get()
        else:
            dims = self.get_full_dimensions(fieldname).keys()
            newslices = []
            # fill in slices with None values
            for ind, slc in enumerate(slices):
                i0, i1, step = slc.start, slc.stop, slc.step
                if i0 is None:
                    i0 = 0
                if i1 is None:
                    i1 = self.get_dimsize(dims[ind])
                if step is None:
                    step = 1
                newslices.append(slice(i0, i1, step))
            # Added conversion to int as get does not support long values.
            slstart = [int(s.start) for s in newslices]
            slstop = [int(s.stop - s.start) for s in newslices]
            slstride = [int(s.step) for s in newslices]
            values = var.get(start=tuple(slstart),
                             count=tuple(slstop),
                             stride=tuple(slstride))
        attrs = self.read_field_attributes(fieldname)
        if '_FillValue' in attrs:
            fill_value = attrs['_FillValue']
        else:
            fill_value = None
        if not fill_value is None:
            values = numpy.ma.array(values, fill_value=fill_value)
        else:
            values = numpy.ma.array(values)
        if 'scale_factor' in attrs:
            values = values * attrs['scale_factor']
        if 'add_offset' in attrs:
            values = values + attrs['add_offset']
        return values

    def read_global_attributes(self):
        return self.get_handler().attributes()

    def read_global_attribute(self, attr):
        """
        """
        return self.read_global_attributes()[attr]

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
        """Returns the minimum date of the file temporal coverage"""
        raise NotImplementedError

    def get_end_time(self):
        """
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
