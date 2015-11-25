# -*- coding: utf-8 -*-
"""
.. module::cerbere.mapper.gribfile

Mapper classs for grib format


note: there is already a gribfile.py mapper but I wnated to test the fact that parameter in grib files are splited into many messages
(I suspected that becasue the soft panoply gave me lots of parameters(114) and pygrib (106) messages and the cerbere mapper only 65 parameters but at the end the cerbere mapper seems to be ok)  
"""

import logging
import datetime
import copy
import collections

import numpy
import pygrib
import pdb
from netCDF4 import date2num

from .. import READ_ONLY, WRITE_NEW, DEFAULT_TIME_UNITS
from .abstractmapper import AbstractMapper
from ..datamodel.field import Field
from ..datamodel.variable import Variable


class GribFile(AbstractMapper):
    """Mapper class to read grib files"""

    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        """
        Initialize a Grib file mapper
        """
        super(GribFile, self).__init__(url=url, mode=mode, **kwargs)
        self._messages = collections.OrderedDict()
        self._param_id = collections.OrderedDict()
        return

    def get_geolocation_field(self, fieldname):
        return {'time': 'time', 'lon': 'lon', 'lat': 'lat'}[fieldname]

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

        See:
        see :func:`get_matching_dimname` for the reverse operation
        """
        return dimname

    def get_matching_dimname(self, dimname):
        """
        Return the equivalent name in the current data mapper for a standard
        feature dimension, or None if not found.

        To be derived when creating an inherited data mapper class.
        """
        return dimname
    
    def get_paramid(self,fieldname):
        for gg in self._param_id.keys():
            if self._param_id[gg] == fieldname:
                print 'gg',gg
                res = gg
        return gg
    
    #@profile
    def _open_grib(self):
        """Open a Grib file"""
        url = self._url
        logging.debug("opening in mode : %s", self._mode)
        logging.debug("url : %s", url)
#         print "START", url
        self._handler = pygrib.open(url)
#         print "OPEN"
        # get list of fields
        self._handler.seek(0)
        for mm,msg in enumerate(self._handler):
#             print msg.shortName [msg.paramId
            self._param_id[mm+1] = msg.shortName
            self._messages[mm+1] = msg
#         pdb.set_trace()
#             if msg.shortName in self._messages:
#                 self._messages[msg.shortName].append(msg)
#             else:
#                 self._messages[msg.shortName] = [msg]
        #define an ordered longitudes indices array
#         firstfield = self._messages[self._messages.keys()[0]][0]
        firstfield = self._messages.values()[0]
        values = firstfield.longitudes
        #check whether USA is on the left
        if (values>300).any():
            self._need_flip_longitudes = True
        else:
            self._need_flip_longitudes = False
        return self._handler

    def open(self, view=None, datamodel=None, datamodel_geolocation_dims=None):
        """Open the file (or any other type of storage)

        Args:
            view (dict, optional): a dictionary where keys are dimension names
                and values are slices. A view can be set on a file, meaning
                that only the subset defined by this view will be accessible.
                This view is expressed as any subset (see :func:`get_values`).
                For example:

                view = {'time':slice(0,0), 'lat':slice(200,300),
                'lon':slice(200,300)}

            datamodel (str): type of feature read or written. Internal argument
                only used by the classes from :package:`cerbere.datamodel`
                package. Can be 'Grid', 'Swath', etc...

            datamodel_geolocation_dims (list, optional): list of the name of the
                geolocation dimensions defining the data model to be read in
                the file. Optional argument, only used by the datamodel
                classes, in case the mapper class can store different types of
                data models.
        """
        if self._url:
            handler = self._open_grib()
        else:
            if self._mode == WRITE_NEW:
                raise Exception('Can not open file time series in write mode')
            return self._urlseries
        if datamodel_geolocation_dims:
            for dim in datamodel_geolocation_dims:
                match = self.get_matching_dimname(dim)
                print dim, match
                self._geodim_model2mapper[dim] = match
                self._geodim_mapper2model[match] = dim
        return handler

    def close(self):
        """Close file"""
        self._url = None  # à défaut de suppression
        if self._handler:
            self._handler.close()
            self._handler = None

    def read_values(self, fieldname, slices=None):
        """
        """
        if fieldname == 'lat':
            firstfield = self._messages[self._messages.keys()[0]][0]
            if firstfield.gridType == 'regular_ll':
                values = firstfield.distinctLatitudes
            elif firstfield.gridType == 'regular_gg':
                values = firstfield.latitudes
                nblat = firstfield.Nj
                nblon = firstfield.Ni
                values = numpy.reshape(values,(nblat,nblon))
                values = values[:,0] #all the values are identical so keep only a vector
#                 if self._need_flip_longitudes==True:
#                     values = numpy.roll(values, (values.shape[1] / 2)-1, axis=1)
        elif fieldname == 'lon':
            firstfield = self._messages[self._messages.keys()[0]][0]
            if firstfield.gridType == 'regular_ll':
                values = firstfield.distinctLongitudes
            elif firstfield.gridType == 'regular_gg':
                values = firstfield.longitudes
                nblat = firstfield.Nj
                nblon = firstfield.Ni
                values = numpy.reshape(values,(nblat,nblon))
                values = values[0,:]
            #put the USA before EURASIA in longitudes vector
            if self._need_flip_longitudes==True:
                values = numpy.roll(values, (len(values)/ 2)-1, axis=0)
#                 values = numpy.roll(values, (values.shape[1] / 2)-1, axis=1)
            ind = numpy.ma.where(values >= 180.)
            values[ind] = values[ind] - 360.
        elif fieldname == 'time':
            firstfield = self._messages[self._messages.keys()[0]][0]
            values = numpy.array([date2num(firstfield.analDate,
                                          DEFAULT_TIME_UNITS
                                          )]
                                 )
            return values
        else:
            paramid = self.get_paramid(fieldname)
            fieldmsg = self._messages[paramid][0]
            if len(self._messages[paramid]) > 1:
                values  = fieldmsg.values
                #put the USA before EURASIA in longitudes vector
                if self._need_flip_longitudes==True:
                    values = numpy.roll(values, (values.shape[1] / 2)-1, axis=1)
                    
            else:
                values = fieldmsg.values
                
                values = numpy.roll(values, (values.shape[1] / 2)-1, axis=1)
            values = numpy.ma.array(values)
            values.mask = numpy.zeros(values.shape)
        if slices:
            return values[slices]
        else:
            return values

    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        # Variable
        if not fieldname in ['time', 'lat', 'lon']:
            attrs = self.read_field_attributes(fieldname)
        else:
            attrs = {}
        if 'name' in attrs:
            descr = attrs['name']
        else:
            descr = None
        var = Variable(shortname=fieldname, description=descr)
        dims = collections.OrderedDict()
        if not fieldname in ['time', 'lat', 'lon']:
            if attrs['gridType'] == 'regular_ll' or attrs['gridType'] == 'regular_gg':
                dims['time'] = 1
                dims['y'] = attrs['Nj']
                dims['x'] = attrs['Ni']
            else:
                raise NotImplementedError
            param_id = self.get_paramid(fieldname) 
            datatype = self._messages[param_id][0]['values'].dtype
            field = Field(var, dims, datatype=datatype)
            field.attach_storage(self.get_field_handler(fieldname))
            # MetaData
            field.attributes = {}
            field.units = None
            if 'units' in attrs:
                field.units = attrs['units']
            field.valid_min = None
            field.valid_max = None
            if 'minimum' in attrs:
                field.valid_min = attrs['minimum']
            if 'maximum' in attrs:
                field.valid_max = attrs['maximum']
            return field
        elif fieldname == 'time':
            dims['time'] = 1
            field = Field(var, dims, datatype=numpy.dtype(numpy.float64))
            field.attach_storage(self.get_field_handler(fieldname))
        else:
            firstfield = self._messages[self._param_id.keys()[0]][0]
            if firstfield.gridType == 'regular_ll' or firstfield.gridType == 'regular_gg':
                if fieldname == 'lat':
                    dims['y'] = firstfield.Nj
                    units = 'degrees_north'
                    var.description = 'latitude'
                elif fieldname == 'lon':
                    dims['x'] = firstfield.Ni
                    units = 'degrees_east'
                    var.description = 'longitude'
                field = Field(
                            var,
                            dims,
                            datatype=numpy.dtype(numpy.float64),
                            units=units)

                field.attach_storage(self.get_field_handler(fieldname))
            else:
                raise NotImplementedError
        return field

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

    def get_fieldnames(self):
        """
        """
        return self._param_id.values()

    def get_dimensions(self, fieldname=None):
        """
        """
        raise NotImplementedError

    def read_field_attributes(self, fieldname):
        """
        """
        attrs = {}
        param_id = self.get_paramid(fieldname)
        for k in self._messages[param_id][0].keys():
            if not k in [
                    'values',
                    'distinctLongitudes',
                    'distinctLatitudes',
                    'latitudes',
                    'longitudes',
                    'latLonValues',
                    'codedValues'
                    ]:
                try:
                    val = self._messages[param_id][0][k]
                except:
#                     logging.error('no value for attribute %s of field %s',k,fieldname)
                    val = None
                attrs[k] = val
        return attrs


    def read_global_attributes(self):
        """
        """
        return {}

    def write_global_attributes(self, attrs):
        """
        write the storage (file) global attributes
        """
        raise NotImplementedError

    def read_global_attribute(self, name):
        """
        """
        return None

    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage"""
        param_id = self.get_paramid(self.get_fieldnames()[0])
        return self._messages[param_id][0].analDate

    def get_end_time(self):
        """
        """
        param_id = self.get_paramid(self.get_fieldnames()[0])
        return self._messages[param_id][0].analDate

    def get_product_version(self):
        """return the product version"""
        raise NotImplementedError

    def get_bbox(self):
        """returns the bounding box of the feature

        Return:
        tuple (lonmin, latmin, lonmax, latmax)
        """
        return (-180., -90., 180., 90.)
