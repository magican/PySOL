# -*- coding: utf-8 -*-
"""
cerbere.mapper.ghrsstncfile
===========================

Mapper class for GHRSST netcdf files

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
import datetime
import os
import logging
from dateutil import parser

import numpy
import netCDF4
from collections import OrderedDict

from .. import READ_ONLY
from .ncfile import NCFile

import re


class GHRSSTNCFile(NCFile):
    """
    """
    def __init__(self, url=None, **kwargs):
        super(GHRSSTNCFile, self).__init__(url, **kwargs)
        self.__collection_id = None
        self.__reversed = None
        self.__is_l2p = None

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
        # url needs to be opened first in order to guess default datamodel
        handler = super(GHRSSTNCFile, self).open(view, datamodel,
                                                 datamodel_geolocation_dims)
        if datamodel is None:
            if self.is_l2():
                feature = 'Swath'
            else:
                feature = 'Grid'
            self._feature_type = feature
        return handler

    def __is_reversed(self):
        """Return True if the file has reversed row/cell axes.

        Internal function only. Should no be called by a user.
        """
        if self.__reversed is None:
            if ('NAVO' in self.get_collection_id() and
                'VIIRS' not in self.get_collection_id()) or\
                    self.get_collection_id() == 'USA-RSS-AMSRE-MW-L2-SST':
                self.__reversed = True
            else:
                self.__reversed = False
        return self.__reversed

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
        dims = self.get_handler().dimensions.keys()
        dim_matching = {
            'time': 'time',
            'cell': 'ni',
            'row': 'nj',
            'x': ['lon', 'ni'],
            'y': ['lat', 'nj'],
        }
        if dimname not in dim_matching:
            return dimname
        match = dim_matching[dimname]
        if type(match) is list:
            for dim in match:
                if dim in dims:
                    choice = dim
            match = choice
        if self.__is_reversed():
            return {'nj': 'ni', 'ni': 'nj'}[match]
        else:
            return match

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
        if self.is_l2():
            matching = {
                'time': 'time',
                'nj': 'row',
                'ni': 'cell',
                'lon': 'x',
                'lat': 'y',
            }
        else:
            matching = {
                'time': 'time',
                'ni': 'x',
                'nj': 'y',
                'lon': 'x',
                'lat': 'y',
            }
        if dimname in matching:
            match = matching[dimname]
        else:
            return dimname
        if self.__is_reversed() and match in ['row', 'cell']:
            return {'row': 'cell', 'cell': 'row'}[match]
        return match

    def get_dimensions(self, fieldname=None):
        """
        """
        if self.is_l2():
            if ('channel' in super(GHRSSTNCFile, self)
                    .get_dimensions(fieldname)):
                # for SLSTR WST product
                return ('channel', 'row', 'cell')
            else:
                return ('row', 'cell',)
        else:
            return super(GHRSSTNCFile, self).get_dimensions(fieldname)

    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature

        Excludes the geolocation fields in the returned list
        '''
        geophyvars = super(GHRSSTNCFile, self).get_fieldnames()
        if 'sst_dtime' in geophyvars:
            geophyvars.remove('sst_dtime')
        return geophyvars

    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        field = super(GHRSSTNCFile, self).read_field(fieldname)
        if fieldname == 'time' and\
                'sst_dtime' in self.get_handler().variables:
            field.dimensions = self.get_full_dimensions('sst_dtime')
            field.variable.description = "time of sst"
            if 'comment' in field.attributes:
                field.attributes.pop('comment')
            if 'long_name' in field.attributes:
                field.attributes.pop('long_name')
        return field

    def read_times(self, slices=None):
        """
        Read time values of a file or file series.

        Takes into account the fact that time information is splitted in two
        netcdf variables for a GHRSST file (time and sst_dtime).
        """
        if self._urlseries is not None:
            raise Exception("Not yet implemented")
        elif self._url is not None:
            var = self.get_handler().variables['time']
            vardt = self.read_values('sst_dtime', slices)
            times = vardt[:] + var[0]
            times = netCDF4.num2date(
                times,
                self.get_handler().variables['time'].units
                )
        return times

    def get_collection_id(self):
        """return the identifier of the product collection"""
        if self.__collection_id is None:
            attrs = self.get_handler().ncattrs()
            if 'DSD_entry_id' in attrs:
                self.__collection_id = self.get_handler().DSD_entry_id
            elif 'id' in attrs:
                self.__collection_id = self.get_handler().id
            else:
                pass
        return self.__collection_id

    def is_l2(self):
        """Return True if the product is a L2P"""
        if self.__is_l2p is None:
            self.__is_l2p = ('L2' in self.get_collection_id())
        return self.__is_l2p

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
        try:
            if fieldname == 'time' and \
                    'sst_dtime' in self.get_handler().variables:
                ncvarname = 'sst_dtime'
            else:
                ncvarname = fieldname
#             if fieldname not in ['lat', 'lon']:
#                 # GHRSST products have a time dimension which is not part of
#                 # cerbere data model
#                 if slices is not None and len(slices) == 2 and\
#                         'bnds' not in fieldname:
#                     slices.insert(0, 0)
            if self.__is_reversed():
                # case where row and cells are reversed
                if slices:
                    if len(slices) == 2:
                        newslices = [slices[1], slices[0]]
                    else:
                        newslices = slices
                    newslices = tuple(newslices)
                else:
                    newslices = slices
                values = super(GHRSSTNCFile, self).read_values(ncvarname,
                                                               newslices)
                if len(values.shape) == 2:
                    values = numpy.ma.transpose(values)
            elif (fieldname == 'lon' or fieldname == 'lat') and \
                    len(self.get_handler().variables[fieldname].shape) == 3:
                # case where the lat/lon have a time dimension (ARC ATSR files)
                if slices is not None:
                    new_slices = slices
                    new_slices.insert(0, 0)
                    values = super(GHRSSTNCFile, self).read_values(ncvarname,
                                                                   tuple(new_slices))
                else:
                    values = super(GHRSSTNCFile, self).read_values(ncvarname,
                                                                   slices)
                    if values.ndim == 3:
                        values = values[0, :, :]
            else:
                values = super(GHRSSTNCFile, self).read_values(
                    ncvarname,
                    slices)
            if self.__is_reversed()\
                    and (fieldname == 'lon' or fieldname == 'lat'):
                # NAVO lat/lon contains 0 values
                if 'NAVO' in self.get_collection_id():
                    values = numpy.ma.masked_equal(values, 0, atol=0.0001,
                                                   copy=False)
            # homogenize longitudes
            if fieldname == 'lon':
                ind = numpy.ma.where(values >= 180.)
                values[ind] = values[ind] - 360.
                ind = numpy.ma.where(values < -180.)
                values[ind] = values[ind] + 360.
            if ncvarname == 'sst_dtime':
                # pixel time is time + sst_dtime
                var = self.get_handler().variables['time']
                values = var[0] + values
            # fix invalid values
            if (self.get_collection_id() == 'JPL-L2P-MODIS_A' or
                    self.get_collection_id() == 'JPL-L2P-MODIS_T'):
                # some lats have nan
                if fieldname == 'lat':
                    values = numpy.ma.fix_invalid(values, copy=False)
                # and corresponding lons have 0.0.... !
                elif fieldname == 'lon':
                    lats = self.read_values('lat', slices=slices)
                    lats = numpy.ma.fix_invalid(lats, copy=False)
                    values = numpy.ma.array(values, mask=lats.mask,
                                            copy=False)
            return values
        except:
            logging.error("Could not read values for fieldname: {}".format(
                fieldname
            ))
            raise

    def get_cycle_number(self):
        return None

    def get_orbit_number(self):
        if self.get_collection_id() == 'USA-RSS-AMSRE-MW-L2-SST':
            fname = os.path.basename(self.get_url())
            orbit = int(fname.split('_r')[1].split('.dat')[0])
            return orbit

    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage"""
        handler = self.get_handler()
        attrs = handler.ncattrs()
        if 'time_coverage_start' in attrs:
            attrdate = handler.getncattr('time_coverage_start')
            # case of Pathfinder - crappy times
            if 'T24' in attrdate:
                # not sure this is what we should do here
                logging.warning("Strange start time %s", attrdate)
                attrdate = attrdate.replace('T24', 'T00')
            return parser.parse(attrdate)
#         if "arc-upa-" in self.get_collection_id().lower():
#             start_time = handler.getncattr('time_coverage_start')
#             return datetime.datetime.strptime(
#                 start_time, "%Y-%m-%d %H:%M:%SZ"
#             )

        elif 'start_date' in attrs:
            attrdate = handler.getncattr('start_date').replace(' UTC', '')
            if 'start_time' in attrs:
                attrtime = handler.getncattr('start_time')
                attrdate = attrdate + 'T' + attrtime.replace(' UTC', '')
            if '.' in attrdate:
                return datetime.datetime.strptime(
                    attrdate, "%Y-%m-%dT%H:%M:%S.%f"
                    )

            else:
                return datetime.datetime.strptime(
                    attrdate, "%Y-%m-%dT%H:%M:%S"
                    )

        elif "start_time" in attrs:
            attrdate = handler.getncattr('start_time')
            if re.match(r"""^\d{8}T\d{6}Z$""", attrdate):
                return datetime.datetime.strptime(
                    attrdate, "%Y%m%dT%H%M%SZ"
                )
        else:
            pass

        return None

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage"""
        handler = self.get_handler()
        attrs = self.get_handler().ncattrs()
        if 'time_coverage_end' in attrs:
            attrdate = handler.getncattr('time_coverage_end')
            # case of Pathfinder - crappy times
            if 'T24' in attrdate:
                # not sure this is what we should do here
                logging.warning("Strange end time %s", attrdate)
                attrdate = attrdate.replace('T24', 'T00')
            return parser.parse(attrdate)
#         if "arc-upa-" in self.get_collection_id().lower():
#             end_time = handler.getncattr('time_coverage_end')
#             return datetime.datetime.strptime(
#                 end_time, "%Y-%m-%d %H:%M:%SZ"
#             )
        elif 'stop_date' in attrs:
            attrdate = handler.getncattr('stop_date').replace(' UTC', '')
            if 'stop_time' in attrs:
                attrtime = handler.getncattr('stop_time')
                attrdate = attrdate + 'T' + attrtime.replace(' UTC', '')
            if '.' in attrdate:
                return datetime.datetime.strptime(
                    attrdate, "%Y-%m-%dT%H:%M:%S.%f"
                    )
            elif re.match(r"""^\d{8}T\d{6}Z$""", attrdate):
                return datetime.datetime.strptime(
                    attrdate, "%Y%m%dT%H%M%SZ"
                )
            else:
                return datetime.datetime.strptime(
                    attrdate, "%Y-%m-%dT%H:%M:%S"
                    )
        elif "stop_time" in attrs:
            attrdate = handler.getncattr('stop_time')
            if re.match(r"""^\d{8}T\d{6}Z$""", attrdate):
                return datetime.datetime.strptime(
                    attrdate, "%Y%m%dT%H%M%SZ"
                )
        else:
            pass

        return None
