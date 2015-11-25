# encoding: utf-8
"""
cerbere.mapper.safeocnncfile
============================

Mapper class for ESA SAFE OCN NetCDF files (for Sentinel-1)

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import os
import logging
import copy
import re
from collections import OrderedDict
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime

import numpy
from netCDF4 import date2num
import gdal
from gdalconst import GA_ReadOnly
from osr import SpatialReference

from .. import READ_ONLY
from .ncfile import NCFile
from ..datamodel.field import Field
from ..datamodel.variable import Variable

# type of product retrieved.
WIND = 'WIND'
WAVE = 'WAVE'
DOPPLER = 'DOPPLER'

# time units
TIME_CONVENTION = 'seconds since 2001-01-01T00:00:00'

class SAFEOCNNCFile(NCFile):
    '''
    Generic storage class for Sentinel-1 OCN files.

    S-1 OCN files are actually netcdf files with specific conventions (they
    are not CF compliant) and additional metadata provided in joined XML
    files (in `annotation` folder)
    '''

    def __init__(self, url=None, mode=READ_ONLY, product=WIND, **kwargs):
        """
        Class constructor

        Kwargs
        ======
        :type product: enum (WIND, WAVE, DOPPLER)
        :keyword product: OCN products contain different fiels that may have different resolution
            and therefore different attached lat/lon.

        Usage
        =====
        OCN products contain different fields at different resolution and
        therefore different attached lat/lon. Fields with different resolution
        are considered as different products and must be opened separately.

        Example :
        -------
        Reading product with all fields at the same resolution :
        >>> f = SAFEOCNNCFile(url='asa-is2-ocn-hh-20051127t235932-20051128t002731-019580-000000-048.nc',
                              product = WIND)

        Reading a field at its native resolution :
        >>> f = SAFEOCNNCFile(url='asa-is2-ocn-hh-20051127t235932-20051128t002731-019580-000000-048.nc',
                              product = WAVE)

        """
        NCFile.__init__(self, url=url, mode=mode, **kwargs)
        self._attributes = OrderedDict()
        if product not in [WIND, WAVE, DOPPLER]:
            raise Exception("Unknown product type (must be WIND, WAVE or DOPPLER) : %s",
                            product)
        if ('IW' in url or 'EW' in url) and product=='WAVE':
            raise Exception("IW and EW product do not contains any WAVE data : %s",
                            product)
#         print "PPPP", product
        logging.debug('product sought %s',product)
        self._product = product
        return

    def close(self):
        self._handler = None
        return

    def get_geolocation_field(self, fieldname):
        """
        return the equivalent field name in the file format for a standard
        geolocation field (lat, lon, time).

        Used for internal purpose and should not be called directly.

        :param fieldname: name of the standard geolocation field (lat, lon
            or time)
        :type fieldname: str

        :rtype: str or None
        :return: name of the corresponding field in the native file format.
            Returns None if no matching is found
        """
        if fieldname == 'time':
            return 'time'
        elif fieldname in ['lat', 'lon']:
            GEO_MATCH = {
                DOPPLER: {'lat': 'rvlLat', 'lon': 'rvlLon'},
                WIND: {'lat': 'owiLat', 'lon': 'owiLon'},
                WAVE: {'lat': 'oswLat', 'lon': 'oswLon'},
                }
            return GEO_MATCH[self._product][fieldname]
        else:
            return None

    def get_matching_dimname(self, dimname):
        """
        Return the equivalent name in the native format for a standard
        dimension.

        This is a translation of the standard names to native ones. It is used
        for internal purpose only and should not be called directly.

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        :arg dimname: standard dimension name
        :type dimname: str

        :rtype: str
        :return: return the native name for the dimension.
            return `dimname` if the input dimension has no standard name.

        *See*

        see :func:`get_standard_dimname` for the reverse operation
        """
        if dimname == 'time':
            return 'time'
        if dimname in ['row', 'cell']:
            DIMS = {
                DOPPLER: {'row': 'rvlAzSize', 'cell': 'rvlRaSize'},
                WIND: {'row': 'owiAzSize', 'cell': 'owiRaSize'},
                WAVE: {'row': 'oswAzSize', 'cell': 'oswRaSize'},
                }
            return DIMS[self._product][dimname]
        else:
            return dimname

    def get_standard_dimname(self, dimname):
        """
        Returns the equivalent standard dimension name for a
        dimension in the native format.

        This is a translation of the native names to standard ones. It is used
        for internal purpose and should not be called directly.

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        :arg dimname: native dimension name
        :type dimname: str

        :rtype: str
        :return: return the (translated) standard name for the dimension.
            return `dimname` if the input dimension has no standard name.

        *See*

        see :func:`get_matching_dimname` for the reverse operation
        """
        STDDIMS = {'rvlAzSize': 'row',
                   'rvlRaSize': 'cell',
                   'owiAzSize': 'row',
                   'owiRaSize': 'cell',
                   'oswAzSize': 'row',
                   'oswRaSize': 'cell',
                   }
        if dimname in STDDIMS:
            return STDDIMS[dimname]
        else:
            return dimname

    def get_fieldnames(self):
        """
        Returns the list of geophysical fields stored for the feature.

        Fields that are not at the same resolution (and therefore have
        different lat/lon) than the opened product (ALL, WIND, WAVE,
        DOPPLER) are not returned. Only consistent fields with the
        `product` specified in :func:`open` call are returned.

        :return: list of field names
        :rtype: list<str>
        """
        allfields = self.get_handler().variables.keys()
        # remove geospatial and temporal fields
        allfields.remove('rvlLon')
        allfields.remove('rvlLat')
        allfields.remove('oswLon')
        allfields.remove('oswLat')
        allfields.remove('owiLon')
        allfields.remove('owiLat')
        # keep only fields related to the opened product
        for field in copy.copy(allfields):
            first_dim = self.get_handler().variables[field].dimensions[0]
            if 'WV' not in self._url:#WV have same grid for the 3 products thus not need to remove fields
                if (self._product == WAVE and not 'osw' in field) or \
                    (self._product == WIND and not 'owi' in field ) or \
                    (self._product == DOPPLER and not 'rvl' in field):
                    allfields.remove(field)
                #add to patch the problem of variable such as owiWindSeaHs(oswAzSize,oswRaSize)
                elif (self._product == WAVE and not 'osw' in first_dim) or \
                    (self._product == WIND and not 'owi' in first_dim) or \
                    (self._product == DOPPLER and not 'rvl' in first_dim):
                    allfields.remove(field)
        return allfields

    def get_dimensions(self, fieldname=None):
        """
        Return the dimension names of a file or a field in the file

        :keyword fieldname: the field from which to get the dimension names.
            For a geolocation field, use the cerbere standard name 
            (time, lat, lon), though native field name will work too.
        :type fieldname: str

        :return: the dimensions of the field or file.
        :rtype: tuple of strings
        """
        if fieldname is None:
            filedims = copy.copy(self.get_handler().dimensions.keys())
            filedims.remove('owiRaSize')
            filedims.remove('owiAzSize')
            filedims.remove('oswRaSize')
            filedims.remove('oswAzSize')
            filedims.remove('rvlRaSize')
            filedims.remove('rvlAzSize')
            filedims.insert(0, 'cell')
            filedims.insert(0, 'row')
            filedims.insert(0, 'time')
            return tuple(filedims)
        else:
            dims = self.get_handler().variables[fieldname].dimensions
            # convert dims to standard names
            newdims = []
            for d in dims:
                newdims.append(self.get_standard_dimname(d))
            return tuple(newdims)

    def read_field(self, fieldname):
        """
        Return the :class:`cerbere.field.Field` object corresponding to
        the requested fieldname.

        The :class:`cerbere.field.Field` class contains all the metadata
        describing a field (equivalent to a variable in netCDF).

        :param fieldname: name of the field
        :type fieldname: str

        :return: the corresponding field object
        :rtype: :class:`cerbere.field.Field`
        """
        # create virtual fields for time
        if fieldname == 'time':
            variable = Variable(
                        shortname='time',
                        description='time of measurement',
                        authority=None,
                        standardname=None
                        )
            field = Field(
                variable,
                OrderedDict([('time', 1)]),
                datatype=numpy.dtype(numpy.float64)
                )
            field.attach_storage(self.get_field_handler(fieldname))
            field.units = TIME_CONVENTION
            return field
        else:
            return super(SAFEOCNNCFile, self).read_field(fieldname)

    def read_values(self, fieldname, slices=None):
        """
        Read the values of a field.

        `slices` is optional. When provided, give for each dimension the 
        corresponding python slice object to subset this dimension. Only the
        dimensions to be subsetted need to be provided, by default the full
        dilmension length is read.

        .. code-block:: python

            # extracting a subset of a field with slices
            data = fd.read_values('owiWindSpeed', slices={'row':slice(10,20), 'cell':slice(30, 40)})

        :param fieldname: name of the field
        :type fieldname: str
        :param slices: dimensions slices, when reading a subset only for some
            dimensions
        :type slices: Dict<str, slice>
        """
        if fieldname == 'time':
            # accounts for virtual variable
            timeval = self.get_start_time()
            return numpy.ma.array([date2num(
                                    timeval,
                                    TIME_CONVENTION)
                                  ])
        else:
            return super(SAFEOCNNCFile, self).read_values(fieldname,
                                                          slices=slices)

    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage"""
        try:
            res = datetime.strptime(
                    self.get_handler().getncattr('firstMeasurementTime'),
                    '%Y-%m-%dT%H:%M:%S.%fZ')
        except:
            res= datetime.strptime(
                    self.get_handler().getncattr('firstMeasurementTime'),
                    '%Y-%m-%dT%H:%M:%SZ')
        return res

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage"""
        try:
            res = datetime.strptime(
                    self.get_handler().getncattr('lastMeasurementTime'),
                    '%Y-%m-%dT%H:%M:%S.%fZ')
        except:
            res = datetime.strptime(
                    self.get_handler().getncattr('lastMeasurementTime'),
                    '%Y-%m-%dT%H:%M:%SZ')
        return res

#     def get_bbox(self):
#         '''
#         returns the bounding box of the feature, as a tuple
#          (lonmin, latmin, lonmax, latmax)
#         '''
#         dimlonname = self.get_geolocation_field('lon')
#         dillatname = self.get_geolocation_field('lat')
#         vlon = self.get_handler().variables[dimlonname][:]
#         vlat = self.get_handler().variables[dillatname][:]
#         if vlon.shape[0] == 1:
#             lonmin = lonmax = vlon[0, 0]
#             latmin = latmax = vlat[0, 0]
#         else:
#             lonmin = min(vlon[0, 0], vlon[0, -1], vlon[-1, -1], vlon[-1, 0])
#             lonmax = max(vlon[0, 0], vlon[0, -1], vlon[-1, -1], vlon[-1, 0])
#             latmin = min(vlat[0, 0], vlat[0, -1], vlat[-1, -1], vlat[-1, 0])
#             latmax = max(vlat[0, 0], vlat[0, -1], vlat[-1, -1], vlat[-1, 0])
#         return (lonmin, latmin, lonmax, latmax)
    
    def get_bbox(self):
        """
        returns the bounding box of the feature, as a tuple
         (lonmin, latmin, lonmax, latmax)
        """
        parent_dur = os.path.abspath(os.path.join(os.path.dirname(self._url),os.pardir))
        manifestpath = os.path.join(parent_dur,'manifest.safe')
        imagette_number = int(os.path.basename(self._url)[-6:-3])
        logging.debug('imagette number %s , filename %s',imagette_number,self._url)
        xmlcontent = minidom.parse(manifestpath)
        hits = xmlcontent.getElementsByTagName('gml:coordinates')
        tmp = hits[imagette_number-1].firstChild.nodeValue
        lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4  = re.split('\s|,',tmp)
        lat1 = float(lat1)
        lon1 = float(lon1)
        lat2 = float(lat2)
        lon2 = float(lon2)
        lat3 = float(lat3)
        lon3 = float(lon3)
        lat4 = float(lat4)
        lon4 = float(lon4)
        lonmax = max(lon1,lon2,lon3,lon4)
        lonmin = min(lon1,lon2,lon3,lon4)
        latmax = max(lat1,lat2,lat3,lat4)
        latmin = min(lat1,lat2,lat3,lat4)
#         logging.debug('lonmin %s, latmin %s, lonmax %s, latmax %s',lonmin, latmin, lonmax, latmax)
        return (lonmin, latmin, lonmax, latmax)