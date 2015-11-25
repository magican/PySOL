# -*- coding: utf-8 -*-
"""
cerbere.mapper.cersaticessmincfile
============================

Mapper class for CERSAT ice concentration SSM/I netcdf files

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor::  <agrouaze@ifremer.fr>
.. codeauthor::  <agrouaze@ifremer.fr>
"""
import datetime
import os
import logging
from dateutil import parser
import calendar
import numpy
import netCDF4
from collections import OrderedDict
from ..datamodel.field import Field
from ..datamodel.variable import Variable
from .. import READ_ONLY
from .ncfile import NCFile
import pdb
import sys
sys.path.append('/home/losafe/users/agrouaze/git/pathfinders_acidification/src/lib')
from read_dat_file import read_official_nsidc_grid_north
import re

def GetPolarGrid(pole,original=False):
    '''
    return lons and lats or NSIDC grid arctic/antarctic at 25km res
    original (bool) True=>get the full array 25km
    pole (str) north or south
    note: originally we used a sterepolar grid but it was not the official NSIDC one (shift of 7km), now the official grid is retrieved from    
    psn25lats_v3.dat: 304 columns x 448 rows, range = [31.102, 89.836]
    '''
#     file = '/home/cercache/users/agrouaze/temporaire/grid_'+pole+'_12km.nc'
#     nc = netCDF4.Dataset(file,'r')
#     tlon = nc.variables['longitude'][:]
#     tlat = nc.variables['latitude'][:]
#     
#     nc.close()
#     if original:
#         lons = tlon
#         lats = tlat
#     else:
#         
    lats,lons = read_official_nsidc_grid_north()
#     lons = lons[slice(0,lons.shape[0],2),slice(0,lons.shape[1],2)]
#     lats = lats[slice(0,lats.shape[0],2),slice(0,lats.shape[1],2)]
#     lons[lons>180] = n[tlon>180]-360.
    varlon = Variable(
                        shortname='lon',
                        description='longitude',
                        standardname='longitude'
                        )
    varlat = Variable(
                        shortname='lat',
                        description='latitude',
                        standardname='latitude'
                        )
    logging.info('GetPolarGrid | lon: %s lat:%s',lons.shape,lats.shape)
#     dimensions = ['cell','row']
    dimensions = ['y','x']
#     dimensions = ['lat','lon']
    dims = OrderedDict(
            zip(dimensions,
                list(lons.shape)
                )
            )
    fieldlon = Field(
                    varlon,
                    dims,
                    values=lons
                    )
    fieldlat = Field(
                    varlat,
                    dims,
                    values=lats
                    )
    fieldlon.units = 'degrees_east'
    fieldlat.units = 'degrees_north'
    return fieldlon,fieldlat

class CERSATIceSSMINCFile(NCFile):
#     def __init__(self, url=None,pole):
#         super(CERSATIceSSMINCFile, self).__init__(url=url)
#         self._pole = pole
#         return
    
    def get_fieldnames(self):
        fieldnames = super(CERSATIceSSMINCFile, self).get_fieldnames()
        fieldnames.extend(["lon", "lat"])
        return fieldnames
    
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

        dim_matching = {
                    'time': 'time',
                    'lon': 'ni',
                    'lat': 'nj',
                    'row':'nj',
                    'cell':'ni',
                    }
        if geodimname in dim_matching:
            res= dim_matching[geodimname]
        else:
            res = geodimname
        return res

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
        dim_matching = {
                'time': 'time',
                'ni': 'lon',
                'nj': 'lat',
                'nj':'row',
                'ni':'cell',
                }
        if geodimname in dim_matching:
            res= dim_matching[geodimname]
        else:
            res = geodimname

        return res
    
    def get_geolocation_field(self, fieldname):
#         allfields = self.get_fieldnames()
#         if fieldname in allfields:
#             res = fieldname
#         else:
#             res = None
        return fieldname
    

    def read_field(self, fieldname):
        
        if fieldname == 'lon':
            pole = self.read_global_attribute('pole')
            fieldlon,fieldlat = GetPolarGrid(pole,original=True)
            field = fieldlon
        elif fieldname == 'lat':
            pole = self.read_global_attribute('pole')
            fieldlon,fieldlat = GetPolarGrid(pole,original=True)
            field = fieldlat
        elif fieldname == 'concentration':
            srcfield = super(CERSATIceSSMINCFile, self).read_field(fieldname)
            field = srcfield.clone()
            #meta info field: unit,fillvalue,min,max,comments
            meta = ('1', -128., 0., 1., [])#copied on osisaf metadata
            field.set_metadata(meta)
        else:
            field = super(CERSATIceSSMINCFile, self).read_field(fieldname)
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
        pole = self.read_global_attribute('pole')
        fieldlon,fieldlat = GetPolarGrid(pole,original=True)
        if fieldname == 'longitude':
            values = fieldlon.values
        elif fieldname == 'latitude':
            values = fieldlat.values
        elif fieldname == 'concentration':
            values = super(CERSATIceSSMINCFile, self).read_values(fieldname,slices=slices)
            values = values.astype(numpy.float32)
        else:
            values = super(CERSATIceSSMINCFile, self).read_values(fieldname,slices=slices)
        return values
            
    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage"""
        time = self.read_values('time')
        return time

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage"""
        time = self.read_values('time')[0]
        attrs = self.read_field_attributes('time')
        unit = attrs['units']
        time_dt = netCDF4.num2date(time,unit)
        year = time_dt.year
        month = time_dt.month
        nbdays = calendar.monthrange(year, month)[1]
        res = time_dt + datetime.timedelta(days=nbdays-1)
#         res = res - datetime.timedelta(days=1)
        return res
