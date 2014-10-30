#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 12:03:29 2012

@author: mag
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime

from numpy import meshgrid, arange, ma

from createMapsEtopo1 import  findSubsetIndices

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 5, 8)
__modified__ = datetime.datetime(2012, 5, 11)
__version__  = "1.0"
__status__   = "Development"

def seawinds(argv=None):
# the output array to write will be nlats x nlons
nlats = 6; nlons = 12
# open a new netCDF file for writing.
ncfile = Dataset('sfc_pres_temp.nc','w') 
# latitudes and longitudes of grid
lats_out = -25.0 + 5.0*arange(nlats,dtype='float32')
lons_out = -125.0 + 5.0*arange(nlons,dtype='float32')
# output data.
press_out = 900. + arange(nlats*nlons,dtype='float32') # 1d array
press_out.shape = (nlats,nlons) # reshape to 2d array
temp_out = 9. + 0.25*arange(nlats*nlons,dtype='float32') # 1d array
temp_out.shape = (nlats,nlons) # reshape to 2d array
# create the lat and lon dimensions.
ncfile.createDimension('latitude',nlats)
ncfile.createDimension('longitude',nlons)
# Define the coordinate variables. They will hold the coordinate
# information, that is, the latitudes and longitudes.
lats = ncfile.createVariable('latitude',dtype('float32').char,('latitude',))
lons = ncfile.createVariable('longitude',dtype('float32').char,('longitude',))
# Assign units attributes to coordinate var data. This attaches a
# text attribute to each of the coordinate variables, containing the
# units.
lats.units = 'degrees_north'
lons.units = 'degrees_east'
# write data to coordinate vars.
lats[:] = lats_out
lons[:] = lons_out
# create the pressure and temperature variables 
press = ncfile.createVariable('pressure',dtype('float32').char,('latitude','longitude'))
temp = ncfile.createVariable('temperature',dtype('float32').char,('latitude','longitude'))
# set the units attribute.
press.units =  'hPa'
temp.units = 'celsius'
# write data to variables.
press[:] = press_out
temp[:] = temp_out
# close the file.
ncfile.close()
print '*** SUCCESS writing example file sfc_press_temp.nc!'