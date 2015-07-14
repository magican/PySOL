'''
Created on 30 juil. 2012

@author: jfpiolle
'''
import logging
import os

import numpy
import matplotlib.pyplot as plt

import datamodel.grid
import mapper.ncfile

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s', datefmt='%d/%m/%Y %I:%M:%S')

logging.debug('Create coordinates')
glons = numpy.arange(-180,180,0.1)
glats = numpy.arange(-90,90,0.1)

logging.debug('instantiate target grid')        
targetGrid = datamodel.grid.Grid(longitudes=(-180,180,0.1),latitudes=(-90,90,0.1))

logging.debug('read input grid')

ncf = mapper.ncfile.NCFile( URL = "./data/20120109000000-OSISAF-L3C_GHRSST-SSTsubskin-SEVIRI_SST-ssteqc_meteosat09_20120109_000000-v02.0-fv01.0.nc")
inputGrid = datamodel.grid.Grid()
inputGrid.load( ncf )

targetGrid.resample( inputGrid, variables=['sea_surface_temperature'],useSrcTime=True )


if os.path.exists('./test.nc'):
    os.remove( './test.nc' )
targetncf = mapper.ncfile.NCFile( URL = "./test.nc", mode=mapper.abstractmapper.WRITE_NEW )
targetGrid.save( targetncf )
