'''
Created on 1 aout 2012

@author: jfpiolle
'''

import logging

import numpy
import matplotlib.pyplot as plt

import datamodel.swath
import datamodel.grid
import mapper.ncfile




logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s', datefmt='%d/%m/%Y %I:%M:%S')



hf = mapper.ncfile.NCFile( URL = "data/20120323-ATS_NR_2P-UPA-L2P-ATS_NR__2PNPDK20120323_185109_000062203113_00042_52642_2614-v01.nc" )


swath = datamodel.swath.Swath()
swath.load( hf )

logging.debug('instantiate target grid')

targetGrid = datamodel.grid.Grid(longitudes=(-180,180,1.),latitudes=(-90,90,1.),projectionId='regular')
targetGrid = targetGrid.resample(swath, variables=['sea_surface_temperature'],neighbours=1)

swath.display_map('sea_surface_temperature', palette=plt.cm.YlOrRd, pretty=True, range=[273.15,305.15,0.25], output='toto.png')
