'''
Created on 18 juin 2012

Test reading/displaying time series file

@author: jfpiolle
'''
import logging

import matplotlib.pyplot as plt
import matplotlib.colors

from cerbere.datamodel.pointtimeseries import PointTimeSeries
from cerbere.mapper.ncfile import NCFile




logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s', datefmt='%d/%m/%Y %I:%M:%S')


ff = '/home/cercache/project/globwave/data/globwave/oceansites/ocean_temperature_sensor/2015/05/WMO61280_20150501T0000_20150515T0100_Lat_40.69N_Lon_1.48E.nc'
ncf = NCFile( url = ff)
traj = PointTimeSeries()
traj.load( ncf )


#val = traj.getData('greyscale_threshold')
traj.display_timeseries('sea_surface_temperature', palette=plt.cm.gray, pretty=True, range=[273,305,0.25], output='toto.png')
