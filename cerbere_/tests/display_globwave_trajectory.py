'''
Created on 18 juin 2012

Test reading/displaying swath file

@author: jfpiolle
'''
import logging


import datamodel.trajectory
import mapper.ncfile




logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s', datefmt='%d/%m/%Y %I:%M:%S')



ncf = mapper.ncfile.NCFile( URL = "./data/GW_L2P_ALT_ENVI_GDR_20110101_061209_20110101_070215_098_078.nc" )
traj = datamodel.trajectory.Trajectory()
traj.load( ncf )


traj.display_map('swh_calibrated', pretty=True, range=[0,8.,0.05], output='trajectory.png')