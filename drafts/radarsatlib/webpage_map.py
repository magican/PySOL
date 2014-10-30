#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 00:05:56 2011

@author: sat kumar tomer (http://civil.iisc.ernet.in/~satkumar/)

this program make the map for the www.ambhas.com
for the HH, 0.5*(HV+VH), and VV data
"""

from mpl_toolkits.basemap import Basemap
import numpy as np
from gis import utm2deg
from pylab import *
import datetime as dt
import gdal
from gdalconst import *

fidSigma = gdal.Open("/home/tomer/RADARSAT/04Mar10/100m//filtered",GA_ReadOnly)
# read the projection details, because the projection details of the backscattering image
#filtered and IncidenceAngle are same, they are read only from one file
GT = fidSigma.GetGeoTransform()
    
# read the HH, HV, VH, VV 
HH = fidSigma.GetRasterBand(1).ReadAsArray()
HV = fidSigma.GetRasterBand(2).ReadAsArray()
VH = fidSigma.GetRasterBand(3).ReadAsArray()
VV = fidSigma.GetRasterBand(4).ReadAsArray()

# make the co ordinate for the kabini
x_ber = np.linspace(664050, 684950, 210)
y_ber = np.linspace(1308950, 1294050, 150)
XI, YI = np.meshgrid(x_ber, y_ber)
Lon,Lat = utm2deg(XI,YI)

#************************ HH ***********************************************
# make the base map
m = Basemap(projection='merc',llcrnrlat=11.72,urcrnrlat=11.825,\
            llcrnrlon=76.51,urcrnrlon=76.67,lat_ts=20,resolution='i')
m.drawcoastlines()
# draw parallels and meridians.
m.drawparallels(np.arange(11.7,11.9,.05),labels=[1,0,0,0])
m.drawmeridians(np.arange(76.4,76.8,.05),labels=[0,0,0,1])

# read the shapefile archive
s = m.readshapefile('/home/tomer/Berambadi/shapefiles/watershedfnl1','berambadi',
    color='k',linewidth=1.5)

# compute native map projection coordinates of lat/lon grid.
x, y = m(Lon,Lat)

# contour data over the map
cs = m.pcolor(x,y,HH,cmap=plt.cm.jet)
cb = plt.colorbar(cs, shrink=0.8, extend='both')

title(" 04/03/2010  (resolution = 100 M)")
savefig('/home/tomer/webpage/AMBHAS/html/images/HH.png')
close()

#************************ 0.5(HV+VH)********************************************
# make the base map
m = Basemap(projection='merc',llcrnrlat=11.72,urcrnrlat=11.825,\
            llcrnrlon=76.51,urcrnrlon=76.67,lat_ts=20,resolution='i')
m.drawcoastlines()
# draw parallels and meridians.
m.drawparallels(np.arange(11.7,11.9,.05),labels=[1,0,0,0])
m.drawmeridians(np.arange(76.4,76.8,.05),labels=[0,0,0,1])

# read the shapefile archive
s = m.readshapefile('/home/tomer/Berambadi/shapefiles/watershedfnl1','berambadi',
    color='k',linewidth=1.5)

# compute native map projection coordinates of lat/lon grid.
x, y = m(Lon,Lat)

# contour data over the map
cs = m.pcolor(x,y,0.5*(HV+VH),cmap=plt.cm.jet)
cb = plt.colorbar(cs, shrink=0.8, extend='both')

title(" 04/03/2010  (resolution = 100 M)")
savefig('/home/tomer/webpage/AMBHAS/html/images/HV.png')
close()

#************************ VV********************************************
# make the base map
m = Basemap(projection='merc',llcrnrlat=11.72,urcrnrlat=11.825,\
            llcrnrlon=76.51,urcrnrlon=76.67,lat_ts=20,resolution='i')
m.drawcoastlines()
# draw parallels and meridians.
m.drawparallels(np.arange(11.7,11.9,.05),labels=[1,0,0,0])
m.drawmeridians(np.arange(76.4,76.8,.05),labels=[0,0,0,1])

# read the shapefile archive
s = m.readshapefile('/home/tomer/Berambadi/shapefiles/watershedfnl1','berambadi',
    color='k',linewidth=1.5)

# compute native map projection coordinates of lat/lon grid.
x, y = m(Lon,Lat)

# contour data over the map
cs = m.pcolor(x,y,VV,cmap=plt.cm.jet)
cb = plt.colorbar(cs, shrink=0.8, extend='both')

title(" 04/03/2010  (resolution = 100 M)")
savefig('/home/tomer/webpage/AMBHAS/html/images/VV.png')
close()