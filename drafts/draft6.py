#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 15:29:19 2012

@author: mag
"""

from mpl_toolkits.basemap import Basemap
import numpy as np

from scipy.interpolate import interp2d, RectBivariateSpline, griddata
from scipy import ogrid, sin, mgrid, ndimage, array

from pylab import *
import datetime as dt
import gdal
from xml.dom import minidom
from gdalconst import GA_ReadOnly

from gis import utm2deg, deg2utm
from ll2UTM import LLtoUTM

pn = '/media/data/data/OTHER/RS2 Agulhas and Lion/RS2_FQA_1xQGSS20101218_173930_00000005/'

# default values
xOff=0
yOff=0
xS=None
yS=None
xBufScale=None
yBufScale=None
s = 'abs'

dataset = gdal.Open("RADARSAT_2_CALIB:SIGMA0:" + pn + "product.xml")
gcpproj = dataset.GetGCPProjection()
RasterXSize = dataset.RasterXSize
RasterYSize = dataset.RasterYSize

dataset = gdal.Open(pn + 'temp1.tif', GA_ReadOnly)

# read the projection details, because the projection details of the backscattering image
#filtered and IncidenceAngle are same, they are read only from one file
GT = dataset.GetGeoTransform()
    
# read the HH, HV, VH, VV 
S_VV = dataset.GetRasterBand(2).ReadAsArray(xoff=xOff, yoff=yOff, \
    win_xsize=xS, win_ysize=yS, \
    buf_xsize=xBufScale, buf_ysize=yBufScale)
S_VV_ABS = absolute(S_VV)**2

xmldoc = minidom.parse(pn + "product.xml")
latGrid = xmldoc.getElementsByTagName('latitude')
lonGrid = xmldoc.getElementsByTagName('longitude')
lineGrid = xmldoc.getElementsByTagName('line')
pixelGrid = xmldoc.getElementsByTagName('pixel')
lat = zeros( (latGrid.length, 1) )
lon = zeros( lat.shape )
line = zeros( lat.shape )
pixel = zeros( lat.shape )
for n in range(latGrid.length):
    lat[n] = latGrid[n].childNodes[0].data
    lon[n] = lonGrid[n].childNodes[0].data
    line[n] = lineGrid[n].childNodes[0].data
    pixel[n] = pixelGrid[n].childNodes[0].data
       
ind = find(pixel == 0)
pixel = reshape(pixel, (ind.size, latGrid.length/ind.size))
line = reshape(line, (ind.size, latGrid.length/ind.size))
lat = reshape(lat, (ind.size, latGrid.length/ind.size))
lon = reshape(lon, (ind.size, latGrid.length/ind.size))

ll_lat = lat.min()
ur_lat = lat.max()
ll_lon = lon.min()
ur_lon = lon.max()


currtime = time()
gx = np.linspace(0, RasterXSize-1, round(RasterXSize/1))
gy = np.linspace(0, RasterYSize-1, round(RasterYSize/1))

flat = interp2d(pixel[0,:], line[:,0], lat, kind='linear')
flon = interp2d(pixel[0,:], line[:,0], lon, kind='linear')

lat_new = flat(gx, gy)
lon_new = flon(gx, gy)
print 'Parallel: time elapsed: %f' % ( time() - currtime )

## make the co ordinate for the kabini
#x_ber = np.linspace(664050, 684950, 210)
#y_ber = np.linspace(1308950, 1294050, 150)
#XI, YI = np.meshgrid(x_ber, y_ber)
#Lon,Lat = utm2deg(XI,YI)

(UTMZone, UTMEasting, UTMNorthing) = LLtoUTM(23,lat[-1,-1],lon[-1,-1])

X,Y = deg2utm(lon, lat, UTMZone[0:2])
x_lion = np.linspace(X.min(), X.max(), round(RasterXSize/8))
y_lion = np.linspace(Y.min(), Y.max(), round(RasterYSize/8))
XI, YI = np.meshgrid(x_lion, y_lion)
Lon,Lat = utm2deg(XI, YI, UTMZone[0:2])

## make the co ordinate for the region
#lon = np.linspace(GT[0]+GT[1]/2, GT[0]+GT[1]*(S_VV.shape[1]-0.5), S_VV.shape[1])
#lat = np.linspace(GT[3]+GT[5]/2, GT[3]+GT[5]*(S_VV.shape[0]-0.5), S_VV.shape[0])
#Lon, Lat = np.meshgrid(lon, lat)

#************************ HH ***********************************************
# make the base map
m = Basemap(llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
            llcrnrlon=ll_lon, urcrnrlon=ur_lon, \
            area_thresh=1000, \
            rsphere=(6378137.00,6356752.3142451793),\
            resolution='h',projection='lcc',\
            lat_1=lat.min(), lat_0=lat.mean(),lon_0=lon.mean())
#
#m = Basemap(projection='cyl',llcrnrlat=11.72,urcrnrlat=11.825,\
#            llcrnrlon=76.51,urcrnrlon=76.67,lat_ts=20,resolution='i')
#
#m = Basemap(projection='merc', llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
#            llcrnrlon=ll_lon, urcrnrlon=ur_lon, resolution='h')

m.drawcoastlines(linewidth=1.25,color='k')
m.drawcountries(linewidth=1.25,color='k')
# draw parallels and meridians.
m.drawparallels(np.arange(ll_lat,ur_lat,.05),labels=[1,0,0,0])
m.drawmeridians(np.arange(ll_lon,ur_lon,.05),labels=[0,0,0,1])

## read the shapefile archive
#s = m.readshapefile('/home/tomer/Berambadi/shapefiles/watershedfnl1','berambadi',
#    color='k',linewidth=1.5)

# compute native map projection coordinates of lat/lon grid.
x, y = m(Lon,Lat)

x, y = m(lon_new,lat_new)

# contour data over the map
cs = m.imshow(S_VV_ABS[::8,::8], cmap=plt.cm.gray, origin='upper')

cs = m.pcolormesh(x,y,S_VV_ABS, cmap=plt.cm.gray)

plt.clim(0, 0.1)
cb = plt.colorbar(cs, shrink=0.8, extend='both')

title(" 04/03/2010  (resolution = 100 M)")

savefig('S_VV_ABS.png')
close()