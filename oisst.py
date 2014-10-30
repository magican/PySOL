#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 20:10:34 2012

@author: mag
"""

import gzip
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime

from numpy import meshgrid, arange, ma

from createMapsEtopo1 import  findSubsetIndices

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 5, 7)
__modified__ = datetime.datetime(2012, 5, 11)
__version__  = "1.0"
__status__   = "Development"

base = '/media/SOLabNFS/store/satellite/oisst/NetCDF/2010/AVHRR-AMSR/'
fn = 'amsr-avhrr-v2.20101218.nc.gz'
outdir = '/home/mag/'

refDate = datetime.datetime(1978,1,1,0,0,0)
#wantedDate = datetime.datetime(2010,12,18,17,39,0)


decompresseddata=gzip.open(base + fn, "rb").read()
f_out=open(outdir + fn[:-3], "wb")
f_out.write(decompresseddata)
f_out.close()
del decompresseddata

cdf = Dataset(outdir + fn[:-3])

# Read lats/lons and time
lats = cdf.variables["lat"][:]
lons = cdf.variables["lon"][:]
time=cdf.variables["time"][:]

## Find which date fits us best
#for t in range(len(time)):
#	currentTime=refDate + datetime.timedelta(hours=int(time[t]))
#	#print currentTime
#	if currentTime.year in wantedYears and currentTime.month in wantedMonths:
#	print "Current time inspecting: %s"%(currentTime)
#	myIndex=t; dateString="%s"%(currentTime)
            
# find subset not to import all the data
res = findSubsetIndices(41,43,2,4,lats,lons)

lon,lat=meshgrid(lons[res[0]:res[1]],lats[res[2]:res[3]])
sst = cdf.variables["sst"][:,:,int(res[2]):int(res[3]),int(res[0]):int(res[1])]
ice = cdf.variables["ice"][:,:,int(res[2]):int(res[3]),int(res[0]):int(res[1])]
print "Extracted data for area: (%s,%s) to (%s,%s)"%(lon.min(),lat.min(),lon.max(),lat.max())

# Squeezing
sst = ma.squeeze(sst)
ice = ma.squeeze(ice)

print "Preparing Basemap"

# setup nsper basemap
# Lat/Lon coords of image corners
ll_lat = lat.min()
ur_lat = lat.max()
ll_lon = lon.min()
ur_lon = lon.max()
cent_lat = lat.mean()
cent_lon = lon.mean()
m = Basemap(llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
				llcrnrlon=ll_lon, urcrnrlon=ur_lon, \
				resolution='h', projection='nsper', \
				satellite_height=798000, \
				lat_0=cent_lat,lon_0=cent_lon)

x, y = m(lon,lat)

plt.figure()
# maximizing figure
mng = plt.get_current_fig_manager()
mng.resize(1920,1080)
#CS1 = m.contourf(x,y,sst)
CS1 = m.pcolor(x,y,sst)

CS1.axis='tight'
plt.colorbar(CS1,orientation='vertical',extend='both', shrink=0.8)
m.drawmeridians(arange(lon.min(),lon.max(),1),labels=[0,0,0,1])
m.drawparallels(arange(lat.min(),lat.max(),1),labels=[1,0,0,0])
m.drawcoastlines()

plt.draw()

plt.savefig('/home/mag/oisst.png', \
		 facecolor='w', edgecolor='w', pad_inches=0.05, dpi=100, \
            bbox_inches="tight", pad_inches=1.75)	
#plt.show()
#plt.close()
