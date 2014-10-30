#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 13:36:16 2012

@author: mag
"""

import bz2
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime
from os import path

from numpy import meshgrid, arange, squeeze, floor, ceil

from createMapsEtopo1 import  findSubsetIndices

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 5, 7)
__modified__ = datetime.datetime(2012, 5, 11)
__version__  = "1.0"
__status__   = "Development"

def g1sst(lonStart=3, lonEnd=4, latStart=41, latEnd=42, \
            wantedDate = datetime.datetime(2010,12,18,17,39,0), \
            m=None, name='g1sst', contour=None):

    base = '/media/SOLabNFS/store/satellite/jpl_ourocean-l4uhfnd_glob_g1sst/G1SST/'
    outdir = '/home/mag/'

    #~ wantedDate = datetime.datetime(2010,12,18,17,39,0)
    wantedDay = wantedDate.timetuple().tm_yday

    fn = wantedDate.strftime('%Y%m%d') + \
                            '-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.nc.bz2'

    # Check if the data already decompressed
    if not path.isfile(outdir + fn[:-3]):
        print "Extracting archive to:\n" , outdir + fn[:-3]
        decompresseddata=bz2.BZ2File(base + wantedDate.strftime('%Y') + '/' \
                                     + str(wantedDay) + '/' + fn, "rb").read()
        f_out=open(outdir + fn[:-3], "wb")
        f_out.write(decompresseddata)
        f_out.close()
        del decompresseddata

    cdf = Dataset(outdir + fn[:-3])

    # Read lats/lons and time
    lats = cdf.variables["lat"][:]
    lons = cdf.variables["lon"][:]

    ## Find which date fits us best
    #for t in range(len(time)):
    #	currentTime=refDate + datetime.timedelta(hours=int(time[t]))
    #	#print currentTime
    #	if currentTime.year in wantedYears and currentTime.month in wantedMonths:
    #	print "Current time inspecting: %s"%(currentTime)
    #	myIndex=t; dateString="%s"%(currentTime)

    # find subset not to import all the data
    res = findSubsetIndices(latStart,latEnd,lonStart,lonEnd,lats,lons)

    lon,lat=meshgrid(lons[res[0]:res[1]],lats[res[2]:res[3]])
    sst = cdf.variables["analysed_sst"][:,int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    print "Extracted data for area: (%s,%s) to (%s,%s)"%(lon.min(),lat.min(),lon.max(),lat.max())

    # Squeezing and converting to Celsius
    sst = squeeze(sst) - 273.15

    if lonStart< 0 and lonEnd < 0:
        lon_0= - (abs(lonEnd)+abs(lonStart))/2.0
    else:
        lon_0=(abs(lonEnd)+abs(lonStart))/2.0

    if latStart< 0 and latEnd < 0:
        lat_0= - (abs(latEnd)+abs(latStart))/2.0
    else:
        lat_0=(abs(latEnd)+abs(latStart))/2.0

    print 'Center longitude ', lon_0
    print 'Center latitude ', lat_0

    if m is None:
        print ( "Using default NSPER Basemap \n")
#        m = Basemap(llcrnrlat=latStart,urcrnrlat=latEnd,\
#                llcrnrlon=lonStart,urcrnrlon=lonEnd,\
#                rsphere=(6378137.00,6356752.3142),\
#                resolution='h',area_thresh=1000.,projection='lcc',\
#                lat_1=latStart,lon_0=lon_0)
        m = Basemap(llcrnrlat=latStart,urcrnrlat=latEnd,\
                llcrnrlon=lonStart,urcrnrlon=lonEnd,\
                satellite_height=798000, \
                resolution='h',area_thresh=1000.,projection='nsper',\
                lat_1=latStart,lon_0=lon_0,lat_0=lat_0)

    color = 'black'

    x, y = m(lon,lat)

    # maximizing figure
    mng = plt.get_current_fig_manager()
    mng.resize(1920,1080)

    if contour is 'fill':
        CS1 = m.contourf(x, y, sst)
    else:
        CS1 = m.pcolor(x, y, sst)

    CS1.axis='tight'
    plt.jet()
    # Set the pad=0.23 so 2 colorbars won't intersect
    cb = plt.colorbar(CS1, shrink=0.95, extend='both', \
                      orientation='horizontal', pad=0.1)
    # A working example (for any value range) with 5 ticks along the bar is:
    clm=(floor(sst.min()), ceil(sst.max()))
    m0=(clm[0])            # colorbar min value
    m4=(clm[1])             # colorbar max value
    m1=(1*(m4-m0)/4.0 + m0)               # colorbar mid value 1
    m2=(2*(m4-m0)/4.0 + m0)               # colorbar mid value 2
    m3=(3*(m4-m0)/4.0 + m0)               # colorbar mid value 3
    cb.set_ticks([m0,m1,m2,m3,m4])
    cb.set_ticklabels([m0,m1,m2,m3,m4])
    cb.update_ticks()
    cb.set_label('G1SST [deg C]')
    plt.draw()

    # set the step of Lat/Lon to plot
    stepLon = round((lon.max()-lon.min())/3, 1)
    stepLat = round((lat.max()-lat.min())/3, 1)
    m.drawmeridians(arange(round(lon.min(),1),round(lon.max(),1), stepLon), \
                    labels=[0,0,0,1], color=color, dashes=[5,5])
    m.drawparallels(arange(round(lat.min(),1),round(lat.max(),1), stepLat), \
                    labels=[1,0,0,0], color=color, dashes=[5,5], rotation=90)
#    m.drawmeridians(arange(lon.min(),lon.max(),0.2),labels=[0,0,1,1],linewidth=0)
#    m.drawparallels(arange(lat.min(),lat.max(),0.2),labels=[1,1,0,0],linewidth=0)
    m.drawcoastlines()

    plt.draw()

    plt.savefig('/home/mag/' + name + wantedDate.strftime('%Y%m%d') + '.png', \
    		 facecolor='w', edgecolor='w', dpi=100, bbox_inches="tight", pad_inches=1.75)
    #plt.show()
    #plt.close()

if __name__ == "__main__":

    g1sst(lonStart=3, lonEnd=4, \
            latStart=41, latEnd=42, \
            m=None, name='g1sst', contour=None)
