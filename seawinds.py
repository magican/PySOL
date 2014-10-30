#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 20:10:34 2012

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

def seawinds(lonStart=3, lonEnd=4, latStart=41, latEnd=42, \
            wantedDate = datetime.datetime(2010,12,18,17,39,0), \
            m=None, name='seawinds', contour=None):

    refDate = datetime.datetime(1978,1,1,0,0,0)
    #wantedDate = datetime.datetime(2010,12,18,17,39,0)

    base = '/media/SOLabNFS/store/satellite/seawinds/SI/uv/6hrly/netcdf/2000s/'
    fn = 'uv' + wantedDate.strftime('%Y%m%d') + '.nc'

    cdf = Dataset(base + fn)

    # Read lats/lons and time
    lats = cdf.variables["lat"][:]
    lons = cdf.variables["lon"][:]
    time=cdf.variables["time"][:]
    
    # Find which date fits us best
    resolution = [ datetime.timedelta(hours=0), datetime.timedelta(hours=6), \
                datetime.timedelta(hours=12), datetime.timedelta(hours=18) ]
    roundTime = datetime.timedelta(seconds=wantedDate.second,\
                        minutes=wantedDate.minute, hours=wantedDate.hour)
    closestTime = sorted(resolution, key=lambda t: abs(roundTime - t))
    t = closestTime.index(min(closestTime))
    closestTime = refDate + datetime.timedelta(hours=int(time[t]))
    print "Closest time: ", closestTime
    
    # find subset not to import all the data
    res = findSubsetIndices(latStart,latEnd,lonStart,lonEnd,lats,lons)

    lon,lat=meshgrid(lons[res[0]:res[1]],lats[res[2]:res[3]])
    u = cdf.variables["u"][t,:,int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    v = cdf.variables["v"][t,:,int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    print "Extracted data for area: (%s,%s) to (%s,%s)"%(lon.min(),lat.min(),lon.max(),lat.max())

    # Squeezing
    u = ma.squeeze(u)
    v = ma.squeeze(v)

    # Get the speed
    w = ma.sqrt(u**2+v**2)

    print "Preparing Basemap"

    # setup nsper basemap
    # Lat/Lon coords of image corners
    if lonStart< 0 and lonEnd < 0:
        lon_0= - (abs(lonEnd)+abs(lonStart))/2.0
    else:
        lon_0=(abs(lonEnd)+abs(lonStart))/2.0

    if latStart< 0 and latEnd < 0:
        lat_0= - (abs(latEnd)+abs(latStart))/2.0
    else:
        lat_0=(abs(latEnd)+abs(latStart))/2.0

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

    x, y = m(lon,lat)

    #!!!Check!!!
    #Compared Q and QQ - no difference
    #uproj,vproj,xx,yy = \
    #m.transform_vector(u[t],v[t],lon[-1,:],lat[:,-1],13,13,returnxy=True,masked=True)
    #Q = m.quiver(xx,yy,uproj,vproj)
    #QQ = m.quiver(x,y,u[t],v[t])

    plt.close('all')
    plt.figure()
    # maximizing figure
    # mng = plt.get_current_fig_manager()
    # mng.resize(1920,1080)
    if contour is 'fill':
        CS1 = m.contourf(x,y,w)
    else:
        CS1 = m.pcolor(x,y,w)
    Q = m.quiver(x,y,u,v) #Checked in matlab should be right
    # make quiver key
    plt.quiverkey(Q, 0.07, -0.03, 10, r'$10 \frac{m}{s}$', labelpos='W')

    CS1.axis='tight'
    plt.jet()
    cb = plt.colorbar(CS1,orientation='vertical',extend='both', shrink=0.8)
    cb.set_label('Seawinds [m/s]')
    m.drawmeridians(arange(lon.min(),lon.max(),0.5),labels=[0,0,0,1])
    m.drawparallels(arange(lat.min(),lat.max(),0.5),labels=[1,0,0,0])
    m.drawcoastlines()

    plt.draw()

    plt.savefig('/home/mag/' + name + closestTime.strftime('%Y%m%d') + '.png', \
             facecolor='w', edgecolor='w', dpi=100, bbox_inches="tight", pad_inches=1.75)
    #plt.show()
    #plt.close()

if __name__ == "__main__":

    seawinds(lonStart=3, lonEnd=4, \
            latStart=41, latEnd=42, \
            m=None, name='seawinds', contour=None)
