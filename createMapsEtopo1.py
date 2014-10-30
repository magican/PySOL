
import datetime
import numpy as np
#from numpy import arange, round
#from numpy import arange
from netCDF4 import Dataset
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap

from gmtColormap import gmtColormap

import laplaceFilter
from mpl_util import LevelColormap

#~ import ipdb

# Originally created by Trond Kristiansen,
# modified  datetime.datetime(2012, 5, 3) by Alexander Myasoedov

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2008, 8, 15)
__modified__ = datetime.datetime(2009, 7, 21)
__version__  = "1.0"
__status__   = "Development"

def findSubsetIndices(min_lat,max_lat,min_lon,max_lon,lats,lons):
    """Array to store the results returned from the function"""
    res=np.zeros((4),dtype=np.float64)
    minLon=min_lon; maxLon=max_lon

    distances1 = []; distances2 = []
    indices=[]; index=1

    for point in lats:
        s1 = max_lat-point # (vector subtract)
        s2 = min_lat-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index-1))
        index=index+1

    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])

    distances1 = []; distances2 = []; index=1

    for point in lons:
        s1 = maxLon-point # (vector subtract)
        s2 = minLon-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index-1))
        index=index+1

    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])

    # Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices
    minJ=indices[1][2]
    maxJ=indices[0][2]
    minI=indices[3][2]
    maxI=indices[2][2]

    res[0]=min(minI,maxI); res[1]=max(minI,maxI); res[2]=min(minJ,maxJ); res[3]=max(minJ,maxJ)
    return res

def makeMap(lonStart=1, lonEnd=5, \
            latStart=37, latEnd=47, \
            m=None, name='etopo1map', contour=None, cb=None):

    # Get the etopo1 data
    etopo1name='/media/SOLabNFS/store/auxdata/dem/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gmt4.grd'
    etopo1 = Dataset(etopo1name,'r')

    lons = etopo1.variables["x"][:]
    lats = etopo1.variables["y"][:]

    res = findSubsetIndices(latStart,latEnd,lonStart,lonEnd,lats,lons)

    lon,lat=np.meshgrid(lons[res[0]:res[1]],lats[res[2]:res[3]])
    print "Extracted data for area %s : (%s,%s) to (%s,%s)"%(name,lon.min(),lat.min(),lon.max(),lat.max())
    bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    bathy = laplaceFilter.laplace_filter(bathy,M=None)

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

    x, y = m(lon,lat)

    if contour is None:
        levels = [-6000,-5000,-3000, -2000, -1500, -1000, -500, \
                  -400, -300, -250, -200, -150, -100, -75, -65, -50, -35, \
                  -25, -15, -10, -5, \
                  0, 5, 10, 15, 25, 35, 50, 65, 75, 100, \
                  150, 200, 250, 300, 400, 500, 1000, \
                  2000, 3000, 5000, 6000]

        cm = gmtColormap(fileName='GMT_globe', \
        GMTPath='/home/mag/Documents/repos/solab/PySOL/cmap_data/')
        cs = m.contourf(x,y,bathy,levels,
                           cmap=LevelColormap(levels,cmap=cm),
                           alpha=1.0,
                           extend='both')
        cs.axis='tight'
        # new axis for colorbar.
        ax = plt.gca()
        l,b,w,h=ax.get_position().bounds
        if cb is not None:
            cax = plt.axes([l-0.25, b+h*0.1, 0.025, h*0.8]) # setup colorbar axes
            cb = plt.colorbar(cs, cax, orientation='vertical')
            cb.set_label('Height [m]')
#           cax = plt.axes([l+w*0.1, b-0.05, w*0.8, 0.025]) # setup colorbar axes
#           cb = plt.colorbar(cs, cax, orientation='horizontal')
            plt.axes(ax)  # make the original axes current again
    elif contour is 'land':
        levels = [0, 5, 10, 15, 25, 35, 50, 65, 75, 100, \
                  150, 200, 250, 300, 400, 500, 1000, \
                  2000, 3000, 5000, 6000]

        cm = gmtColormap(fileName='PySOL_land', \
                GMTPath='/home/mag/Documents/repos/solab/PySOL/cmap_data/', \
                frm='mid')
        cs = m.contourf(x,y,bathy,levels,
                           cmap=LevelColormap(levels,cmap=cm),
                           alpha=1.0,
                           extend='max')
        cs.axis='tight'
        # new axis for colorbar.
        ax = plt.gca()
        if cb is not None:
            l,b,w,h=ax.get_position().bounds
            cax = plt.axes([l-0.25, b+h*0.1, 0.025, h*0.8]) # setup colorbar axes
            cb = plt.colorbar(cs, cax, orientation='vertical')
            cb.set_label('Height [m]')
#            cax = plt.axes([l+w*0.1, b-0.05, w*0.8, 0.025]) # setup colorbar axes
#            cb = plt.colorbar(cs, cax, orientation='horizontal')
            plt.axes(ax)  # make the original axes current again
#        cb.ax.yaxis.set_ylabel_position('left')
    elif contour is 'ocean':
        levels = [-6000,-5000,-3000, -2000, -1500, -1000, -500, \
                  -400, -300, -250, -200, -150, -100, -75, -65, -50, -35, \
                  -25, -15, -10, -5, 0]

        levels = [-2000, -1600, -1200, -800, -100]

        cm = gmtColormap(fileName='GMT_ocean', \
                GMTPath='/home/mag/Documents/repos/solab/PySOL/cmap_data/', \
                to='mid')
        cs2 = m.contour(x,y,bathy,levels,
                           cmap=LevelColormap(levels,cmap=cm),
                           alpha=1.0,
                           extend='min')
        cs2.axis='tight'
        plt.clabel(cs2, levels, fmt = '%i', colors = 'k', \
                   manual=0, inline=0)
#    m.fillcontinents(color='coral',lake_color='aqua')
#    m.drawmeridians(arange(round(lons.min(),1),round(lons.max(),1), 0.5), \
#                    labels=[0,0,0,1], color='k')
#    m.drawparallels(arange(round(lats.min(),1),round(lats.max(),1), 0.5), \
#                    labels=[1,0,0,0], color='k')
    m.drawcoastlines(linewidth=1.25,color='k')
    m.drawcountries(linewidth=1.25,color='k')

    # maximizing figure
    mng = plt.get_current_fig_manager()
    mng.resize(1920,1080)
    
    plt.draw()
    #~ plt.show()
    #~ ipdb.set_trace()

    name = name.split(' ')[0] # split the name if it has spaces
    plotfile='/home/mag/'+str(name)+'_ETOPO1.tiff'
    plt.savefig(plotfile, facecolor='w', edgecolor='w', \
                dpi=300, bbox_inches="tight", pad_inches=0.1)
#    plt.show()

if __name__ == "__main__":

    makeMap(lonStart=1, lonEnd=5, \
            latStart=37, latEnd=47, \
            m=None, name='etopo1map', contour=None, cb=None)
