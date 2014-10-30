#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 14:54:13 2012

@author: mag
"""

from numpy import linspace, arange, round, floor, ceil
from scipy.interpolate import interp2d

import matplotlib.pyplot as plt
#from matplotlib import rc as mplrc
from mpl_toolkits.basemap import Basemap

from imcrop import imcrop

import datetime

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 5, 1)
__modified__ = datetime.datetime(2012, 5, 11)
__version__  = "1.0"
__status__   = "Development"

def interpLL(l, pixel, line, gx, gy):
    """
    Interpolating lat/lon to image size for future pcolormesh
    """
    fl = interp2d(pixel[0,:], line[:,0], l, kind='linear')
    l_new = fl(gx, gy)
    return l_new

def plotRS2(data, lat, lon, pixel, line, scale=8, m=None, clm=None, \
            label='SigmaVVwnr [linear units]', ptsCrop=None):

    if clm is None:
        clm=(floor(data.min()), ceil(data.max()))

    if ptsCrop is not None:
        # Scaling the data
        data = imcrop(ptsCrop, data)
        data = data[::scale,::scale]
        # Interpolating lat/lon to image size for future pcolormesh
        RasterXSize = round((-ptsCrop[0,0]+ptsCrop[-1,0])/scale)
        RasterYSize = round((-ptsCrop[0,-1]+ptsCrop[-1,-1])/scale)
        gx = linspace(ptsCrop[0,0], ptsCrop[-1,0], RasterXSize+1)
        gy = linspace(ptsCrop[0,-1], ptsCrop[-1,-1], RasterYSize+1)
    else:
        # Scaling the data
        data = data[::scale,::scale]
        # Interpolating lat/lon to image size for future pcolormesh
        RasterXSize = round(pixel[0,-1]/scale)
        RasterYSize = round(line[-1,0]/scale)
        gx = linspace(0, pixel[0,-1], RasterXSize+1)
        gy = linspace(0, line[-1,0], RasterYSize+1)

    lat_new = interpLL(lat, pixel, line, gx, gy)
    lon_new = interpLL(lon, pixel, line, gx, gy)

    if m is None:
        print ( "Using default NSPER Basemap \n")

        # Lat/Lon coords of image corners
        ll_lat = lat_new.min()
        ur_lat = lat_new.max()
        ll_lon = lon_new.min()
        ur_lon = lon_new.max()
        cent_lat = lat_new.mean()
        cent_lon = lon_new.mean()

        # use major and minor sphere radii from WGS84 ellipsoid.
        # setup nsper basemap.
        # lat_1 is first standard parallel.
        # lat_2 is second standard parallel (defaults to lat_1).
        # lon_0,lat_0 is central point.
        # rsphere=(6378137.00,6356752.3142) specifies WGS4 ellipsoid
        # area_thresh=1000 means don't plot coastline features less
        # than 1000 km^2 in area.

        m = Basemap(llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
                    llcrnrlon=ll_lon, urcrnrlon=ur_lon, \
                    resolution='i', projection='nsper', \
                    satellite_height=798000, \
                    lat_0=cent_lat,lon_0=cent_lon)

    color = 'black'

    # compute native map projection coordinates of lat/lon grid.
    x, y = m(lon_new,lat_new)

#    # Set the figures fonts (see ~/.matplotlib/matplotlibrc)
#    fonts = {'family' : 'sans', \
#             'style'  : 'normal', \
#             'weight' : 'bold', \
#             'size'   : 17, \
#            }
#    mplrc('font', **fonts)

    # fig = plt.figure(figsize=(RasterXSize/50,RasterYSize/50), dpi=300)
#    fig.add_axes((0.1,0.1,0.8,0.8))
    # maximizing figure
    mng = plt.get_current_fig_manager()
    mng.resize(1920,1080)

    cs = m.pcolormesh(x,y,data)
    cs.axis='tight'
    plt.gray()
    plt.clim(clm)
#    cb = plt.colorbar(cs, shrink=0.8, extend='both', format='%.0e')
#    cb = plt.colorbar(cs, shrink=0.8, extend='both', format='%.2f', \
#                      orientation='vertical')
    cb = plt.colorbar(cs, shrink=0.5, extend='both', \
                      orientation='horizontal', pad=0.1, aspect=33)
    # A working example (for any value range) with 5 ticks along the bar is:
    m0=(clm[0])                      # colorbar min value
    m5=(clm[1])                      # colorbar max value
    m1=round((1*(m5-m0)/5.0 + m0),2) # colorbar mid value 1
    m2=round((2*(m5-m0)/5.0 + m0),2) # colorbar mid value 2
    m3=round((3*(m5-m0)/5.0 + m0),2) # colorbar mid value 3
    m4=round((4*(m5-m0)/5.0 + m0),2) # colorbar mid value 4
    cb.set_ticks([m0,m1,m2,m3,m4,m5])
    cb.set_ticklabels([m0,m1,m2,m3,m4,m5])
    cb.update_ticks()
    cb.set_label(label)
    plt.draw()

    # set the step of Lat/Lon to plot
    stepLon = round((lon_new.max()-lon_new.min())/3, 1)
    stepLat = round((lat_new.max()-lat_new.min())/1.5, 1)
    # fool proofing, so that ronded value is not 0, but 0.05 in case of small area
    if stepLat == 0: stepLat=0.05
    if stepLon == 0: stepLon=0.05
    m.drawmeridians(arange(round(lon_new.min(),1),round(lon_new.max(),1), stepLon), \
                    labels=[0,0,0,1], color=color, dashes=[5,5], linewidth=0)
    m.drawparallels(arange(round(lat_new.min(),1),round(lat_new.max(),1), stepLat), \
                    labels=[1,0,0,0], color=color, dashes=[5,5], linewidth=0, rotation=90)

#   # draw only labels, not the lines
#    m.drawmeridians(arange(round(lon_new.min(),1),round(lon_new.max(),1), 0.2), \
#                    labels=[0,0,0,1], color=color, linewidth=0)
#    m.drawparallels(arange(round(lon_new.min(),1),round(lon_new.max(),1), 0.2), \
#                    labels=[1,0,0,0], color=color, linewidth=0, rotation=90)
#    # set the 3 steps of Lat/Lon to plot
#    rngLon = zeros((3,))
#    rngLon[0], rngLon[2] = round(lon_new.min(),1), round(lon_new.max(),1)
#    rngLon[1] = round(((rngLon[2]-rngLon[0])/2+rngLon[0]),1)
#    rngLat = zeros((3,))
#    rngLat[0], rngLat[2] = round(lat_new.min(),1), round(lat_new.max(),1)
#    rngLat[1] = round((rngLat[2]-rngLat[0])/2+rngLat[0],1)
#
#    m.drawmeridians(rngLon, \
#                    labels=[0,0,0,1], color=color, dashes=[5,5])
#    m.drawparallels(rngLat, \
#                    labels=[1,0,0,0], color=color, dashes=[5,5], rotation=90)

    a = label.split(' ')[0] # split the name if it has spaces

    plt.tight_layout()
    plt.draw()

    plt.savefig('/home/mag/' + a + '.tiff', facecolor='w', edgecolor='w', \
                dpi=300, bbox_inches="tight", pad_inches=0.1)

if __name__ == "__main__":
    plotRS2(data, lat, lon, pixel, line, scale=8, m=None, clm=None, \
            label='SigmaVVwnr [linear units]', ptsCrop=None)
