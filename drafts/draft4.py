#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 00:06:38 2012

@author: mag
"""

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
# setup lambert conformal basemap.
# lat_1 is first standard parallel.
# lat_2 is second standard parallel (defaults to lat_1).
# lon_0,lat_0 is central point.
# rsphere=(6378137.00,6356752.3142) specifies WGS4 ellipsoid
# area_thresh=1000 means don't plot coastline features less
# than 1000 km^2 in area.
m = Basemap(width=12000000,height=9000000,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',projection='cyl')
figure()
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='aqua') 




from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
# setup Lambert Conformal basemap.
m = Basemap(width=12000000,height=9000000,projection='lcc',
            resolution='l',lat_1=45.,lat_2=55,lat_0=50,lon_0=-107.)

m = Basemap(llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
            llcrnrlon=ll_lon, urcrnrlon=ur_lon, \
            rsphere=(6378137.00,6356752.3142),\
            resolution='i',projection='lcc',\
            lat_0=lat.mean(),lon_0=lon.mean())
# draw a boundary around the map, fill the background.
# this background will end up being the ocean color, since
# the continents will be drawn on top.
m.drawmapboundary(fill_color='aqua') 
# fill continents, set lake color same as ocean color. 
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
# label parallels on right and top
# meridians on bottom and left
parallels = np.arange(40.,47,0.2)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(0.,5.,0.2)
m.drawmeridians(meridians,labels=[True,False,False,True])
# plot blue dot and label it as such.
lon1, lat1 = 3.4, 42.0 # Location of Boulder
# convert to map projection coords. 
# Note that lon,lat can be scalars, lists or numpy arrays.
xpt,ypt = m(lon,lat) 
xpt1,ypt1 = m(lon1,lat1) 
# convert back to lat/lon
lonpt, latpt = m(xpt,ypt,inverse=True)
m.contourf(xpt,ypt,S_VV_ABS)
m.plot(xpt1,ypt1,'bo')  # plot a blue dot there
# put some text next to the dot, offset a little bit
# (the offset is in map projection coordinates)
plt.text(xpt1+1000,ypt1+1000,'Test (%5.1fW,%3.1fN)' % (lon1,lat1))
plt.savefig('plotboulder.png')



lats = lat
lons = lon
# shift lats, lons so values represent edges of grid boxes
# (as pcolor expects).
delon = lons[1]-lons[0]; delat = lats[1]-lats[0]
lons = (lons - 0.5*delon).tolist()
lons.append(lons[-1]+delon)
lons = np.array(lons,np.float64)
lats = (lats - 0.5*delat).tolist()
lats.append(lats[-1]+delat)
lats = np.array(lats,np.float64)
# creat figure, axes instances.
fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
# create Basemap instance for Robinson projection.
# coastlines not used, so resolution set to None to skip
# continent processing (this speeds things up a bit)
m = Basemap(projection='robin',lon_0=lons.mean(),resolution='l')
# compute map projection coordinates of grid.
x, y = m(*np.meshgrid(lons, lats))
# draw line around map projection limb.
# color background of map projection region.
# missing values over land will show up this color.
m.drawmapboundary(fill_color='0.3')
# plot sst, then ice with pcolor
im2 = m.pcolor(x,y,S_VV_ABS,shading='flat',cmap=plt.cm.gist_gray)
# draw parallels and meridians, but don't bother labelling them.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
# add colorbar
cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
# add a title.
ax.set_title('SST and ICE analysis for %s'%date)
plt.savefig('plotsst.png')
