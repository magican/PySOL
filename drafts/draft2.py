# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 17:00:13 2012

@author: mag
"""

from numpy import absolute, reshape, array, meshgrid, float64, arange
from scipy import linspace
from osgeo import gdal
from xml.dom import minidom
from matplotlib.mlab import find
from scipy.misc import imresize

import matplotlib.mpl


from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

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
        geotransform = dataset.GetGeoTransform()
        gcps = dataset.GetGCPs()
        gcpproj = dataset.GetGCPProjection()
        RasterXSize = dataset.RasterXSize
        RasterYSize = dataset.RasterYSize
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

#minLon = lon.min()
#maxLon = lon.max()
#minLat = lat.min()
#maxLat = lat.max()



# use major and minor sphere radii from WGS84 ellipsoid.
# setup lambert conformal basemap.
# lat_1 is first standard parallel.
# lat_2 is second standard parallel (defaults to lat_1).
# lon_0,lat_0 is central point.
# rsphere=(6378137.00,6356752.3142) specifies WGS4 ellipsoid
# area_thresh=1000 means don't plot coastline features less
# than 1000 km^2 in area.
overlay_color = 'black'
#m = Basemap(llcrnrlat=30, urcrnrlat=50,\
#            llcrnrlon=0, urcrnrlon=10, \
#            area_thresh=1000, \
#            rsphere=(6378137.00,6356752.3142),\
#            resolution='i',projection='lcc',\
#            lat_1=lat.min(), lat_0=lat.mean(),lon_0=lon.mean())

m = Basemap(llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
            llcrnrlon=ll_lon, urcrnrlon=ur_lon, \
            resolution='h', projection='nsper', \
            satellite_height=798000, \
            lat_0=lat.mean(),lon_0=lon.mean())

#ll_lat = lat[-1,0]
#ur_lat = lat[0,-1]
#ll_lon = lon[-1,0]
#ur_lon = lon[0,-1]
cent_lat = lat.mean()
cent_lon = lon.mean()
#
#ll_lat = +33.3434
#ur_lat = +45.0346
#ll_lon = +1.7154 
#ur_lon = +16.3414 
#
#m = Basemap(projection='cyl', llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
#            llcrnrlon=ll_lon, urcrnrlon=ur_lon, resolution='f')

fig = plt.figure(figsize=(7,7))
ax = fig.add_axes((0.1,0.1,0.8,0.8))

#cs = m.imshow(S_VV_ABS[::1,::1], interpolation='nearest', origin='upper')
cs = m.pcolormesh(x,y,S_VV_ABS)
plt.gray()
plt.clim(0, 0.1)
cb = plt.colorbar(cs, shrink=0.8, extend='both')
cb.set_label('SigmaVVwnr [linear units]')

m.drawcoastlines(linewidth=1.25,color='k')
m.drawcountries(linewidth=1.25,color='k')
m.fillcontinents(color='coral',lake_color='aqua')
m.etopo()
m.shadedrelief()

m.drawmeridians(arange(1,4,0.5),labels=[1,0,0,1],fontsize=12, color=overlay_color)
m.drawparallels(arange(40,43,0.5),labels=[0,1,1,0],fontsize=12, color=overlay_color)














# read in topo data (on a regular lat/lon grid)
etopo = np.loadtxt('etopo20data.gz')
lons  = np.loadtxt('etopo20lons.gz')
lats  = np.loadtxt('etopo20lats.gz')
# create Basemap instance for Equidistant Cylindrical Projection.
#m = Basemap(projection='cyl',lon_0=0.5*(lons[0]+lons[-1]))
# compute map projection coordinates for lat/lon grid.
xe, ye = m(*np.meshgrid(lons,lats))
# make filled contour plot.
cs = m.contourf(xe,ye,etopo,30,cmap=plt.cm.jet)
m.drawcoastlines() # draw coastlines
m.drawmapboundary() # draw a line around the map region
m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
m.drawmeridians(np.arange(0.,420.,60.),labels=[0,0,0,1]) # draw meridians
plt.title('Equidistant Cylindrical Projection') # add a title

plt.colorbar()

plt.show() 










#AreaDefinition for pyresample

ll_lat = lat[-1,0]
ur_lat = lat[0,-1]
ll_lon = lon[-1,0]
ur_lon = lon[0,-1]
cent_lat = lat.mean()
cent_lon = lon.mean()

m = Basemap(llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
            llcrnrlon=ll_lon, urcrnrlon=ur_lon, \
            resolution='i', projection='lcc', \
            rsphere=(6378137.00,6356752.3142451793), \
            lat_0=cent_lat,lon_0=cent_lon)

x, y = m(lon, lat)

x_ll = x[-1,0]
x_ur = x[0,-1]
y_ll = y[-1,0]
y_ur = y[0,-1]

from pyresample import image, geometry, kd_tree, plot

area_def = geometry.AreaDefinition('areaD', 'Europe (3km, HRV, VTC)', 'areaD', \
   {'a': '6378144.0', \
    'b': '6356759.0', \
    'lat_0': cent_lat, \
    'lon_0': cent_lon, \
    'proj': 'stere'}, \
    round(5924/8), \
    round(7930/8), \
    [-1370912.72, -909968.64000000001, 1029087.28, 1490031.3600000001])

area_def = geometry.AreaDefinition('areaD', 'Europe (3km, HRV, VTC)', 'areaD', \
   {'a': '6378144.0', \
    'b': '6356759.0', \
    'lat_0': cent_lat, \
    'lon_0': cent_lon, \
    'proj': 'stere'}, \
    800, \
    800, \
    [-1370912.72, -909968.64000000001, 1029087.28, 1490031.3600000001])

area_def = geometry.AreaDefinition('areaD', 'Europe (3km, HRV, VTC)', 'areaD', \
   {'a': '6378144.0', \
    'b': '6356759.0', \
    'lat_0': cent_lat, \
    'lon_0': cent_lon, \
    'proj': 'stere'}, \
    800, \
    800, \
    [x_ll, x_ur, y_ll, y_ur])


swath_def = geometry.SwathDefinition(lon, lat)

grid_def = geometry.GridDefinition(lons=lon, lats=lat)

overlap_fraction = swath_def.overlaps(area_def)

overlap_fraction = swath_def.overlap_rate(area_def)

result = kd_tree.resample_nearest(swath_def, S_VV_ABS, area_def, \
    radius_of_influence=200, fill_value=None)

bmap = plot.area_def2basemap(area_def)









            
m = Basemap(projection='ortho',lat_0=lat.mean(),lon_0=lon.mean(),resolution='i')
# compute map projection coordinates of grid.
lons, lats = m.makegrid(RasterYSize+1, RasterXSize+1)
x, y = m(lons, lats)

x, y = m(lon, lat)

S_VV_ABS_new = m.transform_scalar(S_VV_ABS,lon[0,:],lat[0,:],RasterYSize,RasterXSize)


# create figure and axes instances
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    m.fillcontinents(color='coral',lake_color='aqua')
    m.drawmapboundary()
    m.drawcoastlines(linewidth=0.25,color='k')
    m.drawcountries(linewidth=0.25,color='k')
    m.drawmeridians(arange(1,4,0.5),labels=[1,0,0,1],fontsize=12, color=overlay_color)
    m.drawparallels(arange(41,42,0.5),labels=[0,1,1,0],fontsize=12, color=overlay_color)
    
    cs = m.pcolormesh(x,y,S_VV_ABS,shading='flat',cmap=plt.cm.gist_gray, vmin=-0.1, vmax=0.1)
    # add colorbar.
    cbar = m.colorbar(cs,location='bottom',pad="5%")
    cbar.set_label('SigmaVVwnr [linear units]')
    









lats = numpy.resize(lat, (round(RasterXSize/10), round(RasterYSize/10)) )
lons = numpy.resize(lon, (round(RasterXSize/10), round(RasterYSize/10)) )


lons, lats = m.makegrid(RasterXSize, RasterYSize) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats)


    
# setup Lambert Conformal basemap.
# set resolution=None to skip processing of boundary datasets.
m = Basemap(width=12000000,height=9000000,projection='lcc',
            resolution=None,lat_1=45.,lat_2=55,lat_0=50,lon_0=-107.)

m = Basemap(projection='ortho', resolution='None',
            lat_1=min(lat),lat_2=max(lat),lon_1=min(lon),lon_2=max(lon))
m = Basemap(projection='stere', resolution='l',
            lat_0=(max(lat)+min(lat))/2,lon_0=(max(lon)+min(lon))/2)
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='coral',lake_color='aqua')

m.shadedrelief()
plt.savefig('background4.png')


lats = lat
lons = lon
# shift lats, lons so values represent edges of grid boxes
# (as pcolor expects).
delon = lons[1]-lons[0]; delat = lats[1]-lats[0]
lons = (lons - 0.5*delon).tolist()
lons.append(lons[-1]+delon)
lons = array(lons,float64)
lats = (lats - 0.5*delat).tolist()
lats.append(lats[-1]+delat)
lats = array(lats,float64)

# creat figure, axes instances.
fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
# create Basemap instance for Robinson projection.
# coastlines not used, so resolution set to None to skip
# continent processing (this speeds things up a bit)
m = Basemap(projection='robin',lon_0=lons.mean(),resolution=None)

# compute map projection coordinates of grid.
x, y = m(*meshgrid(lons, lats))
# draw line around map projection limb.
# color background of map projection region.
# missing values over land will show up this color.
m.drawmapboundary(fill_color='0.3')
# plot data with pcolor
#im1 = m.pcolor(x,y,SigmaVVmf,shading='flat',cmap=plt.cm.jet)
im2 = m.pcolor(x,y,SigmaVVmf,shading='flat',cmap=plt.cm.gist_gray)
# draw parallels and meridians, but don't bother labelling them.
m.drawparallels(arange(-90.,120.,30.))
m.drawmeridians(arange(0.,420.,60.))
# add colorbar
cb = m.colorbar(im2,"bottom", size="5%", pad="2%")
# add a title.
ax.set_title('SST and ICE analysis for %s'%date)
plt.savefig('plotsst.png')
