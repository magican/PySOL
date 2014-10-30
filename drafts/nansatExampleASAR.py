# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 17:05:57 2012

@author: mag
"""

#import numpy
#from pylab import *
#from pylab import imshow
import matplotlib.pyplot as plt

pn = '/home/mag/data/baltic/finngulfWindCases/'
fn = 'ASA_WSM_1PNPDE20110523_084634_000000983102_00395_48254_2349.N1'
fn = 'ASA_WSM_1PNPDE20111127_085336_000002143109_00079_50955_3727.N1'

import epr

product = epr.open(pn + fn)

proc_dataset = product.get_dataset('MAIN_PROCESSING_PARAMS_ADS')

record = proc_dataset.read_record(0)
record.get_field_names()[:20]

field = record.get_field('range_spacing')

product.get_band_names()


# Get proc_data
pd = product.get_band('proc_data')
proc_data = pd.read_as_array(5000, 5000, xoffset=500, yoffset=9000, \
                             xstep=2, ystep=2)

# Get lat/lon
lat = product.get_band('latitude')
lat = lat.read_as_array(5000, 5000, xoffset=500, yoffset=9000, \
                        xstep=2, ystep=2)
lon = product.get_band('longitude')
lon = lon.read_as_array(5000, 5000, xoffset=500, yoffset=9000, \
                        xstep=2, ystep=2)



# Get incident_angle
ia = product.get_band('incident_angle')
incident_angle = ia.read_as_array(5000, 5000, xoffset=500, yoffset=9000, \
                                  xstep=2, ystep=2)


proc_dataset = product.get_dataset('GEOLOCATION_GRID_ADS')
record = proc_dataset.read_record(0)
record.get_field_names()[:20]
angles = record.get_field('last_line_tie_points.angles')
angles = angles.get_elems()

# Calibrate



###------------------------------------------------------------------------###
# USING GDAL

import gdal, struct
from scipy import array,empty,uint16
from pylab import imshow

f=gdal.Open(pn + fn)
a=f.GetRasterBand(1)
img=empty((a.YSize,a.XSize),dtype=uint16)

for yi in xrange(a.YSize):
    scanline=struct.unpack("H"*a.XSize,a.ReadRaster(0,yi,a.XSize,1,a.XSize,1,gdal.GDT_UInt16))
    img[yi,:]=array(scanline).astype(uint16)

imshow(img,vmin=0,vmax=3000,interpolation='nearest',origin='lower')



###------------------------------------------------------------------------###
# USING nansat
from nansat import Nansat
from domain import Domain
import matplotlib.pyplot as plt
from scipy.io import savemat
from numpy import deg2rad

import distancelib
from numpy import array

iPath = '/home/mag/data/baltic/finngulfWindCases/'
oPath = '/home/mag/'
fileName = 'ASA_WSM_1PNPDE20110523_084634_000000983102_00395_48254_2349.N1'

oFileName = oPath + fileName

# create Nansat object
n = Nansat(iPath + fileName, mapperName='ASAR')
#n = Nansat(iPath + fileName)

# list bands and georeference of the object
print n

# get dictionary with all bands metadata
print n.bands()

# get size of the object (Y and X dimensions, to follow Numpy style)
print n.shape()

# get list with coordinates of the object corners
print n.get_corners()

# get lists with coordinates of the object borders
print n.get_border()

raw_counts = n[1]
inc_angle = n[2]

#~ sigma0 = n[3]

sigma0 = raw_counts**2.0 * sin(deg2rad(inc_angle))
sigma0 = 10*log10(sigma0)
n.add_band(bandID=4, array=sigma0)

# 1. Remove speckle noise (using Lee-Wiener filter)
speckle_filter('wiener', 7)

# Reprojected image into Lat/Lon WGS84 (Simple Cylindrical) projection
# 1. Cancel previous reprojection
# 2. Get corners of the image and the pixel resolution
# 3. Create Domain with stereographic projection, corner coordinates and resolution 1000m
# 4. Reproject
# 5. Write image
n.reproject() # 1.
lons, lats = n.get_corners() # 2.
pxlRes = distancelib.getPixelResolution(array(lats), array(lons), n.shape(), units="deg")
srsString = "+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs"
#~ extentString = '-lle %f %f %f %f -ts 3000 3000' % (min(lons), min(lats), max(lons), max(lats))
extentString = '-lle %f %f %f %f -tr %f %f' % (min(lons), min(lats), \
                max(lons), max(lats), pxlRes[1], pxlRes[0])
d = Domain(srs=srsString, ext=extentString) # 3.
n.reproject(d) # 4.

# get array with watermask (landmask) b 
# it must be done after reprojection!
# 1. Get Nansat object with watermask
# 2. Get array from Nansat object. 0 - land, 1 - water
#wm = n.watermask(mod44path='/media/magDesk/media/SOLabNFS/store/auxdata/coastline/mod44w/')
wm = n.watermask(mod44path='/media/data/data/auxdata/coastline/mod44w/')
wmArray = wm[1]

n.write_figure(oFileName + '_proj.png', clim='hist', bands=[3], \
               mask_array=wmArray, mask_lut={0: [204, 153, 25]}) # 5.

# write figure with lat/lon grid
# 1. Get lat/lon arrays from Nansat object (may take some time)
# 2. Make figure with lat/lon grids
lonGrid, latGrid = n.get_geolocation_grids()
n.write_figure(oFileName + '_latlon.png', bands=[1], \
               latGrid=latGrid, lonGrid=lonGrid, \
               latlonGridSpacing=10, latlonLabels=10, \
               mask_array=wmArray, mask_lut={0: [0, 0, 0]})

# make KML file with image borders (to be opened in Googe Earth)
n.write_kml(kmlFileName=oFileName + '_preview.kml')

# make image with map of the file location
n.write_map(oFileName + '_map.png')

# Reprojected image into stereographic projection
# 1. Cancel previous reprojection
# 2. Get corners of the image
# 3. Create Domain with stereographic projection, corner coordinates and resolution 1000m
# 4. Reproject
# 5. Write image
n.reproject() # 1.
lons, lats = n.get_corners() # 2.
meanLon = sum(lons, 0.0) / 4.
meanLat = sum(lats, 0.0) / 4.
srsString = "+proj=stere +lon_0=%f +lat_0=%f +k=1 +ellps=WGS84 +datum=WGS84 +no_defs" % (meanLon, meanLat)
extentString = '-lle %f %f %f %f -tr 75 75' % (min(lons), min(lats), max(lons), max(lats))
dStereo = Domain(srs=srsString, ext=extentString) # 3.
dStereo.write_map(oFileName + '_stereo_map.png')
print dStereo
n.reproject(dStereo) # 4.

n.write_figure(oFileName + '_pro_stereo.png', bands=[3], clim=[-20,0]) # 5.


# Make image reprojected onto map of Northern Europe
# 1. Create Domain object. It describes the desired grid of reprojected image:
# projection, resolution, size, etc. In this case it is geographic projection;
# -10 - 30 E, 50 - 70 W; 2000 x 2000 pixels
# 2. Reproject the Nansat object
# 3. Make simple image
srsString = "+proj=latlong +ellps=WGS84 +no_defs"
extentString = '-lle %f %f %f %f -ts 1000 1000' % (min(lons), min(lats), max(lons), max(lats))

dLatlong = Domain(srs=srsString, ext=extentString)
dLatlong.write_map(oFileName + '_latlong_map.png')
print dLatlong
n.reproject(dLatlong)

n.write_figure(oFileName + '_pro_latlon.tif')













lonVec, latVec = n.get_corners()

lonmin = min(lonVec)
lonmax = max(lonVec)
latmin = min(latVec)
latmax = max(latVec)

lonmin, latmin, lonmax, latmax





d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs", \
           "-lle 20.01092086474455, 54.46124717815684, 33.75349240659604, 67.79294542296294 -ts 2000 2000")

print d

n.reproject(d)


n.resize(0.1)
n.write_geotiffimage(oPath + fileName + '_pro_latlon.tif', bandID=1)
n.resize()

#n.write_figure(oPath + fileName + '_preview.png', clim='hist', ratio=0.95, cmapName='gray')







"""
Figures
"""

row = proc_data.shape[0]`
col = proc_data.shape[1]

fig = plt.figure(num=33)
#plt.imshow(proc_data, cmap=plt.cm.gist_gray, vmin=0, vmax=3000)
plt.imshow(proc_data, cmap=plt.cm.bone, vmin=0, vmax=3000)
plt.title(pd.description)
plt.colorbar()
#plt.show()

# Now check everything with the defaults:
res = fig.get_dpi() # default res is 80dpi
print "res:", res
defaultSize = fig.get_size_inches()
print "Default size in Inches", defaultSize
print "Which should result in a %i x %i Image"%(res*defaultSize[0], res*defaultSize[1])
# the default is 80dpi for savefig:
fig.savefig("test1.png")

# Now make the image twice as big, while keeping the fonts and all the
# same size
fig.set_size_inches( defaultSize*2 )
Size = fig.get_size_inches()
print "Size in Inches", Size
print "Which should result in a %i x %i Image"%(res*Size[0], res*Size[1])
fig.savefig("test2.png")

fig.savefig('test.png', figsize=Size, dpi=res, facecolor='w', edgecolor='k', bbox_inches='tight', pad_inches=0)
