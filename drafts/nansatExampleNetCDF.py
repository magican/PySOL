# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 21:54:28 2012

@author: mag
"""

"""
Remarks to nansat
when "print n" units are "m/m"

"""

from nansat import Nansat
from domain import Domain
import matplotlib.pyplot as plt
from scipy.io import savemat

import distancelib
from numpy import array

iPath = '/media/data/data/gulfstream/Gulfstream_20100401_synergy/'
oPath = '/home/mag/'
fileName = 'MER_FRS_1PNEPA20100401_154249_000001972088_00140_42278_0585.N1.nc'

# create Nansat object
n = Nansat(iPath + fileName)

# list bands and georeference of the object
print n

# Reprojected image into Lat/Lon WGS84 (Simple Cylindrical) projection
# 1. Cancel previous reprojection
# 2. Get corners of the image and the pixel resolution
# 3. Create Domain with stereographic projection, corner coordinates and resolution 1000m
# 4. Reproject
# 5. Write image
n.reproject() # 1.
lons, lats = n.get_corners() # 2.
pxlRes = distancelib.getPixelResolution(array(lats), array(lons), n[1])
pxlRes = array(pxlRes)*360/40000 # great circle distance
srsString = "+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs"
#~ extentString = '-lle %f %f %f %f -ts 3000 3000' % (min(lons), min(lats), max(lons), max(lats))
extentString = '-lle %f %f %f %f -tr %f %f' % (min(lons), min(lats), \
                max(lons), max(lats), pxlRes[1], pxlRes[0])
d = Domain(srs=srsString, ext=extentString) # 3.
print d
n.reproject(d) # 4.

# get array with watermask (landmask) b 
# it must be done after reprojection!
# 1. Get Nansat object with watermask
# 2. Get array from Nansat object. 0 - land, 1 - water
#wm = n.watermask(mod44path='/media/magDesk/media/SOLabNFS/store/auxdata/coastline/mod44w/')
wm = n.watermask(mod44path='/media/data/data/auxdata/coastline/mod44w/')
wmArray = wm[1]

figureName = oPath + fileName + '_proj.png'
n.write_figure(fileName=figureName, clim=[3,133], \
               mask_array=wmArray, mask_lut={0: [204, 153, 25]}) # 5.

#~ # make KML file with image borders (to be opened in Googe Earth)
#~ n.write_kml(kmlFileName=oPath + fileName + '_preview.kml')

# make KML image file with image borders (to be opened in Googe Earth)
n.write_image_kml(kmlFileName=oPath + fileName + '.kml', kmlFigureName=figureName)
