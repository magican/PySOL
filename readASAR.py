# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 21:54:28 2012

@author: mag
"""

"""
Remarks to nansat
when "print n" units are "m/m"

"""

#~ import matplotlib.pyplot as plt
#~ from scipy.io import savemat

import ipdb

import Image

from nansat import Nansat
from domain import Domain

import distancelib
from numpy import deg2rad, array, log10, sin, exp, log

import os, sys, getopt

import gc 
#~ Actually calling gc.collect() yourself at the end of a loop can help avoid fragmenting memory, which in turn helps keep performance up. I've seen this make a significant difference (~20% runtime IIRC)

#~ oPath = '/home/mag/Documents/repos/solab/'

#~ iPath = '/media/SOLabNFS2/store/satellite/asar/2012/001/'
#~ fileName = 'ASA_WSM_1PNPDE20120101_204124_000003303110_00158_51465_5541.N1'

#~ iPath = '/media/data/data/baltic/finngulfWindCases/'
#~ fileName = 'ASA_WSM_1PNPDE20111127_085336_000002143109_00079_50955_3727.N1'

#~ iPath = '/media/data/data/baltic/finngulfWindCases/'
#~ fileName = 'ASA_WSM_1PNPDE20110523_084634_000000983102_00395_48254_2349.N1'

def main( argv=None ):

    year = '2012'
    useMask = False

    if argv is None:
        argv = sys.argv

    if argv is None:
        print ( "Please specify the path/year to the asar folder! \n")
        return

    # Parse arguments
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["year=","oPath=","iPath=","useMask="])
    except getopt.GetoptError:
        print 'readASAR.py -year <year> ...'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'readASAR.py -year <year> ...'
            sys.exit()
        elif opt in ("-year", "--year"):
            year = arg
        elif opt in ("-oPath", "--oPath"):
            oPath = arg
        elif opt in ("-iPath", "--iPath"):
            iPath = arg
        elif opt in ("-useMask", "--useMask"):
            useMask = arg

    oPath = '/media/SOLabNFS2/tmp/roughness/' + year + '/'
    iPath = '/media/SOLabNFS2/store/satellite/asar/' + year + '/'

    if not os.path.exists(oPath):
        os.makedirs(oPath)

    dirNames=os.listdir(iPath)
    for dirName in dirNames:
        fileNames=os.listdir(iPath+dirName)
        for fileName in fileNames:
            figureName = oPath + fileName[0:27] + '/' + fileName + '_proj.png'
            kmlName = oPath + fileName[0:27] + '/' + fileName + '.kml'
            if not os.path.exists(oPath + fileName[0:27] + '/'):
                os.makedirs(oPath + fileName[0:27] + '/')

            if os.path.isfile(kmlName):
                print "%s already processed" % (fileName)
                continue
            else:
                print "%s" % (fileName)

            # try to create Nansat object
            try:
                n = Nansat(iPath + dirName + '/' + fileName, mapperName='asar', logLevel=27)
            except Exception as e:
                print "Failed to create Nansat object:"
                print str(e)
                os.rmdir(oPath + fileName[0:27] + '/' )
                continue
                

            #~ Get the bands
            raw_counts = n[1]
            inc_angle = n[2]

            #~ NICE image (roughness)
            pol = n.bands()[3]['polarization']
            if pol == 'HH':
                ph = (2.20495, -14.3561e-2, 11.28e-4)
                sigma0_hh_ref = exp( ( ph[0]+inc_angle*ph[1]+inc_angle**2*ph[2])*log(10) )
                roughness = n[3]/sigma0_hh_ref
            elif pol == 'VV':
                pv = (2.29373, -15.393e-2, 15.1762e-4)
                sigma0_vv_ref = exp( ( pv[0]+inc_angle*pv[1]+inc_angle**2*pv[2])*log(10) )
                roughness = n[3]/sigma0_vv_ref

            #~ Create new band
            n.add_band(bandID=4, array=roughness, \
               parameters={'name':'roughness', \
               'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave', \
               'dataType': 6})

            # Reproject image into Lat/Lon WGS84 (Simple Cylindrical) projection
            # 1. Cancel previous reprojection
            # 2. Get corners of the image and the pixel resolution
            # 3. Create Domain with stereographic projection, corner coordinates 1000m
            # 4. Reproject
            # 5. Write image
            n.reproject() # 1.
            lons, lats = n.get_corners() # 2.
            # Pixel resolution
            #~ pxlRes = distancelib.getPixelResolution(array(lats), array(lons), n.shape())
            #~ pxlRes = array(pxlRes)*360/40000 # great circle distance
            pxlRes = array(distancelib.getPixelResolution(array(lats), array(lons), n.shape(), 'deg'))
            
            
            ipdb.set_trace()
            
            
            if min(lats) >= 65 and max(lats) >= 75 and max(lats)-min(lats) >= 13:
               pxlRes = array([0.00065, 0.00065])*2 # make the resolution 150x150m
            #~ pxlRes = pxlRes*7 # make the resolution worser
            srsString = "+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs"
            #~ extentString = '-lle %f %f %f %f -ts 3000 3000' % (min(lons), min(lats), max(lons), max(lats))
            extentString = '-lle %f %f %f %f -tr %f %f' % (min(lons), min(lats), \
                            max(lons), max(lats), pxlRes[1], pxlRes[0])
            d = Domain(srs=srsString, ext=extentString) # 3.
            n.reproject(d) # 4.

            if useMask:
               # get array with watermask (landmask) b
               # it must be done after reprojection!
               # 1. Get Nansat object with watermask
               # 2. Get array from Nansat object. 0 - land, 1 - water
               #wm = n.watermask(mod44path='/media/magDesk/media/SOLabNFS/store/auxdata/coastline/mod44w/')
               wm = n.watermask(mod44path='/media/data/data/auxdata/coastline/mod44w/')
               wmArray = wm[1]

               #~ ОШИБКА numOfColor=255 не маскирует, потому что в figure.apply_mask: availIndeces = range(self.d['numOfColor'], 255 - 1)
               #~ n.write_figure(fileName=figureName, bands=[3], \
                            #~ numOfColor=255, mask_array=wmArray, mask_lut={0: 0},
                            #~ clim=[0,0.15], cmapName='gray', transparency=0) # 5.
               n.write_figure(fileName=figureName, bands=[4], \
                                 mask_array=wmArray, mask_lut={0: [0,0,0]},
                                 clim=[0,2], cmapName='gray', transparency=[0,0,0]) # 5.
            else:
               n.write_figure(fileName=figureName, bands=[1], \
                              clim=[0,2], cmapName='gray', transparency=[0,0,0]) # 5.

            # open the input image and convert to RGBA for further tiling with slbtiles
            input_img = Image.open(figureName)
            output_img = input_img.convert("RGBA")
            output_img.save(figureName)

            # make KML image
            n.write_kml_image(kmlFileName=kmlName, kmlFigureName=figureName)

            #~ Change the file permissions
            os.chmod(oPath, 0777)
            os.chmod(oPath + fileName[0:27] + '/', 0777)
            os.chmod(kmlName, 0777)
            os.chmod(figureName, 0777)

            #~ Change the owner and group
            #~ os.chown(oPath, 1111, 1111)
            #~ os.chown(oPath + fileName[0:27] + '/', 1111, 1111)
            #~ os.chown(kmlName, 1111, 1111)
            #~ os.chown(figureName, 1111, 1111)
            
            #~ garbage collection
            gc.collect()

if __name__ == "__main__":
    main(sys.argv[1:])
