#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 18:08:52 2012

@author: mag
"""

import datetime

__author__   = 'Alexander Myasoedov'
__email__= 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 9, 10)
__modified__ = datetime.datetime(2012, 10, 18)
__version__  = "1.0"
__status__   = "Development"

import gdal
import glob
from scipy.misc import imresize
from numpy import float64, mean
from matplotlib.mlab import find
from pylab import imshow

def Usage():
    print( "Usage: " + \
            "      " )
    return 1

def main( argv=None ):

    # default values
    pn = '/media/SOLabNFS2/store/satellite/modis/blacksea/allData/5/MYD02QKM/2012/258/'
    fn = 'MYD02QKM.A2012258.1035.005.2012259155603.hdf'

    if argv is None:
        argv = sys.argv

    if argv is None:
        print ( "Please specify the path to the RS2 folder! \n" + \
                "See USAGE for more details \n")
        return Usage()

    # Parse arguments
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["pn=","fn="])
    except getopt.GetoptError:
        print 'readRS2.py -pn <inputfile> ...'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'readRS2.py -pn <inputfile> ...'
            sys.exit()
        elif opt in ("-pn", "--pn"):
            pn = arg
        elif opt in ("-fn", "--fn"):
            fn = arg

    # To check which fields exist
    #ds=gdal.Open(pn + fn)
    #sds_md = ds.GetMetadata('SUBDATASETS')
    #print 'keys'
    #for i in sds_md.iterkeys():
    #    print i,sds_md[i]


    # Reading 250m resolution bands to G1, G2
    mod = gdal.Open('HDF4_EOS:EOS_SWATH:"' + pn + fn + '":MODIS_SWATH_Type_L1B:EV_250_RefSB')
    B = mod.GetRasterBand(1).ReadAsArray()

    # Checking if other supplementary files exist
    pn35_L2 = pn[0:-15] + '35_L2' + pn[-10:]
    fn35_L2 = fn[0:-41] + '35_L2' + fn[-36:-17] + '*'
    check = str(glob.glob(pn35_L2 + fn35_L2))
    if len(check)<=2:
        # check in the same dir
        check = str(glob.glob(pn + fn35_L2))
        if len(check)<=2:
            print "No '35_L2' exist"
            fn35_L2 = []
            pn35_L2 = []
        else: fn35_L2 = check[-46:-2]; pn35_L2 = pn
    else: fn35_L2 = check[-46:-2]

    pn03 = pn[0:-15] + '03' + pn[-10:]
    fn03 = fn[0:-40] + '3' + fn[-36:-17] + '*'
    check = str(glob.glob(pn03 + fn03))
    if len(check)<=2:
        # check in the same dir
        check = str(glob.glob(pn + fn03))
        if len(check)<=2:
            print "No '03' exist"
            fn03 = []
            pn03 = []
        else: fn03 = check[-46:-2]; pn03 = pn
    else: fn03 = check[-43:-2]

    # Reading geolocation
    if fn03:
        print "   ...Using MYD03/MOD03 geolocation"
        Latitude  = gdal.Open('HDF4_SDS:UNKNOWN:"' + pn03 + fn03 + '":0').ReadAsArray()
        Longitude = gdal.Open('HDF4_SDS:UNKNOWN:"' + pn03 + fn03 + '":1').ReadAsArray()
    else:
        print "   ...Using internal file geolocation"
        Latitude  = gdal.Open('HDF4_SDS:UNKNOWN:"' + pn + fn + '":0').ReadAsArray()
        Longitude = gdal.Open('HDF4_SDS:UNKNOWN:"' + pn + fn + '":1').ReadAsArray()

    Latitude  = imresize(Latitude, B.shape)
    Longitude = imresize(Longitude, B.shape)

    # Masking
    if fn03:
        print "   ...Using MYD035/MOD035 land/cloud mask"
        Masks  = gdal.Open('HDF4_SDS:UNKNOWN:"' + pn35_L2 + fn35_L2 + '":7')
        indGlitter  = Masks.GetRasterBand(1).ReadAsArray()
        indGlitter = imresize(indGlitter, B.shape)
        indGlitter = indGlitter == 47
        indLandOrig = Masks.GetRasterBand(1).ReadAsArray()
        indLandOrig = imresize(indLandOrig, B.shape)
        indLandOrig = indLandOrig > 127
        indCloudOrig = Masks.GetRasterBand(6).ReadAsArray()
        indCloudOrig = imresize(indCloudOrig, B.shape)
        indCloudOrig = indCloudOrig == 0
    else:
        print "   ...Using mod44w land mask"
#        pnmod44w = '/media/SOLabNFS/store/auxdata/coastline/mod44w/'
#        fnmod44w = 'MOD44W_Water_2000_KJ3334_lzw.tif'
#        mod44w = gdal.Open(pnmod44w + fnmod44w).ReadAsArray()
#        Need to reproject

    # Enlarge the land/cloud masks


    # Crtopping the glitter part only
    indGlitter = find( (mean(indGlitter,axis=0)) > 0.01 )
    indGlitter = xrange(min(indGlitter), max(indGlitter))
    indCloudOrig = indCloudOrig[:,indGlitter]
    indLandOrig  = indLandOrig[:,indGlitter]

    maskedAll = float64(B[:,indGlitter])
    maskedAll[indCloudOrig] = float64('nan')
    maskedAll[indLandOrig] = float64('nan')

    imshow(maskedAll,vmin=0,vmax=1300,interpolation='nearest',origin='upper')
