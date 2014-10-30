#!/usr/bin/env python
# coding: utf-8
"""
Created on Fri Mar 23 17:05:57 2012
Creating this script with the help from
http://benjamindeschamps.ca/blog/2009/11/12/\
processing-radarsat-2-imagery-reading-raw-data-and-saving-rgb-composites/
@author: mag
"""
# pn = '/media/data/data/OTHER/RS2 Agulhas and Lion/RS2_FQA_1xQGSS20101218_173930_00000005/'
# pn = '/media/data/data/OTHER/RS2 Agulhas and Lion/RS2_SQA_1xQGSS20091224_164846_00000004/'
# import the required library
import sys, getopt
from numpy import rad2deg
from scipy.io import savemat
#import matplotlib.pyplot as plt

from calib import calib, incidence_angle

#global xOff, yOff, xS, yS, xBufScale, yBufScale, s, inpath

def Usage():
    print( "Usage: " + \
            "      " )
    return 1

def main( argv=None ):

    # default values
    xOff=0
    yOff=0
    xS=None
    yS=None
    xBufScale=None
    yBufScale=None
    s = 'abs'
    filter_name = 'wiener'
    ws = 7
    
    if argv is None:
        argv = sys.argv
        
    if argv is None:
        print ( "Please specify the path to the RS2 folder! \n" + \
                "See USAGE for more details \n")
        return Usage()
    
    # Parse arguments        
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["pn=","xoff="])
    except getopt.GetoptError:
        print 'readRS2.py -pn <inputfile> ...'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'readRS2.py -pn <inputfile> ...'
            sys.exit()
        elif opt in ("-pn", "--pn"):
            pn = arg
        elif opt in ("-xoff", "--xoff"):
            xOff = arg
        elif opt in ("-yoff", "--yoff"):
            yOff = arg
        elif opt in ("-xS", "--xS"):
            xS = arg
        elif opt in ("-yS", "--yS"):
            yS = arg
        elif opt in ("-xScale", "--xScale"):
            xBufScale = arg
        elif opt in ("-yScale", "--yScale"):
            yBufScale = arg
            
    print "Path to the RS2 file:\n", pn

    # calibrate the scene
    print "Calibrating"
    SigmaHH, SigmaVV, SigmaHV, SigmaVH = \
        calib(pn, xOff, yOff, xS, yS, xBufScale, yBufScale, s)
    
    # get the incidence angle
    print "Incidence Angle and Lats/Lons"
    IncidenceAngle = incidence_angle(pn)
    IncidenceAngle = rad2deg(IncidenceAngle)
#    
#    # get the lats/lons
#    (lat, lon, row, col) = calib.xml2geo()
#    
#    # filter the image using median filter
#    print "Filtering"
#    (SigmaHHmf, SigmaHVmf, SigmaVHmf, SigmaVVmf) = \
#        calib.speckle_filter('median', 7)
#    (SigmaHHmf, SigmaHVmf, SigmaVHmf, SigmaVVmf) = \
#        calib.speckle_filter('wiener', 7)
    
#    # save to MAT file
#    print "Exporting"
#    expFn = 'calibrated.mat'
#    savemat(pn + expFn, mdict={'SigmaHH':calib.SigmaHH, \
#        'SigmaVV':calib.SigmaVV, 'IncidenceAngle':IncidenceAngle, \
#        'lat':lat, 'lon':lon, 'SigmaHHmf':SigmaHHmf, 'SigmaHVmf':SigmaHVmf, \
#        'SigmaVHmf':SigmaVHmf, 'SigmaVVmf':SigmaVVmf, \
#        }, do_compression=True)
    
    # export to geotiff
#    calib.save_tiff()
    
    #plt.figure()
    #plt.imshow(SigmaHH_VVmf,  plt.cm.gray, vmin=-20, vmax=0)
    #plt.colorbar()
    #plt.clim(-5,7)
    #plt.figure()
    #plt.imshow(calib.SigmaHH,  plt.cm.gray, vmin=-20, vmax=0)
    #plt.figure()
    #plt.imshow(SigmaHHmf,  plt.cm.gray, vmin=-20, vmax=0)
    #plt.colorbar()
    #plt.show()
    
if __name__ == "__main__":
    main(sys.argv[1:])