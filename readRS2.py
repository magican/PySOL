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
# pn = '/media/SAMSUNG/Radarsat2/RS2 nersc/RS2_20120220_073807_0076_SCWA_HHHV_SGF_181759_2925_7176719/RS2_20120220_073807_0076_SCWA_HHHV_SGF_181759_2925_7176719/'
# pn = '/media/SOLabNFS2/store/satellite/RS2/RS2_FQA_1xQGSS20120729_144547_00000005xQ12_16bxx_24139_FC556/'

# import the required library
from numpy import rad2deg, log10, linspace, tile, conj, array, uint, \
                  real, imag, absolute, ones

import sys, getopt
from time import time
from os import sysconf, path, system
from scipy.signal import wiener
#from scipy.io import savemat

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from plotRS2 import plotRS2
from createMapsEtopo1 import makeMap

from g1sst import g1sst

from imcrop import imcrop, imzoom, intCrpLL

from skimage.filter import threshold_otsu
from skimage.filter import median_filter

#from cmod import rcs2wind

import datetime

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 5, 1)
__modified__ = datetime.datetime(2012, 5, 11)
__version__  = "1.0"
__status__   = "Development"

core_ct = sysconf('SC_NPROCESSORS_ONLN')
print "Running %s cores " % core_ct
if core_ct < 4:
    from calibRS2 import calibRS2
elif core_ct >= 4:
    from calibRS2par import calibRS2par
    print "Using parallel coding for image processing"


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

    if core_ct < 4:
        # calibrate the scene
        print "Calibrating"
        currtime = time()
        calib = calibRS2(pn, xOff, yOff, xS, yS, xBufScale, yBufScale, 'sigma0')
        # get the incidence angle
        print "Incidence Angle and Lats/Lons"
        IncidenceAngle = calib.incidence_angle()
        IncidenceAngle = rad2deg(IncidenceAngle)
        # get the lats/lons
        (lat, lon, line, pixel) = calib.xml2geo()
        print "Filtering"
        # 1. Remove speckle noise (using Lee-Wiener filter)
        calib.speckle_filter('wiener', 7)
        print 'Serial: time elapsed:', time() - currtime
    elif core_ct >= 4:
        # calibrate the scene
        print "Calibrating"
        currtime = time()
        calibPar = calibRS2par(pn, xOff, yOff, xS, yS, xBufScale, yBufScale, s)
        # get the incidence angle
        print "Incidence Angle and Lats/Lons"
        IncidenceAngle = calibPar.incidence_angle()
        IncidenceAngle = rad2deg(IncidenceAngle)
        # get the lats/lons
        (lat, lon, line, pixel) = calibPar.xml2geo()
        print "Filtering"
        # 1. Remove speckle noise (using Lee-Wiener filter)
        calibPar.speckle_filter('wiener', 7)
        print 'Parallel: time elapsed: %f' % ( time() - currtime )
        # export to geotiff
#        calibPar.save_tiff()

    print "Some Sigma calculations..."

    SigmaHHwnr = calibPar.SigmaHHwnr
    SigmaVVwnr = calibPar.SigmaVVwnr
#    SigmaHVwnr = calibPar.SigmaHVwnr

    # Free up some memory
    del calibPar.SigmaHHwnr, calibPar.SigmaVVwnr, calibPar.SigmaHVwnr

#    # Calculate some wind
#    # testing that with sar=-0.387 wind speed is 10m/s
#    w = rcs2wind(sar=-0.387*ones((1,1)), cmdv=4, windir=0, theta=20*ones((1,1)))
#    print "Testing CMOD4 passed, Wind =", w
#
#    sar=10*log10(SigmaVVwnr[::10,::10])
##    sar=10*log10(SigmaVVwnr)
#    w = rcs2wind(sar, cmdv=4, windir=130, theta=33*ones(sar.shape))
#    print "Mean CMOD4 Wind =", w.mean()

    # calculate the conjugate
    SigmaVVConjHH = calibPar.S_VV*conj(calibPar.S_HH)
    SigmaVVConjHHWnr = wiener(SigmaVVConjHH, mysize=(7,7), noise=None)
    # all sigmas in linear units
    SigmaVVConjHHdelta = SigmaVVConjHHWnr/(SigmaVVwnr - SigmaHHwnr)

#    # Check that imagenary part is much less then the real one
#    SigmaVVConjHHr = real(SigmaVVConjHH)
#    SigmaVVConjHHi = imag(SigmaVVConjHH)
#    SigmaVVConjHHr.mean()
#    SigmaVVConjHHi.mean()

    # Free up some memory
    del SigmaVVConjHH, calibPar.S_HH, calibPar.S_VV

    # 2. Find Poolarization ratio P in dB and delta - linear
    # delta - is a Bragg component, which describes wind field impact
    # P = 10*log10(SigmaHHwnr) - 10*log10(SigmaVVwnr)
#    # P must not be > 0, be careful with masking
    #P[P>0] = 0
    P = SigmaHHwnr/SigmaVVwnr
    delta = SigmaVVwnr - SigmaHHwnr



#    # Masking Land from HV polarization
#    maskedLand = SigmaHVwnr>0.003

    # 3. Linear fit of Sigma and delta: Sigma = A*delta
#    A = sum(SigmaVVwnr[~maskedLand]*delta[~maskedLand])/sum(delta[~maskedLand]^2);
#    B = sum(SigmaHHwnr[~maskedLand]*delta[~maskedLand])/sum(delta[~maskedLand]^2);
    # imposing from theory
    A = linspace(2.3,1.8,IncidenceAngle.size) # dependence on inc angle
    A = tile(A, (line[-1,0]+1, 1))

    # 4. Calculate SigmaVVwb = SigmaVVwb === SigmaHHwb
    # substracting wind field variability contribution, the rest is impact of
    # current and film on the image
    SigmaWb = SigmaVVwnr - ( A*delta )

    # Checking that SigmaWb from VV and HH are the same
    # SigmaVVwb = SigmaVVwnr - ( A*delta )
    # SigmaHHwb = SigmaHHwnr - ( B*delta )
    # SigmaWb = SigmaVVwb
    # mean(SigmaHHwb(:)./SigmaVVwb(:))
    # clear SigmaVVwb SigmaHHwb

    # Cropping the oil slicks

#    проверить почему переворот изображения!!!
#    print ( "Cropping...")
#    plt.figure()
#    plt.imshow(delta)
#    plt.clim(-0.01,0.03)
#    plt.gray()
#    plt.show()
#    pts = imzoom()

    ptsCropSouth = array([[ 3200, 4550],
                          [ 3800, 5100]])
    oilCropSouth = imcrop(ptsCropSouth, delta)
    ptsCropNorth = array([[ 3400, 0],
                          [ 4900,  1650]])
    oilCropNorth = imcrop(ptsCropNorth, delta)

    # reduce the speckle noise
    oilN = wiener(oilCropNorth, mysize=(3,3), noise=None)
    oilS = wiener(oilCropSouth, mysize=(3,3), noise=None)

    # use median filter to smooth
    oilN = median_filter(oilN, radius=37)
    oilS = median_filter(oilS, radius=37)

    # find the global threshold
    global_threshN = threshold_otsu(oilN)
    global_threshS = threshold_otsu(oilS)

    # create a binary mask
    markersN = ones(oilN.shape, dtype=uint)
    markersS = ones(oilS.shape, dtype=uint)
    markersN[oilN < global_threshN*0.84] = 0
    markersS[oilS < global_threshS*0.84] = 0

    # Find the relation of P in/out-side the Slicks
    oilCropNorthP = imcrop(ptsCropNorth, P)
    oilCropPN = oilCropNorthP[markersN==0].mean()/oilCropNorthP[markersS==1].mean()

    oilCropSouthP = imcrop(ptsCropSouth, P)
    oilCropPS = oilCropSouthP[markersS==0].mean()/oilCropSouthP[markersS==1].mean()

    # Checking that calculations are OK
    oilCropNorthKVV = imcrop(ptsCropNorth, SigmaVVwnr)
    KVV = oilCropNorthKVV[markersN==0].mean()/oilCropNorthKVV[markersN==1].mean()
    oilCropNorthKHH = imcrop(ptsCropNorth, SigmaHHwnr)
    KHH = oilCropNorthKHH[markersN==0].mean()/oilCropNorthKHH[markersN==1].mean()
    PSlick = KHH/KVV*oilCropNorthKHH[markersN==1].mean()/oilCropNorthKVV[markersN==1].mean()

    print "Polarization ratio Out/In-side Slick = \
        \n %0.2f for Northern \n %0.2f for Southern" \
        %(1/oilCropPN, 1/oilCropPS)

    print ( "Plotting on a map... \nUsing default NSPER Basemap \n")

    # setup nsper basemap
    # Lat/Lon coords of image corners
    ll_lat = lat.min()
    ur_lat = lat.max()
    ll_lon = lon.min()
    ur_lon = lon.max()
    cent_lat = lat.mean()
    cent_lon = lon.mean()
    m = Basemap(llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
                llcrnrlon=ll_lon, urcrnrlon=ur_lon, \
                resolution='f', projection='nsper', \
                satellite_height=798000, \
                lat_0=cent_lat,lon_0=cent_lon)

#    m = Basemap(llcrnrlat=30, urcrnrlat=50,\
#                llcrnrlon=1, urcrnrlon=20, \
#                resolution='l', projection='nsper', \
#                satellite_height=798000, \
#                lat_0=cent_lat,lon_0=cent_lon)

    print ( "Plotting figures... \n")

    scale = 8 # setting scale factor to resize in 1/scale times
    label1 = 'SigmaVVwnr [dB]'
    label2 = 'SigmaHHwnr [dB]'
    label11 = 'SigmaVVwnr [linear units]'
    label22 = 'SigmaHHwnr [linear units]'
    label3 = 'P [dB]'
    label33 = 'P [linear units]'
    label4 = 'delta [linear units]'
    label5 = 'SigmaWb [linear units]'
    label6 = 'SigmaVVConjHH [linear units]'
    label7 = 'SigmaVVConjHHdelta [linear units]'
    label8 = 'Slick delta [linear units]'
    label88 = 'Slick P [linear units]'
    label9 = 'Slick delta [linear units]'
    label99 = 'Slick P [linear units]'

#    label8 = 'Wind CMOD4 [m/s]'

    import plotRS2
    reload(plotRS2)
    from plotRS2 import plotRS2

#    import mpl_util
#    reload(mpl_util)
#
#    import gmtColormap
#    reload(gmtColormap)
#
#    import createMapsEtopo1
#    reload(createMapsEtopo1)
#    from createMapsEtopo1 import makeMap
#
#    import g1sst
#    reload(g1sst)
#    from g1sst import g1sst

#    plotRS2(10*log10(SigmaVVwnr), lat, lon, pixel, line, \
#            scale=scale, m=m, clm=(-25,-5), label=label1)
#    makeMap(ll_lon, ur_lon, \
#            ll_lat, ur_lat, m, name=label1, contour='land')
#    plt.close('all')
#
#    plotRS2(10*log10(SigmaHHwnr), lat, lon, pixel, line, \
#            scale=scale, m=m, clm=(-25,-5), label=label2)
#    makeMap(ll_lon, ur_lon, \
#            ll_lat, ur_lat, m, name=label2, contour='land')
#    plt.close('all')

    plotRS2(SigmaVVwnr, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(0,0.05), label=label11)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label11, contour='land')
    plt.close('all')

    plotRS2(SigmaHHwnr, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(0,0.05), label=label22)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label22, contour='land')
    plt.close('all')

#    plotRS2(P, lat, lon, pixel, line, \
#            scale=scale, m=m, clm=(-2.7,0), label=label3)
#    makeMap(ll_lon, ur_lon, \
#            ll_lat, ur_lat, m, name=label3, contour='land')
#    plt.close('all')

    plotRS2(P, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(0.3,1.3), label=label33)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label3, contour='land')
    plt.close('all')

    plotRS2(delta, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(-0.01,0.03), label=label4)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label4, contour='land')
    plt.close('all')

    plotRS2(SigmaWb, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(0,0.03), label=label5)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label5, contour='land')
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label5, contour='ocean')
    plt.close('all')

    plotRS2(SigmaVVConjHHWnr, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(0,0.05), label=label6)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label6, contour='land')
    plt.close('all')

    plotRS2(SigmaVVConjHHdelta, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(0,5), label=label7)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label7, contour='land')
    plt.close('all')

    # Plot SST on a new figure
    plt.figure()
    g1sst(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name='g1sst', contour=None)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name='g1sst', contour='land')
    plt.close('all')

    # Plot cropped Oil Spills
    # Northern
    # setup nsper basemap
    # Lat/Lon coords of image corners
    lat_new, lon_new = intCrpLL(lat, lon, pixel, line, scale=1, ptsCrop=ptsCropNorth)
    ll_lat = lat_new.min()
    ur_lat = lat_new.max()
    ll_lon = lon_new.min()
    ur_lon = lon_new.max()
    cent_lat = lat_new.mean()
    cent_lon = lon_new.mean()
    m = Basemap(llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
                llcrnrlon=ll_lon, urcrnrlon=ur_lon, \
                resolution='c', projection='nsper', \
                satellite_height=798000, \
                lat_0=cent_lat,lon_0=cent_lon)
    # Delta
    plotRS2(delta, lat, lon, pixel, line, \
            scale=1, m=m, clm=(-0.01,0.03), label=label8, ptsCrop=ptsCropNorth)
    plt.close('all')
    system("mv " + "/home/mag/Slick.png" + " " + "/home/mag/NorthernDelta.png")
    # P
    plotRS2(P, lat, lon, pixel, line, \
            scale=1, m=m, clm=(0.3,1.3), label=label88, ptsCrop=ptsCropNorth)
    plt.close('all')
    system("mv " + "/home/mag/Slick.png" + " " + "/home/mag/NorthernP.png")
    # Southern
        # setup nsper basemap
    # Lat/Lon coords of image corners
    lat_new, lon_new = intCrpLL(lat, lon, pixel, line, scale=1, ptsCrop=ptsCropSouth)
    ll_lat = lat_new.min()
    ur_lat = lat_new.max()
    ll_lon = lon_new.min()
    ur_lon = lon_new.max()
    cent_lat = lat_new.mean()
    cent_lon = lon_new.mean()
    m = Basemap(llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
                llcrnrlon=ll_lon, urcrnrlon=ur_lon, \
                resolution='c', projection='nsper', \
                satellite_height=798000, \
                lat_0=cent_lat,lon_0=cent_lon)
    # Delta
    plotRS2(delta, lat, lon, pixel, line, \
            scale=1, m=m, clm=(-0.01,0.03), label=label9, ptsCrop=ptsCropSouth)
    plt.close('all')
    system("mv " + "/home/mag/Slick.png" + " " + "/home/mag/SouthernDelta.png")
    # P
    plotRS2(P, lat, lon, pixel, line, \
            scale=1, m=m, clm=(0.3,1.3), label=label99, ptsCrop=ptsCropSouth)
    plt.close('all')
    system("mv " + "/home/mag/Slick.png" + " " + "/home/mag/SouthernP.png")



#    lat_new, lon_new = intCrpLL(lat, lon, pixel, line, scale=1, ptsCrop=ptsCropNorth)
#    ll_lat = lat_new.min()
#    ur_lat = lat_new.max()
#    ll_lon = lon_new.min()
#    ur_lon = lon_new.max()
#    makeMap(ll_lon, ur_lon, ll_lat, ur_lat, name=label8)


#    # Plot RS2 CMOD4 Wind
#    plotRS2(w, lat, lon, pixel, line, \
#            scale=1, m=m, clm=(0,7), label=label8)
#    makeMap(ll_lon, ur_lon, \
#            ll_lat, ur_lat, m, name=label8, contour='land')
#    plt.close('all')

#    # save to MAT file
#    print "Exporting"
#    expFn = 'calibrated.mat'
#    if not path.isfile(pn + expFn):
#        print "Exporting to:\n" , pn + expFn
#        savemat(pn + expFn, mdict={ \
#                'SigmaHH':calibPar.SigmaHH, \
#                'SigmaVV':calibPar.SigmaVV, \
#                'IncidenceAngle':IncidenceAngle, \
#                'lat':lat, 'lon':lon, \
#                'SigmaHHwnr':SigmaHHwnr, 'SigmaVVwnr':SigmaVVwnr,\
#        #        'SigmaHVwnr':SigmaHVwnr, 'SigmaVHwnr':SigmaVHwnr, \
#                'S_HH':calibPar.S_HH, 'S_VV':calibPar.S_VV
#                }, do_compression=False)

#    # calculate the magnitude
#    S_VV_ABS = absolute(calibPar.S_VV)
#    S_HH_ABS = absolute(calibPar.S_HH)
#    # calculate the linear sigma_naught
#    SigmaVV = S_VV_ABS**2
#    SigmaHH = S_HH_ABS**2
#
#    SigmaVV_HH = calibPar.SigmaVV - calibPar.SigmaHH
#    SigmaVVConjHH = calibPar.S_VV*conj(calibPar.S_HH)
#    SigmaVVConjHH = real(SigmaVVConjHH)
#
#    # save to MAT file
#    print "Exporting"
#    expFn = '/home/mag/Public/calibrated.mat'
#
#    if not path.isfile(expFn):
#        print "Exporting to:\n" , expFn
#        savemat(expFn, mdict={ \
#                'SigmaVVConjHH':SigmaVVConjHH, \
#                'SigmaVV_HH':SigmaVV_HH,
#                'IncidenceAngle':IncidenceAngle, \
#                'lat':lat, 'lon':lon, \
#                }, do_compression=False)
#
#    expFn2 = '/home/mag/Public/calibrated_HHVV.mat'
#    if not path.isfile(expFn2):
#        print "Exporting to:\n" , expFn2
#        savemat(expFn2, mdict={ \
#                'SigmaVV':SigmaVV, \
#                'SigmaHH':SigmaHH,
#                'IncidenceAngle':IncidenceAngle, \
#                'lat':lat, 'lon':lon, \
#                }, do_compression=False)

if __name__ == "__main__":
    main(sys.argv[1:])