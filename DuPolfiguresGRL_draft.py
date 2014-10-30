#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  9 13:27:06 2012
Created on Fri Mar 23 17:05:57 2012
Creating this script with the help from
http://benjamindeschamps.ca/blog/2009/11/12/\
processing-radarsat-2-imagery-reading-raw-data-and-saving-rgb-composites/
@author: mag
"""

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

from imcrop import imcrop, intCrpLL

import transect

from skimage.filter import threshold_otsu
from skimage.filter import median_filter

import gtk

#from cmod import rcs2wind

import datetime

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 5, 1)
__modified__ = datetime.datetime(2012, 6, 9)
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
    pn = '/media/data/data/OTHER/RS2 Agulhas and Lion/RS2_FQA_1xQGSS20101218_173930_00000005/'
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

    # 2. Find Poolarization ratio PR in dB and delta - linear
    # delta - is a Bragg component, which describes wind field impact
    # PR = 10*log10(SigmaHHwnr) - 10*log10(SigmaVVwnr)
#    # PR must not be > 0, be careful with masking
    #P[P>0] = 0
    PR = SigmaHHwnr/SigmaVVwnr
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

    # Find the relation of PR in/out-side the Slicks
    oilCropNorthPR = imcrop(ptsCropNorth, PR)
    oilCropPN = oilCropNorthPR[markersN==0].mean()/oilCropNorthPR[markersS==1].mean()

    oilCropSouthP = imcrop(ptsCropSouth, PR)
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

    oilCropNorthVV = imcrop(ptsCropNorth, SigmaVVwnr)
    oilCropNorthDelta = imcrop(ptsCropNorth, delta)
    oilCropNorthSigmaWb = imcrop(ptsCropNorth, SigmaWb)





    # Interpolating lat/lon to image size for future pcolormesh
    from numpy import linspace, arange, round, floor, ceil
    from scipy.interpolate import interp2d
    def interpLL(l, pixel, line, gx, gy):
        """
        Interpolating lat/lon to image size for future pcolormesh
        """
        fl = interp2d(pixel[0,:], line[:,0], l, kind='linear')
        l_new = fl(gx, gy)
        return l_new

    ptsCrop = ptsCropNorth
    RasterXSize = round((-ptsCrop[0,0]+ptsCrop[-1,0]))
    RasterYSize = round((-ptsCrop[0,-1]+ptsCrop[-1,-1]))
    gx = linspace(ptsCrop[0,0], ptsCrop[-1,0], RasterXSize+1)
    gy = linspace(ptsCrop[0,-1], ptsCrop[-1,-1], RasterYSize+1)
    
    RasterXSize = round(pixel[0,-1])
    RasterYSize = round(line[-1,0])
    gx = linspace(0, pixel[0,-1], RasterXSize+1)
    gy = linspace(0, line[-1,0], RasterYSize+1)    
    
    lat_new = interpLL(lat, pixel, line, gx, gy)
    lon_new = interpLL(lon, pixel, line, gx, gy)

    # save to MAT file
    from scipy.io import savemat
    print "Exporting"
    pn = '/home/mag/Documents/MATLAB/Leon_Agulhas_RS2/'
    expFn = 'OilSlickTrns/OilSlickTrns.mat'
    if not path.isfile(pn + expFn):
        print "Exporting to:\n" , pn + expFn
        savemat(pn + expFn, mdict={ \
                'VV':oilCropNorthVV, \
                'PD':oilCropNorthDelta, \
                'SigmaWb':oilCropNorthSigmaWb, \
                'PR':oilCropNorthPR, \
                'lat':lat_new, 'lon':lon_new, \
                }, do_compression=False)

    expFn = 'CurrentFrontTrns/CurrentFrontTrns.mat'
    if not path.isfile(pn + expFn):
        print "Exporting to:\n" , pn + expFn
        savemat(pn + expFn, mdict={ \
                'PD':delta, \
                'PR':PR, \
                'VV':SigmaVVwnr, \
                'SigmaWb':SigmaWb, \
                'lat':lat_new, 'lon':lon_new, \
                }, do_compression=False)

    import h5py
    pn = '/home/mag/Documents/MATLAB/Leon_Agulhas_RS2/'
    expFn = 'CurrentFrontTrns/CurrentFrontTrns.mat'
    expFn = 'OilSlickTrns/OilSlickTrns.mat'
    f = h5py.File(pn + expFn)
    SigmaWb = f["SigmaWb"]
    SigmaWb = flipud(rot90(array(SigmaWb)))
    lat = f["lat"]
    lat = flipud(rot90(array(lat)))
    lon = f["lon"]
    lon = flipud(rot90(array(lon)))
    y = f["tWb"]['transx']
    y = array(y)
    x = f["tWb"]['transy']
    x = array(x)
    pts = zeros((2,2))
    pts[:,0] = x
    pts[:,1] = y
#    transec_init = f["tWb"]['transec_init']
#    transec_init = array(transec_init)

    import distancelib
    reload(distancelib)
    # Get pixel resolution in Col/Row direction
    PxResCol, PxResRow = distancelib.getPixelResolution(lat, lon, SigmaWb)
    # Get the length of the transect. Input coords in (row,col)
    trnsLength = distancelib.getDistancePx([pts[0,1], pts[0,0]], [pts[1,1], pts[1,0]], lat, lon, SigmaWb)

    import transect
    reload(transect)
    transect.plotTrns(SigmaWbtrans, pts, lat, lon, SigmaWb, pn)




    # Making transect
    plt.figure()
    plt.imshow(oilCropNorthDelta)
    plt.gray()
    plt.clim(0, 0.02)
    plt.colorbar()
    pts = transect.imagePts()
    plt.savefig('/home/mag/PDtrns.png', facecolor='w', edgecolor='w', \
                 dpi=100, bbox_inches="tight", pad_inches=1.75)

    # Transect North Oil Slick
    pts = array([[610, 871], [822, 1007]])

    # [1] returns only Mean Transect
    PDtrans = transect.transect(oilCropNorthDelta, pts, 30)[1]
    PRtrans = transect.transect(oilCropNorthPR, pts, 30)[1]
    VVtrans = transect.transect(oilCropNorthVV, pts, 30)[1]
    SigmaWbtrans = transect.transect(oilCropNorthSigmaWb, pts, 30)[1]

#    PDtrans = transect.transect(oilCropNorthDelta, pts, 100)[0]






    plt.figure()
    plt.imshow(CurrentWb)
    plt.gray()
    plt.clim(0, 0.03)
    plt.colorbar()
    pts = transect.imagePts()
    plt.savefig('/home/mag/PDtrns.png', facecolor='w', edgecolor='w', \
                 dpi=100, bbox_inches="tight", pad_inches=1.75)

    # Transect Current Front
    pts = array([[1207., 4318.], [2473., 3763.]])
    
    PDtrans = transect.transect(delta, pts, 30)[1]
    PRtrans = transect.transect(PR, pts, 30)[1]
    VVtrans = transect.transect(SigmaVVwnr, pts, 30)[1]
    
    SigmaWbtrans = transect.transect(SigmaWb, pts, 30, method='cubic')[1]



    # use median filter to smooth
    from scipy.ndimage import median_filter as mflt
    PRtransSmthd = mflt(PRtrans, size=21, mode='nearest')
    PDtransSmthd = mflt(PDtrans, size=21, mode='nearest')
    VVtransSmthd = mflt(VVtrans, size=21, mode='nearest')
    SigmaWbtransSmthd = mflt(SigmaWbtrans, size=21, mode='nearest')

    plt.figure()
#    plt.imshow(SigmaWb)
    plt.imshow(SigmaWb)
    plt.gray()
    plt.clim(0, 0.02)
    plt.colorbar()
    ax = plt.gca().get_xbound()
    ay = plt.gca().get_ybound()
    x = pts[:, 0]
    y = pts[:, 1]
    plt.plot(x, y, '-o')
    plt.text(x[0],y[0],"A",bbox=dict(facecolor='w', alpha=0.7),stretch='expanded',fontsize=21)
    plt.text(x[1],y[1],"B",bbox=dict(facecolor='w', alpha=0.7),stretch='expanded',fontsize=21)
    plt.gca().set_xbound(ax)
    plt.gca().set_ybound(ay)
    plt.draw()
#    window = gtk.Window()
#    screen = window.get_screen()
#    width = screen.get_width()
#    height = screen.get_height()
#    mng = plt.get_current_fig_manager()
#    mng.resize(width, height)
    plt.savefig('/home/mag/PDtrns.png', facecolor='w', edgecolor='w', \
                 dpi=100, bbox_inches="tight", pad_inches=1.75)

    plt.figure()
    plt.plot(PDtransSmthd)
    plt.grid()
#    window = gtk.Window()
#    screen = window.get_screen()
#    width = screen.get_width()
#    height = screen.get_height()
#    mng = plt.get_current_fig_manager()
#    mng.resize(width, height)
    plt.savefig('/home/mag/PDtrnsLine.png', facecolor='w', edgecolor='w', \
                 dpi=100, bbox_inches="tight", pad_inches=1.75)

    plt.figure()
    plt.plot(PRtransSmthd)
    plt.grid()
#    window = gtk.Window()
#    screen = window.get_screen()
#    width = screen.get_width()
#    height = screen.get_height()
#    mng = plt.get_current_fig_manager()
#    mng.resize(width, height)
    plt.savefig('/home/mag/PRtransLine.png', facecolor='w', edgecolor='w', \
                 dpi=100, bbox_inches="tight", pad_inches=1.75)

    plt.figure()
    plt.plot(VVtransSmthd)
    plt.grid()
#    window = gtk.Window()
#    screen = window.get_screen()
#    width = screen.get_width()
#    height = screen.get_height()
#    mng = plt.get_current_fig_manager()
#    mng.resize(width, height)
    plt.savefig('/home/mag/VVtransLine.png', facecolor='w', edgecolor='w', \
                 dpi=100, bbox_inches="tight", pad_inches=1.75)

    plt.figure()
    plt.plot(SigmaWbtransSmthd)
    plt.grid()
#    window = gtk.Window()
#    screen = window.get_screen()
#    width = screen.get_width()
#    height = screen.get_height()
#    mng = plt.get_current_fig_manager()
#    mng.resize(width, height)
    plt.savefig('/home/mag/SigmaWbtransLine.png', facecolor='w', edgecolor='w', \
                 dpi=100, bbox_inches="tight", pad_inches=1.75)

    plt.close('all')


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
    label1 = 'Sigma VV [dB]'
    label2 = 'Sigma HH [dB]'
    label11 = 'Sigma VV [linear units]'
    label22 = 'Sigma HH [linear units]'
    label3 = 'PR [dB]'
    label33 = 'PR [linear units]'
    label4 = 'PD [linear units]'
    label5 = 'Wb contribution [linear units]'
    label6 = 'SigmaVVConjHH [linear units]'
    label7 = 'SigmaVVConjHHdelta [linear units]'

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
#    plotRS2(PR, lat, lon, pixel, line, \
#            scale=scale, m=m, clm=(-2.7,0), label=label3)
#    makeMap(ll_lon, ur_lon, \
#            ll_lat, ur_lat, m, name=label3, contour='land')
#    plt.close('all')

    # Coords to plot Figure number
    x, y = m(3.15, 42.15)

    # Figure 1 a - VV
    plotRS2(SigmaVVwnr, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(0,0.05), label=label11)
    plt.text(x,y,"a",bbox=dict(facecolor='w', alpha=0.5),stretch='expanded',fontsize=27)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label11, contour='land')
    plt.close('all')
    a = label11.split(' ')[0] # split the name if it has spaces
    system("rm " + "/home/mag/" + a + ".tiff")
    system("mv " + "/home/mag/" + a + "_ETOPO1.tiff" + " " + "/home/mag/Figure1a.tiff")
    system("convert -compress lzw /home/mag/Figure1a.tiff /home/mag/Figure1a.tiff")

    # Figure 1 b - HH
    plotRS2(SigmaHHwnr, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(0,0.05), label=label22)
    plt.text(x,y,"b",bbox=dict(facecolor='w', alpha=0.5),stretch='expanded',fontsize=27)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label22, contour='land')
    plt.close('all')
    a = label22.split(' ')[0] # split the name if it has spaces
    system("rm " + "/home/mag/" + a + ".tiff")
    system("mv " + "/home/mag/" + a + "_ETOPO1.tiff" + " " + "/home/mag/Figure1b.tiff")
    system("convert -compress lzw /home/mag/Figure1b.tiff /home/mag/Figure1b.tiff")

    # Figure 1 c - polarization ratio (PR) (always<1)
    plotRS2(PR, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(0.4,1), label=label33)
    plt.text(x,y,"c",bbox=dict(facecolor='w', alpha=0.5),stretch='expanded',fontsize=27)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label3, contour='land')
    plt.close('all')
    a = label33.split(' ')[0] # split the name if it has spaces
    system("rm " + "/home/mag/" + a + ".tiff")
    system("mv " + "/home/mag/" + a + "_ETOPO1.tiff" + " " + "/home/mag/Figure1c.tiff")
    system("convert -compress lzw /home/mag/Figure1c.tiff /home/mag/Figure1c.tiff")

    # Figure 1 d - polarization difference (PD) (always>0)
    plotRS2(delta, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(0,0.02), label=label4)
    plt.text(x,y,"d",bbox=dict(facecolor='w', alpha=0.5),stretch='expanded',fontsize=27)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label4, contour='land')
    plt.close('all')
    a = label4.split(' ')[0] # split the name if it has spaces
    system("rm " + "/home/mag/" + a + ".tiff")
    system("mv " + "/home/mag/" + a + "_ETOPO1.tiff" + " " + "/home/mag/Figure1d.tiff")
    system("convert -compress lzw /home/mag/Figure1d.tiff /home/mag/Figure1d.tiff")


    # Figure 4 - wave breaking contribution
    plotRS2(SigmaWb, lat, lon, pixel, line, \
            scale=scale, m=m, clm=(0,0.03), label=label5)
    plt.text(x,y,"d",bbox=dict(facecolor='w', alpha=0.5),stretch='expanded',fontsize=27)
    makeMap(ll_lon, ur_lon, \
            ll_lat, ur_lat, m, name=label5, contour='land')
    plt.close('all')
    a = label5.split(' ')[0] # split the name if it has spaces
    system("rm " + "/home/mag/" + a + ".tiff")
    system("mv " + "/home/mag/" + a + "_ETOPO1.tiff" + " " + "/home/mag/Figure4.tiff")
    system("convert -compress lzw /home/mag/Figure4.tiff /home/mag/Figure4.tiff")

    import transect
    reload(transect)
    transect.plotTrns(SigmaWbtrans, pts, lat, lon, SigmaWb, pn)


    # Figure 3
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
    # Coords to plot Figure number
    xSl, ySl = m(3.473, 42.173)
    # Figure 3 a - VV
    plotRS2(SigmaVVwnr, lat, lon, pixel, line, \
            scale=1, m=m, clm=(0,0.05), label=label11, ptsCrop=ptsCropNorth)
    plt.text(xSl,ySl,"a",bbox=dict(facecolor='w', alpha=0.5),stretch='expanded',fontsize=27)
    plt.savefig('/home/mag/Figure3a.tiff', facecolor='w', edgecolor='w', \
                 dpi=300, bbox_inches="tight", pad_inches=1.75)
    plt.close('all')
    # Figure 3 b - PR
    plotRS2(PR, lat, lon, pixel, line, \
            scale=1, m=m, clm=(0.5,1), label=label33, ptsCrop=ptsCropNorth)
    plt.text(xSl,ySl,"b",bbox=dict(facecolor='w', alpha=0.5),stretch='expanded',fontsize=27)
    plt.savefig('/home/mag/Figure3b.tiff', facecolor='w', edgecolor='w', \
                 dpi=300, bbox_inches="tight", pad_inches=1.75)
    plt.close('all')
#    # Figure 3 c - PD  (always>0)
#    plotRS2(delta, lat, lon, pixel, line, \
#            scale=1, m=m, clm=(0,0.03), label=label4, ptsCrop=ptsCropNorth)
#    plt.text(xSl,ySl,"c",bbox=dict(facecolor='w', alpha=0.5),stretch='expanded',fontsize=27)
#    plt.savefig('/home/mag/Figure3c.tiff', facecolor='w', edgecolor='w', \
#                 dpi=300, bbox_inches="tight", pad_inches=1.75)
#    plt.close('all')
    system("convert -compress lzw /home/mag/Figure3a.tiff /home/mag/Figure3a.tiff")
    system("convert -compress lzw /home/mag/Figure3b.tiff /home/mag/Figure3b.tiff")
#    system("convert -compress lzw /home/mag/Figure3c.tiff /home/mag/Figure3c.tiff")
#    # Southern
#        # setup nsper basemap
#    # Lat/Lon coords of image corners
#    lat_new, lon_new = intCrpLL(lat, lon, pixel, line, scale=1, ptsCrop=ptsCropSouth)
#    ll_lat = lat_new.min()
#    ur_lat = lat_new.max()
#    ll_lon = lon_new.min()
#    ur_lon = lon_new.max()
#    cent_lat = lat_new.mean()
#    cent_lon = lon_new.mean()
#    m = Basemap(llcrnrlat=ll_lat, urcrnrlat=ur_lat,\
#                llcrnrlon=ll_lon, urcrnrlon=ur_lon, \
#                resolution='c', projection='nsper', \
#                satellite_height=798000, \
#                lat_0=cent_lat,lon_0=cent_lon)
#    # Coords to plot Figure number
#    xSl, ySl = m(3.491, 41.973)
#    # Figure 3 d - VV
#    plotRS2(SigmaVVwnr, lat, lon, pixel, line, \
#            scale=1, m=m, clm=(0,0.05), label=label11, ptsCrop=ptsCropSouth)
#    plt.text(xSl,ySl,"d",bbox=dict(facecolor='w', alpha=0.5),stretch='expanded',fontsize=27)
#    plt.savefig('/home/mag/Figure3d.tiff', facecolor='w', edgecolor='w', \
#                 dpi=300, bbox_inches="tight", pad_inches=1.75)
#    plt.close('all')
#    # Figure 3 e - PR (always<1)
#    plotRS2(PR, lat, lon, pixel, line, \
#            scale=1, m=m, clm=(0.4,1), label=label33, ptsCrop=ptsCropSouth)
#    plt.text(xSl,ySl,"e",bbox=dict(facecolor='w', alpha=0.5),stretch='expanded',fontsize=27)
#    plt.savefig('/home/mag/Figure3e.tiff', facecolor='w', edgecolor='w', \
#                 dpi=300, bbox_inches="tight", pad_inches=1.75)
#    plt.close('all')
#    # Figure 3 f - PD  (always>0)
#    plotRS2(delta, lat, lon, pixel, line, \
#            scale=1, m=m, clm=(0,0.03), label=label4, ptsCrop=ptsCropSouth)
#    plt.text(xSl,ySl,"f",bbox=dict(facecolor='w', alpha=0.5),stretch='expanded',fontsize=27)
#    plt.savefig('/home/mag/Figure3f.tiff', facecolor='w', edgecolor='w', \
#                 dpi=300, bbox_inches="tight", pad_inches=1.75)
#    plt.close('all')
#    system("convert -compress lzw /home/mag/Figure3d.tiff /home/mag/Figure3d.tiff")
#    system("convert -compress lzw /home/mag/Figure3e.tiff /home/mag/Figure3e.tiff")
#    system("convert -compress lzw /home/mag/Figure3f.tiff /home/mag/Figure3f.tiff")
    import transect
    reload(transect)
    transect.plotTrns(PRtrans, pts, lat, lon, PR, pixel, line, pn='/home/mag/', \
                      label='PR [linear units]', clm=(0.5,1), \
                      ptsCrop=ptsCropNorth, m=m, scale=1)

#    plotRS2(SigmaVVConjHHWnr, lat, lon, pixel, line, \
#            scale=scale, m=m, clm=(0,0.05), label=label6)
#    makeMap(ll_lon, ur_lon, \
#            ll_lat, ur_lat, m, name=label6, contour='land')
#    plt.close('all')
#
#    plotRS2(SigmaVVConjHHdelta, lat, lon, pixel, line, \
#            scale=scale, m=m, clm=(0,5), label=label7)
#    makeMap(ll_lon, ur_lon, \
#            ll_lat, ur_lat, m, name=label7, contour='land')
#    plt.close('all')
#
#    # Plot SST on a new figure
#    plt.figure()
#    g1sst(ll_lon, ur_lon, \
#            ll_lat, ur_lat, m, name='g1sst', contour=None)
#    makeMap(ll_lon, ur_lon, \
#            ll_lat, ur_lat, m, name='g1sst', contour='land')
#    plt.close('all')



if __name__ == "__main__":
    main(sys.argv[1:])