#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 19:01:08 2012

@author: mag
"""

from multiprocessing import Process, Queue
from time import time

from osgeo import gdal
from xml.dom import minidom

from numpy import log10, absolute, zeros, reshape
from scipy.signal import medfilt2d, wiener
from math import asin

#add an argument to function for assigning a queue
def sigmaHH(dataset, xOff=0, yOff=0, xS=None, yS=None, \
    xBufScale=None, yBufScale=None, s='abs', que=[]):
    # calculate the sinclair matrix
    S_HH = dataset.GetRasterBand(1).ReadAsArray(xoff=xOff, yoff=yOff, \
    win_xsize=xS, win_ysize=yS, \
    buf_xsize=xBufScale, buf_ysize=yBufScale)
    # calculate the magnitude
    S_HH_ABS = absolute(S_HH)
    # calculate the linear sigma_naught
    if s == 'abs':
        SigmaHH = S_HH_ABS**2
    # calculate the sigma_naught in dB
    if s == 'sigma0': 
        SigmaHH = 2*10*log10(S_HH_ABS)
    #we're putting return value into queue
    que.put(SigmaHH)

def sigmaVV(dataset, xOff=0, yOff=0, xS=None, yS=None, \
    xBufScale=None, yBufScale=None, s='abs', que=[]):
    # calculate the sinclair matrix
    S_VV = dataset.GetRasterBand(2).ReadAsArray(xoff=xOff, yoff=yOff, \
    win_xsize=xS, win_ysize=yS, \
    buf_xsize=xBufScale, buf_ysize=yBufScale)
    # calculate the magnitude
    S_VV_ABS = absolute(S_VV)
    # calculate the linear sigma_naught
    if s == 'abs':
        SigmaVV = S_VV_ABS**2
    # calculate the sigma_naught in dB
    if s == 'sigma0': 
        SigmaVV = 2*10*log10(S_VV_ABS)
#    SigmaVVwnr = wiener(SigmaVV,mysize=(7,7),noise=None)
    #we're putting return value into queue
    que.put(SigmaVV)

def sigmaHV(dataset, xOff=0, yOff=0, xS=None, yS=None, \
    xBufScale=None, yBufScale=None, s='abs', que=[]):
    # calculate the sinclair matrix
    S_HV = dataset.GetRasterBand(3).ReadAsArray(xoff=xOff, yoff=yOff, \
    win_xsize=xS, win_ysize=yS, \
    buf_xsize=xBufScale, buf_ysize=yBufScale)
    # calculate the magnitude
    S_HV_ABS = absolute(S_HV)
    # calculate the linear sigma_naught
    if s == 'abs':
        SigmaHV = S_HV_ABS**2
    # calculate the sigma_naught in dB
    if s == 'sigma0': 
        SigmaHV = 2*10*log10(S_HV_ABS)
    que.put(SigmaHV)

def sigmaVH(dataset, xOff=0, yOff=0, xS=None, yS=None, \
    xBufScale=None, yBufScale=None, s='abs', que=[]):
    # calculate the sinclair matrix
    S_VH = dataset.GetRasterBand(4).ReadAsArray(xoff=xOff, yoff=yOff, \
    win_xsize=xS, win_ysize=yS, \
    buf_xsize=xBufScale, buf_ysize=yBufScale)
    # calculate the magnitude
    S_VH_ABS = absolute(S_VH)
    # calculate the linear sigma_naught
    if s == 'abs':
        SigmaVH = S_VH_ABS**2
    # calculate the sigma_naught in dB
    if s == 'sigma0': 
        SigmaVH = 2*10*log10(S_VH_ABS)
    que.put(SigmaVH)

def incidence_angle(inpath):
    """Getting the incidence angle"""
    RasterXSize = \
    gdal.Open("RADARSAT_2_CALIB:SIGMA0:" + inpath + "product.xml").RasterXSize
    # calculate the incidence angle image
    xmldoc = minidom.parse(inpath+"lutSigma.xml")
    SigmaGains = xmldoc.getElementsByTagName('gains')
    SigmaGains = SigmaGains[0].toxml()
    SigmaGains = SigmaGains[7:-8]
    SigmaGains = SigmaGains.split(' ')
    xmldoc = minidom.parse(inpath+"lutBeta.xml")
    BetaGains = xmldoc.getElementsByTagName('gains')
    BetaGains = BetaGains[0].toxml()
    BetaGains = BetaGains[7:-8]
    BetaGains = BetaGains.split(' ')
    IncidenceAngle = zeros( (1, RasterXSize) )
    for j in range(IncidenceAngle.shape[1]):
        IncidenceAngle[0, j] = \
        asin((float(BetaGains[j])/float(SigmaGains[j]))**2)
    return IncidenceAngle

def calib(inpath, xOff=0, yOff=0, xS=None, yS=None, \
xBufScale=None, yBufScale=None, s='abs'):
    """Calibrating the RS2 image"""

    # read the data, GCPs and projection
    dataset = gdal.Open("RADARSAT_2_CALIB:SIGMA0:" + inpath + "product.xml")
    geotransform = dataset.GetGeoTransform()
    gcps = dataset.GetGCPs()
    gcpproj = dataset.GetGCPProjection()
    RasterXSize = dataset.RasterXSize
    RasterYSize = dataset.RasterYSize

    # define the image crop and quicklook
    if xBufScale is not None:
        xBufScale = RasterXSize/xBufScale
    if yBufScale is not None: 
        yBufScale = RasterYSize/yBufScale
   
    # Get the polarizations        
    xmldoc = minidom.parse(inpath + "product.xml")
    polarizations = xmldoc.getElementsByTagName('polarizations')
    polarizations = polarizations[0].toxml()
    polarizations = polarizations[15:-16]
    print "%s" % polarizations

    if dataset.RasterCount == 4:

#        core_ct = sysconf('SC_NPROCESSORS_ONLN')
        print "Starting main program"
#        print "Running %s cores " % core_ct
        #create a queue object
        queue1 = Queue()
        queue2 = Queue()
        queue3 = Queue()
        queue4 = Queue()
        hilo1 = Process(target=sigmaHH, \
                args=(dataset, xOff, yOff, xS, yS, \
                xBufScale, yBufScale, s, queue1))
        hilo2 = Process(target=sigmaVV, \
                args=(dataset, xOff, yOff, xS, yS, \
                xBufScale, yBufScale, s, queue2))
        hilo3 = Process(target=sigmaHV, \
                args=(dataset, xOff, yOff, xS, yS, \
                xBufScale, yBufScale, s, queue3))
        hilo4 = Process(target=sigmaVH, \
                args=(dataset, xOff, yOff, xS, yS, \
                xBufScale, yBufScale, s, queue4))
        print "Launching threads"
        hilo1.start()
        hilo2.start()
        hilo3.start()
        hilo4.start()
        #and we're assigning the return value to sigma        
        SigmaHH = queue1.get()
        SigmaVV = queue2.get()
        SigmaHV = queue3.get()
        SigmaVH = queue4.get()
        hilo1.join()
        hilo2.join()
        hilo3.join()
        hilo4.join()
        
        return SigmaHH, SigmaVV, SigmaHV, SigmaVH
 
    elif dataset.RasterCount == 1:
        S = dataset.GetRasterBand(1).ReadAsArray()

        # calculate the magnitude
        if s == 'abs':
            S_ABS = absolute(S)

        # calculate the sigma_naught
        if s == 'sigma0':
            Sigma0 = 2*10*log10(S_ABS)
    
if __name__ == '__main__':
    calib(inpath=None, xOff=0, yOff=0, xS=None, yS=None, \
    xBufScale=None, yBufScale=None, s='abs')