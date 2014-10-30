#!/usr/bin/env python
# coding: utf-8
"""
Created on Fri Mar 23 17:05:57 2012
using
http://benjamindeschamps.ca/blog/2009/11/12/\
processing-radarsat-2-imagery-reading-raw-data-and-saving-rgb-composites/
@author: mag
"""
# import the required library
#import numpy
from numpy import log10, absolute, zeros, reshape
from math import asin
from osgeo import gdal
from xml.dom import minidom
from scipy.signal import medfilt2d, wiener
from matplotlib.mlab import find

from multiprocessing import Process, Queue

# Filter warnings
import warnings
warnings.filterwarnings(action = "ignore", category=RuntimeWarning)
warnings.filterwarnings(action = "ignore", category=FutureWarning)

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
    que.put(S_HH)

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
    que.put(S_VV)

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

def medfilt(data, kernel_size=7, que=[]):
    dataMed = medfilt2d(data, kernel_size)
    que.put(dataMed)

def lee(data, mysize=(7,7), noise=None, que=[]):
    dataWnr = wiener(data, mysize, noise)
    que.put(dataWnr)

class calibRS2par:
    """\
    Initial class for calibrating RS2 images and getting the incidence angle
    Usage:
        calib = calibRS2(pn, xOff=0, yOff=0, \
        xS=None, yS=None, xBufScale=None, yBufScale=None, s='abs')
            pn         = pathname
            xOff, yOff = offsets to crop the image
            xS, yS     = sizes to crop
            xBufScale, yBufScale = scale the image
            s='abs' || s = 'sigma0' - sigma0 "linear" and in "dB"
        All the image is read by default
    """
    def __init__(self, inpath, xOff=0, yOff=0, xS=None, yS=None, \
    xBufScale=None, yBufScale=None, s='abs'):
        """Calibrating the RS2 image"""

        self.inpath = inpath

        # read the data, GCPs and projection
        dataset = gdal.Open("RADARSAT_2_CALIB:SIGMA0:" + inpath + "product.xml")
        self.geotransform = dataset.GetGeoTransform()
        self.gcps = dataset.GetGCPs()
        self.gcpproj = dataset.GetGCPProjection()
        self.RasterXSize = dataset.RasterXSize
        self.RasterYSize = dataset.RasterYSize

        # define the image crop and quicklook
        if xBufScale is not None:
            xBufScale = self.RasterXSize/xBufScale
        if yBufScale is not None:
            yBufScale = self.RasterYSize/yBufScale

        # Get the polarizations
        xmldoc = minidom.parse(inpath + "product.xml")
        polarizations = xmldoc.getElementsByTagName('polarizations')
        polarizations = polarizations[0].toxml()
        self.polarizations = polarizations[15:-16]
        print "%s" % self.polarizations

        if dataset.RasterCount == 4:
    #        core_ct = sysconf('SC_NPROCESSORS_ONLN')
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
            self.SigmaHH = queue1.get()
            self.S_HH = queue1.get()
            self.SigmaVV = queue2.get()
            self.S_VV = queue2.get()
            self.SigmaHV = queue3.get()
            self.SigmaVH = queue4.get()
            hilo1.join()
            hilo2.join()
            hilo3.join()
            hilo4.join()
#            return SigmaHH, SigmaVV, SigmaHV, SigmaVH

        elif dataset.RasterCount == 2:
            #~ #Prepare to call a function from string name
            #~ possibles = globals().copy()
            #~ possibles.update(locals())
            #~ method1 = possibles.get('sigma' + self.polarizations[0:2])
            #~ method2 = possibles.get('sigma' + self.polarizations[3:5])
            #~ queue1 = Queue()
            #~ queue2 = Queue()
            #~ hilo1 = Process(target=method1, \
                    #~ args=(dataset, xOff, yOff, xS, yS, \
                    #~ xBufScale, yBufScale, s, queue1))
            #~ hilo2 = Process(target=method2, \
                    #~ args=(dataset, xOff, yOff, xS, yS, \
                    #~ xBufScale, yBufScale, s, queue2))
            #~ print "Launching threads"
            #~ hilo1.start()
            #~ hilo2.start()
            #~ #and we're assigning the return value to sigma
            #~ #creating vars from string polarizations
            #~ setattr(self, 'sigma' + self.polarizations[0:2], queue1.get())
            #~ setattr(self, 'sigma' + self.polarizations[3:5], queue2.get())
            #~ hilo1.join()
            #~ hilo2.join()

            S1 = dataset.GetRasterBand(1).ReadAsArray()
            S2 = dataset.GetRasterBand(2).ReadAsArray()
            # calculate the magnitude
            S_ABS1 = absolute(S1)
            S_ABS2 = absolute(S2)
            if s == 'abs':
				setattr(self, 'sigma' + self.polarizations[0:2], S_ABS1**2)
				setattr(self, 'sigma' + self.polarizations[3:5], S_ABS2**2)
            # calculate the sigma_naught
            if s == 'sigma0':
				setattr(self, 'sigma' + self.polarizations[0:2], 2*10*log10(S_ABS1))
				setattr(self, 'sigma' + self.polarizations[3:5], 2*10*log10(S_ABS2))

        elif dataset.RasterCount == 1:
            S = dataset.GetRasterBand(1).ReadAsArray()

            # calculate the magnitude
            S_ABS = absolute(S)

            if s == 'abs':
                self.Sigma0 = S_ABS**2

            # calculate the sigma_naught
            if s == 'sigma0':
                self.Sigma0 = 2*10*log10(S_ABS)

    def incidence_angle(self):
        """Getting the incidence angle"""
        # calculate the incidence angle image
        xmldoc = minidom.parse(self.inpath+"lutSigma.xml")
        SigmaGains = xmldoc.getElementsByTagName('gains')
        SigmaGains = SigmaGains[0].toxml()
        SigmaGains = SigmaGains[7:-8]
        SigmaGains = SigmaGains.split(' ')
        xmldoc = minidom.parse(self.inpath+"lutBeta.xml")
        BetaGains = xmldoc.getElementsByTagName('gains')
        BetaGains = BetaGains[0].toxml()
        BetaGains = BetaGains[7:-8]
        BetaGains = BetaGains.split(' ')
        IncidenceAngle = zeros( (1, self.RasterXSize) )
        for j in range(IncidenceAngle.shape[1]):
            IncidenceAngle[0, j] = \
            asin((float(BetaGains[j])/float(SigmaGains[j]))**2)
        return IncidenceAngle

    def xml2geo(self):
        """
        Getting the lattitudes and longitudes from the product.xml
        As well as corresponding row and coloumn
        """
        # Get the lats/lons, rows/cols
        xmldoc = minidom.parse(self.inpath + "product.xml")
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
        return lat, lon, line, pixel

    def speckle_filter(self, filter_name='wiener', ws=7):
        """
        Filter the image using 'median' filtering methods
        """
        if filter_name == 'median':
            if len(self.polarizations) == 11:
                # filter the image using median filter
                queue1 = Queue()
                queue2 = Queue()
                queue3 = Queue()
                queue4 = Queue()
                hilo1 = Process(target=medfilt, \
                        args=(self.SigmaHH, ws, queue1))
                hilo2 = Process(target=medfilt, \
                        args=(self.SigmaVV, ws, queue2))
                hilo3 = Process(target=medfilt, \
                        args=(self.SigmaHV, ws, queue3))
                hilo4 = Process(target=medfilt, \
                        args=(self.SigmaVH, ws, queue4))
                print "Launching threads"
                hilo1.start()
                hilo2.start()
                hilo3.start()
                hilo4.start()
                #and we're assigning the return value to sigma
                self.SigmaHHmed = queue1.get()
                self.SigmaVVmed = queue2.get()
                self.SigmaHVmed = queue3.get()
                self.SigmaVHmed = queue4.get()
                hilo1.join()
                hilo2.join()
                hilo3.join()
                hilo4.join()
    #            return self.SigmaHHmf, self.SigmaHVmf, \
    #                self.SigmaVHmf, self.SigmaVVmf
            elif len(self.polarizations) == 5:
                queue1 = Queue()
                queue2 = Queue()
                hilo1 = Process(target=medfilt, \
                        args=(eval('self.sigma' + self.polarizations[0:2]), ws, queue1))
                hilo2 = Process(target=medfilt, \
                        args=(eval('self.sigma' + self.polarizations[3:5]), ws, queue2))
                print "Launching threads"
                hilo1.start()
                hilo2.start()
                #~ creating vars from string polarizations
                setattr(self, 'sigma' + self.polarizations[0:2] + 'med', queue1.get())
                setattr(self, 'sigma' + self.polarizations[3:5] + 'med', queue2.get())
                hilo1.join()
                hilo2.join()
        elif filter_name == 'wiener':
            if len(self.polarizations) == 11:
                # filter the image using wiener filter
                queue1 = Queue()
                queue2 = Queue()
                queue3 = Queue()
                queue4 = Queue()
                hilo1 = Process(target=lee, \
                        args=(self.SigmaHH, (ws,ws), None, queue1))
                hilo2 = Process(target=lee, \
                        args=(self.SigmaVV, (ws,ws), None, queue2))
                hilo3 = Process(target=lee, \
                        args=(self.SigmaHV, (ws,ws), None, queue3))
                hilo4 = Process(target=lee, \
                        args=(self.SigmaVH, (ws,ws), None, queue4))
                print "Launching threads"
                hilo1.start()
                hilo2.start()
                hilo3.start()
                hilo4.start()
                #and we're assigning the return value to sigma
                self.SigmaHHwnr = queue1.get()
                self.SigmaVVwnr = queue2.get()
                self.SigmaHVwnr = queue3.get()
                self.SigmaVHwnr = queue4.get()
                hilo1.join()
                hilo2.join()
                hilo3.join()
                hilo4.join()
            elif len(self.polarizations) == 5:
                queue1 = Queue()
                queue2 = Queue()
                hilo1 = Process(target=lee, \
                        args=(eval('self.sigma' + self.polarizations[0:2]), (ws,ws), None, queue1))
                hilo2 = Process(target=lee, \
                        args=(eval('self.sigma' + self.polarizations[3:5]), (ws,ws), None, queue2))
                print "Launching threads"
                hilo1.start()
                hilo2.start()
                #~ creating vars from string polarizations
                setattr(self, 'sigma' + self.polarizations[0:2] + 'wnr', queue1.get())
                setattr(self, 'sigma' + self.polarizations[3:5] + 'wnr', queue2.get())
                hilo1.join()
                hilo2.join()
        else:
            print 'Please specify the name of the filter: "median" or "wiener"'

    def save_tiff(self):
        # save the data as Geotiff
        driver = gdal.GetDriverByName('GTiff')
        output_dataset = driver.Create(self.inpath + 'sigma.tiff', \
            self.RasterXSize, self.RasterYSize, 4, gdal.GDT_Float32)
        output_dataset.SetGeoTransform(self.geotransform)
        output_dataset.SetGCPs(self.gcps, self.gcpproj)
        output_dataset.GetRasterBand(1).WriteArray(self.SigmaHH, 0, 0)
        output_dataset.GetRasterBand(2).WriteArray(self.SigmaVV, 0, 0)
        output_dataset.GetRasterBand(3).WriteArray(self.SigmaHV, 0, 0)
        output_dataset.GetRasterBand(4).WriteArray(self.SigmaVH, 0, 0)
        output_dataset = None
