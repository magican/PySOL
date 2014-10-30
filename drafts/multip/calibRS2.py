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
from scipy.signal import medfilt2d
from matplotlib.mlab import find

# Filter warnings
import warnings
warnings.filterwarnings(action = "ignore", category=RuntimeWarning)
warnings.filterwarnings(action = "ignore", category=FutureWarning)
        
class calibRS2:
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
        
    def speckle_filter(self, filter_name, ws):
        """
        Filter the image using 'median' filtering methods
        """
        if filter_name == 'median':
            # filter the image using median filter 
            self.SigmaHHmf = medfilt2d(self.SigmaHH, kernel_size=ws)
            self.SigmaHVmf = medfilt2d(self.SigmaHV, kernel_size=ws)
            self.SigmaVHmf = medfilt2d(self.SigmaVH, kernel_size=ws)
            self.SigmaVVmf = medfilt2d(self.SigmaVV, kernel_size=ws)
            return self.SigmaHHmf, self.SigmaHVmf, \
                self.SigmaVHmf, self.SigmaVVmf
            
        else:
            print 'Please specify the name of the filter: "median"'

    def save_tiff(self):
        # save the data as Geotiff
        driver = gdal.GetDriverByName('GTiff')
        output_dataset = driver.Create(self.inpath + 'sigma.tiff', self.RasterXSize, self.RasterYSize,4,gdal.GDT_Float32)
        output_dataset.SetGeoTransform(self.geotransform)
        output_dataset.SetGCPs(self.gcps, self.gcpproj)
        output_dataset.GetRasterBand(1).WriteArray(self.SigmaHH, 0, 0)
        output_dataset.GetRasterBand(2).WriteArray(self.SigmaHV, 0, 0)
        output_dataset.GetRasterBand(3).WriteArray(self.SigmaVH, 0, 0)
        output_dataset.GetRasterBand(4).WriteArray(self.SigmaVV, 0, 0)
        output_dataset = None