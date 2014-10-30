#! /usr/bin/env python
# -*- coding: utf-8 -*-"""
"""
Created on Sun Oct 10 00:42:01 2010

@author: tomer
"""
from scitools.all import *
from gdalconst import *
import os
import gdal

# perform the analysis on all the images
dirName=os.listdir("/home/tomer/RADARSAT")
dirName.remove('scripts')
dirName.remove('documents')

for file in dirName:
    # open the files
    fidSigma = gdal.Open("/home/tomer/RADARSAT/"+file+"/filtered",GA_ReadOnly)
    fidIA=gdal.Open("/home/tomer/RADARSAT/"+file+"/IncidenceAngle",GA_ReadOnly)
    if fidSigma is None:
        print 'file does not exist'
    # set the progress index to zero
    ProgressIndex=0
    # read the projection details, because the projection details of the backscattering image
    #filtered and IncidenceAngle are same, they are read only from one file
    geotransform = fidSigma.GetGeoTransform()
    gcps = fidSigma.GetGCPs()
    gcpproj = fidSigma.GetGCPProjection()
    
    # read the HH, HV, VH, VV 
    HH = fidSigma.GetRasterBand(1).ReadAsArray()
    HV = fidSigma.GetRasterBand(2).ReadAsArray()
    VH = fidSigma.GetRasterBand(3).ReadAsArray()
    VV = fidSigma.GetRasterBand(4).ReadAsArray()
    IA=fidIA.GetRasterBand(1).ReadAsArray()
    
    # we are done with reading files, hence close them
    fidSigma = None
    fidIA = None

pcolor(HH)