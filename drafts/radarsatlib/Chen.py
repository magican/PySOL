#! /usr/bin/env python
# -*- coding: utf-8 -*-"""
"""
Created on Sun Oct 11 12:28:01 2010
Based on the article:
A Simple Model for Retrieving Bare Soil Moisture from Radar-Scattering Coefficients
By K. S. Chen, S. K. Yen, and W. P. Huang

@author: S K Tomer
"""
# import required libraries
import gdal
from gdalconst import *
import numpy
import math
from numpy import exp

# perform the analysis on all the images
import os
dirName = os.listdir("/home/tomer/RADARSAT")
dirName.remove('scripts')
dirName.remove('documents')

# define the coefficients
C1 = -0.09544
C2 = -0.00971
C3 = 0.029238
C4 = -1.74678

# define the frequency
f = 5.4

for file in dirName:
    # open the files
    fidSigma = gdal.Open("/home/tomer/RADARSAT/"+file+"/filtered",GA_ReadOnly)
    fidIA = gdal.Open("/home/tomer/RADARSAT/"+file+"/IncidenceAngle",GA_ReadOnly)
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
    
    # calculate the p and q
    p=HH-VV
        
    # calculate the soil moisture (mv)
    mv=exp(C1*p+C2*IA+C3*f+C4)
               
    # save the soil moisture (mv) as Geotiff
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create("/home/tomer/RADARSAT/"+file+"/Chen",HH.shape[1],HH.shape[0],1,gdal.GDT_Float32)
    output_dataset.SetGeoTransform(geotransform)
    output_dataset.SetGCPs(gcps, gcpproj)
    output_dataset.GetRasterBand(1).WriteArray(mv, 0, 0)
    output_dataset = None
    
    # print the prcocessing
    print numpy.mean(IA)
    print file+ " is done"
