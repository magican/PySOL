#! /usr/bin/env python
# -*- coding: utf-8 -*-"""
"""
Created on Sun Oct 11 12:28:01 2010

@author: S K Tomer
"""
# import required libraries
import gdal
from gdalconst import *
from pca_module import *
from numpy import *

# perform the analysis on all the images
import os
dirName=os.listdir("/home/tomer/RADARSAT")
dirName.remove('scripts')
dirName.remove('documents')

for file in dirName:
    # open the files
    fidSigma = gdal.Open("/home/tomer/RADARSAT/"+file+"/filtered",GA_ReadOnly)
    fidIA=gdal.Open("/home/tomer/RADARSAT/"+file+"/IncidenceAngle",GA_ReadOnly)
    if fidSigma is None:
        print 'file does not exist'
    
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
    
    # coverting the Sigma from dB scale to linear scale
    #HH=10.0**(HH/10.0)
    #HV=10.0**(HV/10.0)
    #VH=10.0**(VH/10.0)
    #VV=10.0**(VV/10.0)
    
    # calculate the p and q
    #p=HH/VV
    #q=HV/VV
    
    # calculate the ep
    lamda=5.6
    ep=(26.5+14*VV-11*HH-255*log10(cos(IA))-130*log10(sin(IA))-21*log10(lamda))/(3.36*tan(IA))
    #ep=numpy.log10(((HH**(0.7857))/VV)*10**(-0.19)*(numpy.cos(IA)**1.82)*(numpy.sin(IA)**0.93)*lamda**0.15)/(-0.024*numpy.tan(IA));
        
    # calculate the soil moisture 
    mv=(-530+292*ep-5.5*ep**2+0.043*ep**3)*1e-4
    
    # save the soil moisture (mv) as Geotiff
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create("/home/tomer/RADARSAT/"+file+"/Dubois",HH.shape[1],HH.shape[0],1,gdal.GDT_Float32)
    output_dataset.SetGeoTransform(geotransform)
    output_dataset.SetGCPs(gcps, gcpproj)
    output_dataset.GetRasterBand(1).WriteArray(mv, 0, 0)
    output_dataset = None
    
    # print the prcocessing
    print file+ " is done"