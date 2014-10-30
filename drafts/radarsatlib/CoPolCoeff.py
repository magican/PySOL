#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 17:06:59 2010
by taking help from
http://benjamindeschamps.ca/blog/2009/11/12/processing-radarsat-2-imagery-reading-raw-data-and-saving-rgb-composites/
@author: sat kumar tomer
"""
# import the required library
import numpy as np
import math
from osgeo import gdal
import os
from xml.dom import minidom
import scipy.signal

dirName=os.listdir("/home/tomer/RADARSAT")
dirName.remove('scripts')
dirName.remove('documents')
for file in dirName:
    # define the paths
    InPath="/home/tomer/RADARSAT/"+file+"/rawdata/"
    OutPath="/home/tomer/RADARSAT/"+file+"/"
    # read the data, GCPs and projection
    dataset = gdal.Open("RADARSAT_2_CALIB:SIGMA0:" + InPath + "product.xml")
    geotransform = dataset.GetGeoTransform()
    gcps = dataset.GetGCPs()
    gcpproj = dataset.GetGCPProjection()

    # calculate the sinclair matrix
    S_HH = dataset.GetRasterBand(1).ReadAsArray()
    S_HV = dataset.GetRasterBand(3).ReadAsArray()
    S_VH = dataset.GetRasterBand(4).ReadAsArray()
    S_VV = dataset.GetRasterBand(2).ReadAsArray()

    # calculate the co-polarized correlation coefficient
    nume = S_HH*np.conj(S_HV)
    deno = np.sqrt((np.absolute(S_HH)*np.absolute(S_VV))**2)
    rho = np.absolute(nume/deno)
        
    # filter the image using median filter of window (7X7)
    rho=scipy.signal.medfilt2d(rho, kernel_size=7)
    
    # filter the image using wiener filter of window (7X7)
    #rho=scipy.signal.wiener(rho,mysize=(7,7),noise=None)
    rho[rho>1]=1
    rho[rho<0]=0
    # save the data as Geotiff
    if os.path.exists(OutPath + "CoPolCoeff"):
        os.remove(OutPath + "CoPolCoeff")

    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(OutPath + "CoPolCoeff",dataset.RasterXSize,dataset.RasterYSize,1,gdal.GDT_Float32)
    output_dataset.SetGeoTransform(geotransform)
    output_dataset.SetGCPs(gcps, gcpproj)
    output_dataset.GetRasterBand(1).WriteArray(rho, 0, 0)
    output_dataset = None

    # print the prcocessing
    print file+ " is done"    