#! /usr/bin/env python
# -*- coding: utf-8 -*-"""
"""
Created on Sun Oct 10 00:42:01 2010

@author: tomer
"""
# import required libraries
import gdal
from gdalconst import *
from pca_module import *
from numpy import *

# open the files
fidSigma = gdal.Open("/home/tomer/RADARSAT/15May10/filtered",GA_ReadOnly)
fidIA=gdal.Open("/home/tomer/RADARSAT/15May10/IncidenceAngle",GA_ReadOnly)
if fidSigma is None:
    print 'file does not exist'

# read the projection details, because the projection details of the backscattering image
#filtered and IncidenceAngle are same, they are read only from one file
geotransformSigma = fidSigma.GetGeoTransform()
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

# scale the data (mean=0) and (std=1)
HH=HH-numpy.mean(HH)
HH=HH/numpy.std(HH)
HV=HV-numpy.mean(HV)
HV=HV/numpy.std(HV)
VH=VH-numpy.mean(VH)
VH=VH/numpy.std(VH)
VV=VV-numpy.mean(VV)
VV=VV/numpy.std(VV)

HH=HH.flatten()
HV=HV.flatten()
VH=VH.flatten()
VV=VV.flatten()

data=[HH,HV,VH,VV]
data=numpy.array(data)

T, P, explained_var = PCA_nipals(data, standardize=True,E_matrices=True)