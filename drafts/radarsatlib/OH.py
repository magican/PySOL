#! /usr/bin/env python
# -*- coding: utf-8 -*-"""
"""
Created on Sun Oct 11 12:28:01 2010

@author: S K Tomer
"""
from numpy import *
from itertools import count, izip 
from gdalconst import *
import os
import gdal


qDB=linspace(-25,-2,100)
qLI=10**(qDB/10);
gamma0=linspace(2,30,100)
theta=0.37
QDB,Gamma0 = meshgrid(qDB,gamma0)

pLI=(1-((2.0*theta/pi)**(1.0/(3.0*Gamma0)))*(1-QDB/(0.23*sqrt(Gamma0))))**2
pDB=10*log10(pLI)

def LUT(q,p):  # The OH model
    diffq=(qDB-q)**2
    minivalue, minindexX=min(izip(diffq,count()))
    diffp=(pDB[minindexX][:]-p)**2
    minivalue, minindexY=min(izip(diffp,count()))
    return Gamma0[minindexX][minindexY]

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
    
    # calculate the p and q
    p=HH-VV
    q=HV-VV
    
    # convert the p and q from dB to linear scale
    p=10**(p/10);
    q=10**(q/10);
    
    # calculate the soil moisture (mv)
    mv=zeros(IA.shape)
    for i in range(IA.shape[0]):
        for j in range(IA.shape[1]):
            q0=LUT(q[i][j],p[i][j])
            ep=((1+math.sqrt(q0))/(1-math.sqrt(q0)))**2
            mv[i][j]=(-530+292*ep-5.5*ep**2+0.043*ep**3)*1e-4
        
    # save the soil moisture (mv) as Geotiff
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create("/home/tomer/RADARSAT/"+file+"/OH",HH.shape[1],HH.shape[0],1,gdal.GDT_Float32)
    output_dataset.SetGeoTransform(geotransform)
    output_dataset.SetGCPs(gcps, gcpproj)
    output_dataset.GetRasterBand(1).WriteArray(mv, 0, 0)
    output_dataset = None
    
    # print the prcocessing
    print file+ " is done"   