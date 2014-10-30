#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 16 13:30:19 2010

@author: s k tomer
"""

import xlrd
from numpy import *
from pylab import *
from osgeo import gdal
from osgeo.gdalconst import * 
import os
from scipy import polyfit
from matplotlib.pyplot import *

# extract data from images
def DataExtract(IL,IP,data):
    DataE=zeros(len(IL))
    for i in range(len(IP)):
        DataE[i]=data[int(IL[i])][int(IP[i])]
        Datae=zeros(50)
    for i in range(50):
        Datae[i]=median(DataE[int(ind_pre[i]):int(ind_pre[i+1])])
    return Datae

# regression coefficient for x and y
def regressNan(y,x,xhat):
    A = vstack([ones(len(x)),x]).T
    coeff = linalg.lstsq(A, y)[0]
    yhat = coeff[0]+coeff[1]*xhat
    return yhat

def regress2(y,x1,x2):
    A = vstack([ones(len(x1)),x1,x2]).T
    coeff = linalg.lstsq(A, y)[0]
    yhat = coeff[0]+coeff[1]*x1+coeff[2]*x2
    return yhat


# read the field data
book = xlrd.open_workbook('/home/tomer/FieldData/SoilMoisture.xls')
sheet = book.sheet_by_name('Sheet1')
SM=zeros((50,9))
for i in range(50): 
    for j in range(9):
        try:
            SM[i][j] = sheet.cell_value(i+1,j+13)    
        except:
            SM[i][j]=NaN
book = xlrd.open_workbook('/home/tomer/FieldData/SoilRough.xls')
sheet = book.sheet_by_name('Sheet1')
Soil=zeros((50,11))
for i in range(50): 
    for j in range(11):
        try:
            Soil[i][j] = sheet.cell_value(i+1,j+1)    
        except:
            Soil[i][j]=NaN
# read the co-ordinates
book = xlrd.open_workbook('/home/tomer/RADARSAT/scripts/Co-ord.xls')
sheet = book.sheet_by_name('Sheet1')
Coord=zeros((29588,2))
for i in range(29588): 
    for j in range(2):
        try:
            Coord[i][j] = sheet.cell_value(i+1,j+2)    
        except:
            Coord[i][j]=NaN
sheet = book.sheet_by_name('Sheet2')
ind_pre=zeros((51,1))
for i in range(51): 
    for j in range(1):
        try:
            ind_pre[i][j] = sheet.cell_value(i+1,j)    
        except:
            ind_pre[i][j]=NaN

# perform the analysis on all the images
dirName=os.listdir("/home/tomer/RADARSAT")
dirName.remove('scripts')
dirName.remove('documents')

file = dirName[2]

# assign the corresponding indices of file to soil moisture data to ind
if '04Mar10' in file: ind=4
if '08Feb10' in file: ind=2
if '15Jan10' in file: ind=1
if '15May10' in file: ind=6
if '21Apr10' in file: ind=5
if '22Dec09' in file: ind=0

# open the files
fidSigma = gdal.Open("/home/tomer/RADARSAT/"+file+"/filtered",GA_ReadOnly)
fidIA=gdal.Open("/home/tomer/RADARSAT/"+file+"/IncidenceAngle",GA_ReadOnly)

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
IA = fidIA.GetRasterBand(1).ReadAsArray()

# read the GCP information
ImageLine=zeros(len(gcps))
ImagePixel=zeros(len(gcps))
ImageX=zeros(len(gcps))
ImageY=zeros(len(gcps))
for i in range(len(gcps)):
    a=gcps[i]
    ImageLine[i]=a.GCPLine
    ImagePixel[i]=a.GCPPixel
    ImageX[i]=a.GCPX
    ImageY[i]=a.GCPY

# convert lat, lon into Image co-ordinates
A = vstack([ones(len(gcps)),ImageX,ImageY,ImageX**2,ImageY**2,ImageX*ImageY]).T
cL = linalg.lstsq(A, ImageLine)[0]
cP = linalg.lstsq(A, ImagePixel)[0]
Lat=Coord[:,1]
Lon=Coord[:,0]
IP=cP[0]+cP[1]*Lon+cP[2]*Lat+cP[3]*Lon**2+cP[4]*Lat**2+cP[5]*Lat*Lon
IL=cL[0]+cL[1]*Lon+cL[2]*Lat+cL[3]*Lon**2+cL[4]*Lat**2+cL[5]*Lat*Lon

# extract the data
HHe=DataExtract(IL,IP,HH)
HVe=DataExtract(IL,IP,HV)
VHe=DataExtract(IL,IP,VH)
VVe=DataExtract(IL,IP,VV)

figure(figsize=(15,15))
yhat=regress2(SM[:,ind],HHe,HVe-VVe)
plot(SM[:,ind],yhat,'ro')  
xlim([0, 40]); ylim([0, 40])
xlabel('Measured Soil Moisture (v/v)',fontsize=20); ylabel('Predicted Soil Moisture (v/v)',fontsize=20)
savefig('/home/tomer/RADARSAT/documents/HHq.png',dpi=None, facecolor='w', edgecolor='w',pad_inches=0.05)


#subplot(2,2,1)
#plot(HHe,SM[:,ind]/100,'ro')
#hold
#xhat1=linspace(min(HHe)-1,max(HHe)+1,2)
#plot(xhat1,regressNan(SM[:,ind]/100,HHe,xhat1),'r')
#xlabel('HH (dB)',fontsize=20); ylabel('Soil Moisture (v/v)',fontsize=20)
#r=str(corrcoef(HHe,SM[:,ind])[0,1])
#text(min(xhat1),0.8*max(SM[:,ind])/100,'R = '+r[0:5],fontsize=20)
#subplot(2,2,2)
#plot(HVe,SM[:,ind]/100,'ro')
#hold
#xhat1=linspace(min(HVe)-1,max(HVe)+1,2)
#plot(xhat1,regressNan(SM[:,ind]/100,HVe,xhat1),'r')
#xlabel('HV (dB)',fontsize=20); ylabel('Soil Moisture (v/v)',fontsize=20)
#r=str(corrcoef(HVe,SM[:,ind])[0,1])
#text(min(xhat1),0.8*max(SM[:,ind])/100,'R = '+r[0:5],fontsize=20)
#subplot(2,2,3)
#plot(VHe,SM[:,ind]/100,'ro')
#hold
#xhat1=linspace(min(VHe)-1,max(VHe)+1,2)
#plot(xhat1,regressNan(SM[:,ind]/100,VHe,xhat1),'r')
#xlabel('VH (dB)',fontsize=20); ylabel('Soil Moisture (v/v)',fontsize=20)
#r=str(corrcoef(VHe,SM[:,ind])[0,1])
#text(min(xhat1),0.8*max(SM[:,ind])/100,'R = '+r[0:5],fontsize=20)
#subplot(2,2,4)
#plot(VVe,SM[:,ind]/100,'ro')
#hold
#xhat1=linspace(min(VVe)-1,max(VVe)+1,2)
#plot(xhat1,regressNan(SM[:,ind]/100,VVe,xhat1),'r')
#xlabel('VV (dB)',fontsize=20); ylabel('Soil Moisture (v/v)',fontsize=20)
#r=str(corrcoef(VVe,SM[:,ind])[0,1])
#text(min(xhat1),0.8*max(SM[:,ind])/100,'R = '+r[0:5],fontsize=20)
#savefig('/home/tomer/RADARSAT/documents/HHVV.png',dpi=None, facecolor='w', edgecolor='w',pad_inches=0.05)

#plot(linspace(1,50,50),Soil[:,6],'ro')
#xlabel('Plot No.',fontsize=20); ylabel('Correlation Length (mm)',fontsize=20)
#hold 
#plot(linspace(1,50,50),Soil[:,7],'ro')
#errorbar(linspace(1,50,50),0.5*(Soil[:,6]+Soil[:,7]),yerr=0.5*(Soil[:,6]-Soil[:,7]))
#savefig('/home/tomer/FieldData/corr.png',dpi=None, facecolor='w', edgecolor='w',pad_inches=0.05)

#plot(linspace(1,50,50),Soil[:,4],'ro')
#xlabel('Plot No.',fontsize=20); ylabel('RMSE (mm)',fontsize=20)
#hold 
#plot(linspace(1,50,50),Soil[:,5],'ro')
#errorbar(linspace(1,50,50),0.5*(Soil[:,4]+Soil[:,5]),yerr=0.5*(Soil[:,4]-Soil[:,5]))

#subplot(2,2,1)   
#plot(Soil[:,0],'ro')
#xlabel('Plot No.',fontsize=20); ylabel('Clay (%)',fontsize=20)
#subplot(2,2,2)   
#plot(Soil[:,1],'ro')
#xlabel('Plot No.',fontsize=20)
#ylabel('Silt (%)',fontsize=20)
#subplot(2,2,3)   
#plot(Soil[:,2],'ro')
#xlabel('Plot No.',fontsize=20); ylabel('Sand (%)',fontsize=20)
#subplot(2,2,4)   
#plot(Soil[:,3],'ro')
#xlabel('Plot No.',fontsize=20); ylabel('Gravel (%)',fontsize=20)
#savefig('/home/tomer/FieldData/CSSG.png',dpi=None, facecolor='w', edgecolor='w',pad_inches=0.05)


#figure(figsize=(25,15))
#subplot(2,2,1)   
#plot(Soil[:,1],'.')
#xlabel('Plot No.')
#ylabel('Soil Moisture (v/v)')

#subplot(2,2,2)   
#plot(Soil[:,2],'.')
#xlabel('Plot No.')
#ylabel('Soil Moisture (v/v)')
#text(linspace(1,50,50),Soil[:,2],'aa')
#pyplot.savefig('/home/tomer/FieldData/bp1.png',dpi=None, facecolor='w', edgecolor='w',pad_inches=0.1)

#figure(figsize=(25,5))
#a=nanmean(SM)
#ind=[0,1,2,4,5,7]
#plot(a[ind],'--or')
#xticks( arange(6), ("22/12/2009",'15/01/2010','08/02/2010','04/03/2010','21/04/2010','15/05/2010'),fontsize=20 )
#ylabel('Soil Moisture (v/v)',fontsize=20)
