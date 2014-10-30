#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 01 23:21:47 2010

@author: sat kumar tomer (http://civil.iisc.ernet.in/~satkumar/)
"""

# import required libraries
import os, xlrd
import osgeo.gdal as gdal
from pylab import plot, xlim, ylim,title,savefig,xlabel,ylabel,figure, close, legend
import numpy as np

def utm2image(GT,utm):
    Xpixel = ((utm[:,0] - GT[0])*GT[5] - (utm[:,1] - GT[3])*GT[2])/(GT[1]*GT[5]-GT[4]*GT[2])
    Ypixel = ((utm[:,1] - GT[3])*GT[1] - (utm[:,0] - GT[0])*GT[4])/(GT[1]*GT[5]-GT[4]*GT[2])
    return (np.round(Xpixel)).astype('int'),(np.round(Ypixel)).astype('int')

# extract data from images
def DataExtract(IL,IP,data):
    DataE = np.zeros(len(IL))
    DataE = data[IL,IP]
    Datae = np.zeros(50)
    for i in range(50):
        Datae[i] = np.median(DataE[int(ind_pre[i]):int(ind_pre[i+1])])
    return Datae

def nans(shape, dtype=float):
    a = np.empty(shape, dtype)
    a.fill(np.nan)
    return a

# read the co-ordinates
book = xlrd.open_workbook('/home/tomer/RADARSAT/scripts/Co-ord.xls')
sheet = book.sheet_by_name('Sheet1')
Coord = np.zeros((29588,2))
for i in range(29588): 
    for j in range(2):
        Coord[i,j] = sheet.cell_value(i+1,j)    
sheet = book.sheet_by_name('Sheet2')
ind_pre = np.zeros((51,))
for i in range(51): 
    ind_pre[i] = sheet.cell_value(i+1,0) 
    
# read the field data
book = xlrd.open_workbook('/home/tomer/FieldData/SoilMoisture.xls')
sheet = book.sheet_by_name('Sheet1')
SM = np.zeros((50,9))
jj = range(9)
jj.remove(3); jj.remove(6); jj.remove(8)

for i in range(50): 
    for j in range(9):
        try:
            SM[i][j] = sheet.cell_value(i+1,j+13)    
        except:
            SM[i][j] = np.nan
SM = SM[:,jj]

# read the coefficient from the Envisat analysis
book = xlrd.open_workbook('/home/tomer/ENVISAT/RegreCoeff.xls')
sheet = book.sheet_by_name('Sheet1')
coeff = np.zeros((24,3))
for i in range(24): 
    for j in range(3):
        coeff[i,j] = sheet.cell_value(i+1,j)    

# perform the analysis on all the images
cwd = '/home/tomer/RADARSAT/'
dirs = os.listdir(cwd)
dirs.remove('scripts')
dirs.remove('documents')

sigma_HH = nans((50,6))

for Dir in dirs:
    fid = gdal.Open(cwd+Dir+'/100m/filtered',gdal.GA_ReadOnly)
    HH = fid.GetRasterBand(1).ReadAsArray()
    GT = fid.GetGeoTransform()
    gcps = fid.GetGCPs()
    gcpproj = fid.GetGCPProjection()
    fid = None
    x,y = utm2image(GT,Coord)
    HHe=DataExtract(y,x,HH)
    #print Dir + ' is done'
    
    if '04Mar10' in Dir: ind=3
    if '08Feb10' in Dir: ind=2
    if '15Jan10' in Dir: ind=1
    if '15May10' in Dir: ind=5
    if '21Apr10' in Dir: ind=4
    if '22Dec09' in Dir: ind=0
    
    sigma_HH[:,ind] = HHe
    
figure(figsize=(30,25))
##PlotNo = 10
##PlotNo = 40
##PlotNo = 2
##x = (-7,-10)
#x = (-5,-12)
##x = (-11,-14)
#PlotNo = (40,20)
#y = coeff[coeff[:,0]==PlotNo,1] + x*coeff[coeff[:,0]==PlotNo,2]
#plot(sigma_HH[PlotNo-1,:]-0.12*(20-23),SM[PlotNo-1,:],'ro', label='RADARSAT obs.')
#plot(x,y, label='ENVISAT fitted')
#xlabel('BC_HH (dB)',fontsize=20)
#ylabel('volumetric SM (%)',fontsize=20)
#legend(loc=4)
#ylim(ymin=0)
#title('Plot No. ' + str(PlotNo))
##savefig('/home/tomer/RADARSAT/documents/'+str(PlotNo) +'.png',dpi=None, facecolor='w', edgecolor='w',pad_inches=0.05)


# figure (all in one plot)
#PlotNo = (coeff[:,0]).astype('int')
PlotNo = np.array([1,2,35,36,8,9,10,12,26,27,33,37,38,40,18,19,24,25,29,30,31,32])
#PlotNo = np.array([1,2])
for i in range(len(PlotNo)):
    if i<7: 
        sym = 'd'; MS=9
    elif i<14: 
        sym = 's'; MS=8
    elif i<21:
        sym = 'p'; MS=9
    else: 
        sym = 'h'; MS=9
    plot(sigma_HH[PlotNo[i]-1,:]-0.12*(20-23),SM[PlotNo[i]-1,:],sym, ms=MS, label=str(PlotNo[i]))
xlabel('BC_HH (dB)',fontsize=20)
ylabel('volumetric SM (%)',fontsize=20)
legend(loc=4)
ylim(ymin=0,ymax=45)
xlim((-14,-2))
title('Plot No. ' + str(PlotNo))
savefig('/home/tomer/RADARSAT/documents/all.png',dpi=None, facecolor='w', edgecolor='w',pad_inches=0.05)