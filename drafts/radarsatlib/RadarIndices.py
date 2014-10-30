#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 22:39:54 2010

@author: sat kumar tomer (http://civil.iisc.ernet.in/~satkumar/)
"""

# import required libraries
import os
import osgeo.gdal as gdal
from pylab import imshow,colorbar,title,savefig,cm,figure, close

# perform the analysis on all the images
cwd = '/home/tomer/RADARSAT/21Apr10/100m/'
#cwd = '/home/tomer/RADARSAT/08Feb10/'

fidSigma=gdal.Open(cwd+'filtered',gdal.GA_ReadOnly)
fidRho=gdal.Open(cwd+'CoPolCoeff',gdal.GA_ReadOnly)
HH = fidSigma.GetRasterBand(1).ReadAsArray()
HV = fidSigma.GetRasterBand(2).ReadAsArray()
VH = fidSigma.GetRasterBand(3).ReadAsArray()
VV = fidSigma.GetRasterBand(4).ReadAsArray()
rho = fidRho.GetRasterBand(1).ReadAsArray()
fidSigma = None
q=HV-VV
p=HH-VV

close()
figure(figsize=(35,25))
imshow(q,interpolation='nearest',cmap=cm.hsv)
colorbar(shrink=0.85)
title('q')
savefig(cwd+'q.png', facecolor='w', edgecolor='w',pad_inches=0.05) 

close()
figure(figsize=(35,25))
imshow(p,interpolation='nearest',cmap=cm.hsv)
colorbar(shrink=0.85)
title('p')
savefig(cwd+'p.png', facecolor='w', edgecolor='w',pad_inches=0.05) 

close()
figure(figsize=(35,25))
imshow(rho,interpolation='nearest',cmap=cm.hsv)
colorbar(shrink=0.85)
title('rho')
savefig(cwd+'rho.png', facecolor='w', edgecolor='w',pad_inches=0.05) 


#for dir in dirName:
#    if os.path.isdir(cwd+dir) and 'rawdata' not in dir:
#        print dir
#        fid=gdal.Open(os.path.join(cwd+dir,'filtered'),gdal.GA_ReadOnly)
#        HH = fid.GetRasterBand(1).ReadAsArray()
#        HV = fid.GetRasterBand(2).ReadAsArray()
#        VH = fid.GetRasterBand(3).ReadAsArray()
#        VV = fid.GetRasterBand(4).ReadAsArray()
#        fid = None
#        
#        close()
#        figure(figsize=(35,25))
#        imshow(HH,interpolation='nearest',cmap=cm.hsv)
#        colorbar(shrink=0.85)
#        title('HH (' + dir +' )')
#        savefig(os.path.join(cwd+dir,'HH.png'), facecolor='w', edgecolor='w',pad_inches=0.05)   
#        
#        close()
#        figure(figsize=(35,25))
#        imshow(HV,interpolation='nearest',cmap=cm.hsv)
#        colorbar(shrink=0.85)
#        title('HV (' + dir +' )')
#        savefig(os.path.join(cwd+dir,'HV.png'), facecolor='w', edgecolor='w',pad_inches=0.05)   
#        
#        close()
#        figure(figsize=(35,25))
#        imshow(VH,interpolation='nearest',cmap=cm.hsv)
#        colorbar(shrink=0.85)
#        title('VH (' + dir +' )')
#        savefig(os.path.join(cwd+dir,'VH.png'), facecolor='w', edgecolor='w',pad_inches=0.05)   
#        
#        close()
#        figure(figsize=(35,25))
#        imshow(VV,interpolation='nearest',cmap=cm.hsv)
#        colorbar(shrink=0.85)
#        title('VV (' + dir +' )')
#        savefig(os.path.join(cwd+dir,'VV.png'), facecolor='w', edgecolor='w',pad_inches=0.05)   
#
#    r='15m'# 15 30 60 90 300
#    ensure_dir(cwd+dir,r)
#    Ifile = cwd+dir+'/filtered'
#    Ofile = cwd+dir+'/' + r + '/filtered'
#    T1 = cwd+dir+'/temp1.tif'
#    T2 = cwd+dir+'/temp2.tif'
#    
#    if os.path.exists(T1):
#        os.remove(T1)
#    if os.path.exists(T2):
#        os.remove(T2)
#        
#    # convert from geographical co-ordinates to UTM zone 43
#    SC= 'gdalwarp -r lanczos -t_srs  \'+proj=utm +zone=43 +datum=WGS84\' '+Ifile+' '+T1
#    returncode = call(SC, shell=True)
#    
#    Ber = ' 664000 1309000 685000 1294000 '
#    # cut the area around Berambadi watershed
#    SC='gdal_translate -a_ullr' +Ber+ '-projwin' + Ber + T1+ ' ' + T2
#    returncode = call(SC, shell=True)
#
#    # changing the resolution 
#    SC= 'gdalwarp -r lanczos -tr ' +r[:-1] +' '+ r[:-1] +' ' + T2 + ' ' +Ofile
#    returncode = call(SC, shell=True)
#    print dir
#    
#    if os.path.exists(T1):
#        os.remove(T1)
#    if os.path.exists(T2):
#        os.remove(T2)

