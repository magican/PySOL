#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 22:39:54 2010

@author: sat kumar tomer (http://civil.iisc.ernet.in/~satkumar/)
"""

# import required libraries
import os, xlrd
from subprocess import call

# perform the analysis on all the images
cwd = '/home/tomer/RADARSAT/'
dirName=os.listdir(cwd)
dirName.remove('scripts')
dirName.remove('documents')

def ensure_dir(f,sd):
    d = os.listdir(f)
    if sd not in d:
        os.makedirs(f+'/'+sd)

for dir in dirName:
    r='100m'# 15 30 60 90 300
    ensure_dir(cwd+dir,r)
    Ifile = cwd+dir+'/filtered'
    Ofile = cwd+dir+'/' + r + '/filtered'
    T1 = cwd+dir+'/temp1.tif'
    T2 = cwd+dir+'/temp2.tif'
    
    if os.path.exists(T1):
        os.remove(T1)
    if os.path.exists(T2):
        os.remove(T2)
        
    # convert from geographical co-ordinates to UTM zone 43
    SC= 'gdalwarp -r lanczos -t_srs  \'+proj=utm +zone=43 +datum=WGS84\' '+Ifile+' '+T1
    returncode = call(SC, shell=True)
    
    Ber = ' 664000 1309000 685000 1294000 '
    # cut the area around Berambadi watershed
    SC='gdal_translate -a_ullr' +Ber+ '-projwin' + Ber + T1+ ' ' + T2
    returncode = call(SC, shell=True)

    # changing the resolution 
    SC= 'gdalwarp -r lanczos -tr ' +r[:-1] +' '+ r[:-1] +' ' + T2 + ' ' +Ofile
    returncode = call(SC, shell=True)
    print dir
    
    if os.path.exists(T1):
        os.remove(T1)
    if os.path.exists(T2):
        os.remove(T2)
        
    # ********************* Incidence Angle*************************************
    Ifile = cwd+dir+'/IncidenceAngle'
    Ofile = cwd+dir+'/' + r + '/IncidenceAngle'
    T1 = cwd+dir+'/temp1.tif'
    T2 = cwd+dir+'/temp2.tif'
    
    if os.path.exists(T1):
        os.remove(T1)
    if os.path.exists(T2):
        os.remove(T2)
        
    # convert from geographical co-ordinates to UTM zone 43
    SC= 'gdalwarp -r lanczos -t_srs  \'+proj=utm +zone=43 +datum=WGS84\' '+Ifile+' '+T1
    returncode = call(SC, shell=True)
    
    Ber = ' 664000 1309000 685000 1294000 '
    # cut the area around Berambadi watershed
    SC='gdal_translate -a_ullr' +Ber+ '-projwin' + Ber + T1+ ' ' + T2
    returncode = call(SC, shell=True)

    # changing the resolution 
    SC= 'gdalwarp -r lanczos -tr ' +r[:-1] +' '+ r[:-1] +' ' + T2 + ' ' +Ofile
    returncode = call(SC, shell=True)
    print dir
    
    if os.path.exists(T1):
        os.remove(T1)
    if os.path.exists(T2):
        os.remove(T2)