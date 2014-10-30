#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 18:15:10 2012

@author: mag
"""

import os
import osgeo.gdal as gdal
from pylab import imshow,colorbar,title,savefig,cm,figure, close
from subprocess import call

pn = '/media/data/data/OTHER/RS2\ Agulhas\ and\ Lion/RS2_FQA_1xQGSS20101218_173930_00000005/'

Ifile = pn + 'sigma.tiff'

T1 = pn + 'temp1.tif'
T2 = pn + 'temp2.tif'

if os.path.exists(T1):
    os.remove(T1)
    
if os.path.exists(T2):
    os.remove(T2)

# convert from geographical co-ordinates to UTM zone 31
SC= 'gdalwarp -r lanczos -t_srs  \'+proj=utm +zone=31 +datum=WGS84\' '+Ifile+' '+T1
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
