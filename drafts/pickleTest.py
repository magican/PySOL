#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun May 13 23:38:40 2012

@author: mag
"""

#import pickle

import cPickle as pickle

output = open('data.pkl', 'wb')

# Pickle dictionary using protocol 0.
pickle.dump(SigmaHHwnr, output)

# Pickle the list using the highest protocol available.
pickle.dump(selfref_list, output, -1)

output.close()


import pickle

pkl_file = open('data.pkl', 'rb')

data1 = pickle.load(pkl_file)
data2 = pickle.load(pkl_file)

pkl_file.close()

from osgeo import gdal

driver = gdal.GetDriverByName('GTiff')
output_dataset = driver.Create('sigma.tiff', \
    calibPar.RasterXSize, calibPar.RasterYSize, 4, gdal.GDT_Float64)
output_dataset.SetGeoTransform(calibPar.geotransform)
output_dataset.SetGCPs(calibPar.gcps, calibPar.gcpproj)
output_dataset.GetRasterBand(1).WriteArray(calibPar.SigmaHHwnr, 0, 0)
output_dataset.GetRasterBand(2).WriteArray(calibPar.SigmaVVwnr, 0, 0)
output_dataset.GetRasterBand(3).WriteArray(calibPar.SigmaHVwnr, 0, 0)
output_dataset.GetRasterBand(4).WriteArray(calibPar.SigmaVHwnr, 0, 0)
output_dataset = None

fid=gdal.Open('sigma.tiff',gdal.GA_ReadOnly)
HH = fid.GetRasterBand(1).ReadAsArray()
HV = fid.GetRasterBand(2).ReadAsArray()
VH = fid.GetRasterBand(3).ReadAsArray()
VV = fid.GetRasterBand(4).ReadAsArray()
fid = None