# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 18:55:33 2012

@author: mag
"""

from osgeo import osr, gdal

infile = '/home/mag/data/OTHER/RS2 Agulhas and Lion/RS2_FQA_1xQGSS20101218_173930_00000005/'

# get the existing coordinate system
ds = gdal.Open("RADARSAT_2_CALIB:SIGMA0:" + infile + "product.xml")
old_cs= osr.SpatialReference()
old_cs.ImportFromWkt(ds.GetProjectionRef())

# create the new coordinate system
wgs84_wkt = """
GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.01745329251994328,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]]"""
new_cs = osr.SpatialReference()
new_cs .ImportFromWkt(wgs84_wkt)

# create a transform object to convert between coordinate systems
transform = osr.CoordinateTransformation(old_cs,new_cs) 

#get the point to transform, pixel (0,0) in this case
width = ds.RasterXSize
height = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5] 

#get the coordinates in lat long
latlong = transform.TransformPoint(minx,miny) 




  # Calculates latitude and longitude for the center of a given pixel in the image.
  # From https://svn.osgeo.org/gdal/trunk/gdal/swig/python/samples/tolatlong.py
from osgeo import osr, gdal
from osgeo import *

infile = '/home/mag/data/OTHER/RS2 Agulhas and Lion/RS2_FQA_1xQGSS20101218_173930_00000005/'
pixel = 1
line = 1

# get the existing coordinate system
ds = gdal.Open("RADARSAT_2_CALIB:SIGMA0:" + infile + "product.xml")

ds = gdal.Open(infile + "imagery_HH.tif")

ds = gdal.Open("/home/mag/RS2-SLC-FQ14W-ASC-18-Dec-2010_17.tif")

# Read geotransform matrix and calculate ground coordinates
geomatrix = ds.GetGeoTransform()
X = geomatrix[0] + geomatrix[1] * pixel + geomatrix[2] * line
Y = geomatrix[3] + geomatrix[4] * pixel + geomatrix[5] * line

# Shift to the center of the pixel
X += geomatrix[1] / 2.0
Y += geomatrix[5] / 2.0

# Build Spatial Reference object based on coordinate system, fetched from the
# opened dataset
old_cs = osr.SpatialReference()
old_cs.ImportFromWkt(ds.GetProjection())

old_csLatLong = old_cs.CloneGeogCS()
ct = osr.CoordinateTransformation(old_cs, old_csLatLong)
(lon, lat, height) = ct.TransformPoint(X, Y)
