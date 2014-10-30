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
from PIL import Image
import pylab

def pca(X):
  # Principal Component Analysis
  # input: X, matrix with training data as flattened arrays in rows
  # return: projection matrix (with important dimensions first),
  # variance and mean

  #get dimensions
  num_data,dim = X.shape

  #center data
  mean_X = X.mean(axis=0)
  for i in range(num_data):
      X[i] -= mean_X

  if dim>100:
      #print 'PCA - compact trick used'
      M = dot(X,X.T) #covariance matrix
      e,EV = linalg.eigh(M) #eigenvalues and eigenvectors
      tmp = dot(X.T,EV).T #this is the compact trick
      V = tmp[::-1] #reverse since last eigenvectors are the ones we want
      S = sqrt(e)[::-1] #reverse since eigenvalues are in increasing order
  else:
      print 'PCA - SVD used'
      U,S,V = linalg.svd(X)
      V = V[:num_data] #only makes sense to return the first num_data

  #return the projection matrix, the variance and the mean
  return V,S,mean_X

# perform the PCA analysis on all the images
import os
dirName=os.listdir("/home/tomer/RADARSAT")
dirName.remove('scripts')
dirName.remove('documents')

for file in dirName:
    # open the files
    fidSigma = gdal.Open("/home/tomer/RADARSAT/"+file+"/filtered",GA_ReadOnly)
    fidIA=gdal.Open("/home/tomer/RADARSAT/"+file+"/IncidenceAngle",GA_ReadOnly)
    if fidSigma is None:
        print 'file does not exist'

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
    
    # scale the data (mean=0) and (std=1)
    HH=HH-numpy.mean(HH)
    HH=HH/numpy.std(HH)
    HV=HV-numpy.mean(HV)
    HV=HV/numpy.std(HV)
    VH=VH-numpy.mean(VH)
    VH=VH/numpy.std(VH)
    VV=VV-numpy.mean(VV)
    VV=VV/numpy.std(VV)

    #im = numpy.array(Image.open(imlist[0])) #open one image to get the size
    m = HH.shape[0] #get the size of the images
    n = HH.shape[1] #get the size of the images
    imnbr = 4 #get the number of images
    
    #create matrix to store all flattened images
    immatrix = numpy.array([HH.flatten(),HV.flatten(),VH.flatten(),VV.flatten()])
    # numpy.array([numpy.array(Image.open(imlist[i])).flatten() for i in range(imnbr)],'f')
    
    #perform PCA
    V,S,immean = pca(immatrix)
    
    # save the PCA data as Geotiff
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create("/home/tomer/RADARSAT/"+file+"/PCASigma",HH.shape[1],HH.shape[0],4,gdal.GDT_Float32)
    output_dataset.SetGeoTransform(geotransform)
    output_dataset.SetGCPs(gcps, gcpproj)
    output_dataset.GetRasterBand(1).WriteArray(V[0].reshape(m,n), 0, 0)
    output_dataset.GetRasterBand(2).WriteArray(V[1].reshape(m,n), 0, 0)
    output_dataset.GetRasterBand(3).WriteArray(V[2].reshape(m,n), 0, 0)
    output_dataset.GetRasterBand(4).WriteArray(V[3].reshape(m,n), 0, 0)
    output_dataset = None
    
    # print the prcocessing
    print file+ " is done"