#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 21 19:02:23 2012

@author: mag
"""

oilCropNorthKVV = imcrop(ptsCropNorth, SigmaVVwnr)
KVV = oilCropNorthKVV[markersN==0].mean()/oilCropNorthKVV[markersN==1].mean()

oilCropNorthKHH = imcrop(ptsCropNorth, SigmaHHwnr)
KHH = oilCropNorthKHH[markersN==0].mean()/oilCropNorthKHH[markersN==1].mean()

Psleak = KHH/KVV*oilCropNorthKHH[markersN==1].mean()/oilCropNorthKVV[markersN==1].mean()