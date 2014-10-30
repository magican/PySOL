#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 23:49:30 2012

@author: mag
"""

import matplotlib.pyplot as plt

import g1sst
reload(g1sst)
from g1sst import g1sst
import createMapsEtopo1
reload(createMapsEtopo1)
from createMapsEtopo1 import makeMap

ll_lon=2.5
ur_lon=3
ll_lat=41
ur_lat=42

plt.close('all')

g1sst(ll_lon, ur_lon, \
        ll_lat, ur_lat, m=None, name='g1sst', contour=None)
makeMap(ll_lon, ur_lon, \
        ll_lat, ur_lat, m=None, name='g1sst', contour='land')
makeMap(ll_lon, ur_lon, \
        ll_lat, ur_lat, m=None, name='g1sst', contour='ocean')