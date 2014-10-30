#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 17 15:48:21 2012

@author: mag
"""

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime

from numpy import sqrt

from createMapsEtopo1 import  findSubsetIndices

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 5, 17)
__modified__ = datetime.datetime(2012, 5, 17)
__version__  = "1.0"
__status__   = "Development"

import pygrib
from os import system, path

def ncepgfs(lonStart=3, lonEnd=4, latStart=41, latEnd=42, \
            m=None, name='ncepgfs', contour=None):
    """
    Plot latest NCEP GFS field from \
    http://nomad1.ncep.noaa.gov/pub/gfs_master/
    """

    base = '/media/SOLabNFS/SERVERS/media/hyrax/data/auxdata/model/ncep/gfs/'
    outdir = '/home/mag/'
    
    fn = 'gfs20101218/gfs.t18z.master.grbf00'
#    fn = 'gfs20101218/gfs.t18z.master.grbf03'
    
    refDate = datetime.datetime(1978,1,1,0,0,0)
    wantedDate = datetime.datetime(2010,12,18,17,39,0)

    # Get the file date
    grbs = pygrib.open(base + fn)
    grb = grbs.message(1)
    fileDate = grb.analDate
    filetime = fileDate.strftime("%Y%m%d_%H%M")
    print "File time: ", filetime
    print "File date: ", fileDate
    
    lats, lons = grb.latlons()
    u10 = grbs.message(9)['values']
    v10 = grbs.message(9)['values']
#    u = grbs.message(3)['values']
#    v = grbs.message(4)['values']
    w = sqrt(u10**2+v10**2)
    

if __name__ == "__main__":
    main()