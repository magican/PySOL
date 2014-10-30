#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 18:25:37 2012

@author: mag
"""

import datetime

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 9, 10)
__modified__ = datetime.datetime(2012, 9, 10)
__version__  = "1.0"
__status__   = "Development"


from numpy import log10, absolute, zeros, reshape
from math import asin
from osgeo import gdal
from xml.dom import minidom
from scipy.signal import medfilt2d, wiener
from matplotlib.mlab import find

from multiprocessing import Process, Queue

# Filter warnings
import warnings
warnings.filterwarnings(action = "ignore", category=RuntimeWarning)
warnings.filterwarnings(action = "ignore", category=FutureWarning)

#add an argument to function for assigning a queue
def sigmaHH(dataset, xOff=0, yOff=0, xS=None, yS=None, \
    xBufScale=None, yBufScale=None, s='abs', que=[]):
    # calculate the sinclair matrix
    S_HH = dataset.GetRasterBand(1).ReadAsArray(xoff=xOff, yoff=yOff, \
    win_xsize=xS, win_ysize=yS, \
    buf_xsize=xBufScale, buf_ysize=yBufScale)
    # calculate the magnitude
    S_HH_ABS = absolute(S_HH)
    # calculate the linear sigma_naught
    if s == 'abs':
        SigmaHH = S_HH_ABS**2
    # calculate the sigma_naught in dB
    if s == 'sigma0':
        SigmaHH = 2*10*log10(S_HH_ABS)
    #we're putting return value into queue
    que.put(SigmaHH)

def __init__(self, inpath, xOff=0, yOff=0, xS=None, yS=None, \
xBufScale=None, yBufScale=None, s='abs'):
    """Calibrating the RS2 image"""

    self.inpath = inpath

        core_ct = sysconf('SC_NPROCESSORS_ONLN')
        print "Running %s cores " % core_ct
        #create a queue object
        queue1 = Queue()
        queue2 = Queue()
        queue3 = Queue()
        queue4 = Queue()
        hilo1 = Process(target=sigmaHH, \
                args=(dataset, xOff, yOff, xS, yS, \
                xBufScale, yBufScale, s, queue1))
        hilo2 = Process(target=sigmaVV, \
                args=(dataset, xOff, yOff, xS, yS, \
                xBufScale, yBufScale, s, queue2))
        hilo3 = Process(target=sigmaHV, \
                args=(dataset, xOff, yOff, xS, yS, \
                xBufScale, yBufScale, s, queue3))
        hilo4 = Process(target=sigmaVH, \
                args=(dataset, xOff, yOff, xS, yS, \
                xBufScale, yBufScale, s, queue4))
        print "Launching threads"
        hilo1.start()
        hilo2.start()
        hilo3.start()
        hilo4.start()
        #and we're assigning the return value to sigma
        self.SigmaHH = queue1.get()
        self.S_HH = queue1.get()
        self.SigmaVV = queue2.get()
        self.S_VV = queue2.get()
        self.SigmaHV = queue3.get()
        self.SigmaVH = queue4.get()
        hilo1.join()
        hilo2.join()
        hilo3.join()
        hilo4.join()