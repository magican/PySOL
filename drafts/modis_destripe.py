#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 11:54:59 2012

@author: mag
"""

import datetime

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 10, 2)
__modified__ = datetime.datetime(2012, 10, 2)
__version__  = "1.0"
__status__   = "Development"


from numpy import flipud, rot90, array, zeros, linspace, histogram, cumsum, \
                  where, interp, min, max, fmin, fmax
import h5py
from scipy.io import savemat

def modis_destripe(img, resolution, nbins=100):
    """Returns the destripped image of modis band image im.
    Resolution is the resolution of image im: 1KM, HKM, or QKM.
    """
    res2nsens = {'1KM':10, 'HKM':20, 'QKM':40}
    if resolution not in res2nsens:
        raise ValueError("invalid resolution %s" % resolution)
    nsens = res2nsens[resolution]
    destriped = img.copy()
    lastrow = (img.shape[0]//nsens)*nsens
    goodsen = 0
    strippedsens = range(1, nsens)  # 0 is a good sensor
    D1 = img[goodsen:lastrow:nsens, :]
    m = D1.min()
    M = D1.max()
    H = zeros((2, nbins))
    for sen in strippedsens:
        D = img[sen:lastrow:nsens, :]
        x = linspace(fmin(m, min(D)), fmax(M, max(D)), nbins+1)
        WH, _  = histogram(D1.ravel(), bins=x, normed=True)
        H[0,:] = cumsum(WH)
        WH, _  = histogram(D.ravel(), bins=x, normed=True)
        H[1,:] = cumsum(WH)
        y = zeros(nbins)
        for i in xrange(nbins):
            indL = where(H[0,:] <= H[1,i])[0]
            if len(indL) == 0 or len(indL) == nbins:
                y[i] = x[i]
            else:
                pos = indL.max()
                xL = x[pos]
                fL = H[0, pos]
                xR = x[pos+1]
                fR = H[0, pos+1]
                y[i] = xL + (H[1,i]-fL)*(xR-xL)/float(fR-fL)
        Dall = img[sen::nsens, :]
        B = interp(Dall.ravel(), x[:-1], y)
        destriped[sen::nsens, :] = B.reshape(Dall.shape)
    return destriped
    
pn = '/home/mag/'
expFn = 'B.mat'
f = h5py.File(pn + expFn)
img = f["B"]
img = flipud(rot90(array(img)))

destriped = modis_destripe(img, 'QKM', 100)

pn = '/home/mag/'
expFn = 'destriped.mat'
savemat(pn + expFn, mdict={ \
        'destriped':destriped, \
        }, do_compression=False)