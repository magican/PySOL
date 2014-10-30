#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 16:19:03 2012

@author: mag
"""

from multiprocessing import Pool
from time import time
# time.clock for cputime works incorrectly in the example,
# it doesn't take into account other threads
N = 5000
K = 5000
def cb(r): #optional: callback function
    pass
    #print r
def CostlyFunction(z):
    r = 0
    for k in xrange(1, K+2):
        r += z ** (1 / k**1.5)
    return r
currtime = time()
for i in xrange(N):
    CostlyFunction(i)
print 'serial: time elapsed:', time() - currtime
currtime = time()
po = Pool()
for i in xrange(N):
    po.apply_async(CostlyFunction,(i,),callback=cb)
po.close()
po.join()
print '2: parallel: time elapsed:', time() - currtime