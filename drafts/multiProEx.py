#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 15:18:34 2012

@author: mag
"""

from numpy import arange,sqrt, random, linalg
from datetime import datetime

from multiprocessing import Pool

import sys
import os
import os.path
import Image



#resize_factor = 2
#dest = os.getcwd()

#def resize(x):
#    try:
#        # Attempt to open an image file
#        filepath = x
#        image = Image.open(filepath)
#        print "opening", filepath, ":", x
#    except IOError, e:
#        # Report error, and then skip to the next argument
#        print "Problem opening", filepath, ":", e
#        return
#
#    print "resizing", filepath, ":", x
#    h,w = image.size
#    h,w = (int(h * resize_factor), int(w * resize_factor))
#    # Resize the image
#    image = image.resize((h,w), Image.ANTIALIAS)
#    fname = os.path.basename(filepath)
#
#    # Split our original filename into name and extension
#    (name, extension) = os.path.splitext(fname)
#
#    # Save the thumbnail as "(original_name)_thumb.png"
#    image.save(os.path.join(dest,name + '.jpg'),quality=80)
#    image = None


#global counter
#counter = 0
#
#def cb(r):
#    global counter
#    print counter, r
#    counter +=1
#
#def det(M):
#    return linalg.det(M)

if __name__ == '__main__':
    a = datetime.now()
#    core_ct = os.sysconf('SC_NPROCESSORS_ONLN')
#    pool = Pool(processes=core_ct)
    pool = Pool()
#    pool.map(resize,sys.argv[1:])
    for i in xrange(1,400):
        j = random.normal(1,1,(600,600))
        pool.applyasync(det,(j,),callback=cb)
#        pool.map_async(det,(j,),callback=cb)
    pool.close()
    pool.join()
    b = datetime.now()
    c = b - a
    print counter
    print "Time elapsed = ", c