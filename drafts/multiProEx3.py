#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 16:19:03 2012

@author: mag
"""

from multiprocessing import Pool
from time import time
import os
import os.path
import Image

# time.clock for cputime works incorrectly in the example,
# it doesn't take into account other threads

resize_factor = 0.05
dest = os.getcwd()

def resize(x):
    try:
        # Attempt to open an image file
        filepath = x
        image = Image.open(filepath)
        print "opening", filepath, ":", x
    except IOError, e:
        # Report error, and then skip to the next argument
        print "Problem opening", filepath, ":", e
        return

    print "resizing", filepath, ":", x
    h,w = image.size
    h,w = (int(h * resize_factor), int(w * resize_factor))
    # Resize the image
    image = image.resize((h,w), Image.ANTIALIAS)
    fname = os.path.basename(filepath)

    # Split our original filename into name and extension
    (name, extension) = os.path.splitext(fname)

    # Save the thumbnail as "(original_name)_thumb.png"
    image.save(os.path.join(dest,name + '.jpg'),quality=80)
    image = None

currtime = time()
resize('/home/mag/Documents/repos/solab/PySOL/drafts/gebco.tif')
print 'serial: time elapsed: ', time() - currtime

currtime = time()
core_ct = os.sysconf('SC_NPROCESSORS_ONLN')
po = Pool(processes=core_ct)
po.apply_async(resize,('/home/mag/Documents/repos/solab/PySOL/drafts/gebco.tif',))
po.close()
po.join()
print '%i: parallel: time elapsed: %f' % ( core_ct, time() - currtime )