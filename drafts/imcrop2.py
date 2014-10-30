#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 11:34:20 2012

@author: mag
"""

import numpy as np
import matplotlib.pyplot as plt

import imcrop
reload(imcrop)
from imcrop import imcrop, imzoom

plt.figure()
image = np.random.randn(128,128)
plt.imshow(image)
pts = imzoom()
imCrop = imcrop(pts, image)


import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
for i, label in enumerate(('A', 'B', 'C', 'D')):
    ax = fig.add_subplot(2,2,i+1)
    ax.text(0.05, 0.95, label, transform=ax.transAxes,
      fontsize=16, fontweight='bold', va='top')

plt.show()