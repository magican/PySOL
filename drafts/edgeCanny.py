"""
===================
Canny edge detector
===================

The Canny filter is a multi-stage edge detector. It uses a filter based on the
derivative of a Gaussian in order to compute the intensity of the gradients.The
Gaussian reduces the effect of noise present in the image. Then, potential
edges are thinned down to 1-pixel curves by removing non-maximum pixels of the
gradient magnitude. Finally, edge pixels are kept or removed using hysteresis
thresholding on the gradient magnitude.

The Canny has three adjustable parameters: the width of the Gaussian (the
noisier the image, the greater the width), and the low and high threshold for
the hysteresis thresholding.
"""

import matplotlib.pyplot as plt

import datetime

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 5, 21)
__modified__ = datetime.datetime(2012, 5, 21)
__version__  = "1.0"
__status__   = "Development"


from skimage.filter import threshold_otsu
from skimage.filter import median_filter

ptsCropSouth = array([[ 3000, 4650],
                      [ 3700, 5200]])
oilCropSouth = imcrop(ptsCropSouth, delta)
ptsCropNorth = array([[ 3400, 0],
                      [ 4900,  1650]])
oilCropNorth = imcrop(ptsCropNorth, delta)

image = oilCropNorth

# reduce the speckle noise
image = wiener(image, mysize=(3,3), noise=None)

# use median filter to smooth
image = median_filter(image, radius=37)

# find the global threshold
global_thresh = threshold_otsu(image)

# create a binary mask
markers = ones(image.shape, dtype=uint)
markers[image < global_thresh*0.84] = 0

oilCropNorthMeanP = P[markers==1]
oilCropNorthMeanP = oilCropNorthMeanP.mean()

## display results
#plt.figure()
#plt.subplot(131)
#plt.imshow(image, cmap=plt.cm.gray)
#plt.axis('off')
#plt.title('noisy image', fontsize=20)
#plt.subplot(132)
#plt.imshow(binary_global, cmap=plt.cm.gray)
#plt.axis('off')
#plt.title('Global thresholding', fontsize=20)
#plt.subplot(133)
#plt.imshow(markers, cmap=plt.cm.gray)
#plt.axis('off')
#plt.title('Markers mask', fontsize=20)
#plt.subplots_adjust(wspace=0.02, hspace=0.02, top=0.9,
#                    bottom=0.02, left=0.02, right=0.98)
#plt.show()




"""
Global thresholding
"""
from skimage.filter import threshold_otsu, threshold_adaptive
from skimage.filter import median_filter

ptsCropSouth = array([[ 3000, 4650],
                      [ 3700, 5200]])
oilCropSouth = imcrop(ptsCropSouth, delta)
ptsCropNorth = array([[ 3400, 0],
                      [ 5000,  1650]])
oilCropNorth = imcrop(ptsCropNorth, delta)

image = oilCropNorth

image = median_filter(image, radius=50)

global_thresh = threshold_otsu(image)
binary_global = image > global_thresh

block_size = 300
binary_adaptive = threshold_adaptive(image, block_size)

markers = ones(image.shape, dtype=uint)
markers[image < global_thresh*0.78] = 0


# display results
plt.figure()
plt.subplot(221)
plt.imshow(image, cmap=plt.cm.gray)
plt.axis('off')
plt.title('noisy image', fontsize=20)
plt.subplot(222)
plt.imshow(binary_global, cmap=plt.cm.gray)
plt.axis('off')
plt.title('Global thresholding', fontsize=20)
plt.subplot(223)
plt.imshow(binary_adaptive, cmap=plt.cm.gray)
plt.axis('off')
plt.title('Adaptive thresholding', fontsize=20)
plt.subplot(224)
plt.imshow(markers, cmap=plt.cm.gray)
plt.axis('off')
plt.title('Markers mask', fontsize=20)

plt.subplots_adjust(wspace=0.02, hspace=0.02, top=0.9,
                    bottom=0.02, left=0.02, right=0.98)


plt.show()


"""
Histogram Equalization
"""

from skimage.util.dtype import dtype_range
from skimage import exposure
import numpy as np

img = oilCropNorth

# Contrast stretching
pl = np.percentile(img, 5)
ph = np.percentile(img, 95)
img_rescale = exposure.rescale_intensity(img, in_range=(pl, ph))

# Equalization
img_eq = exposure.equalize(img)


# display results
plt.subplot(131)
plt.imshow(img, cmap=plt.cm.gray)
plt.axis('off')
plt.title('noisy image', fontsize=20)
plt.subplot(132)
plt.imshow(img_rescale, cmap=plt.cm.gray)
plt.axis('off')
plt.title('Contrast stretching', fontsize=20)
plt.subplot(133)
plt.imshow(img_eq, cmap=plt.cm.gray)
plt.axis('off')
plt.title('Histogram equalization', fontsize=20)

plt.subplots_adjust(wspace=0.02, hspace=0.02, top=0.9,
                    bottom=0.02, left=0.02, right=0.98)


plt.show()






import numpy as np
from skimage.transform import hough, probabilistic_hough
from skimage.filter import canny

image = oilCropNorth

# Line finding, using the Probabilistic Hough Transform
edges = canny(image, sigma=1)
lines = probabilistic_hough(edges, threshold=10, line_length=5, line_gap=3)

plt.subplot(131)
plt.imshow(image, cmap=plt.cm.gray)
plt.title('Input image')
plt.subplot(132)
plt.imshow(edges, cmap=plt.cm.gray)
plt.title('Sobel edges')

plt.subplot(133)
plt.imshow(edges * 0)

for line in lines:
    p0, p1 = line
    plt.plot((p0[0], p1[0]), (p0[1], p1[1]))

plt.title('Lines found with PHT')
plt.axis('image')
plt.show()


from skimage.segmentation import random_walker

# Generate noisy synthetic data
data = oilCropNorth
global_thresh = threshold_otsu(data)
markers = np.zeros(data.shape, dtype=np.uint)
markers[data < global_thresh*0.78] = 1
markers[data > global_thresh*0.78] = 0

# Run random walker algorithm
labels = random_walker(data, markers, beta=10, mode='bf')

# Plot results
plt.figure()
plt.subplot(131)
plt.imshow(data, cmap='gray', interpolation='nearest')
plt.axis('off')
plt.title('Noisy data')
plt.subplot(132)
plt.imshow(markers, cmap='hot', interpolation='nearest')
plt.axis('off')
plt.title('Markers')
plt.subplot(133)
plt.imshow(labels, cmap='gray', interpolation='nearest')
plt.axis('off')
plt.title('Segmentation')

plt.subplots_adjust(hspace=0.01, wspace=0.01, top=1, bottom=0, left=0,
                    right=1)
plt.show()



from skimage import filter

def imnorm(oldvalue, newmin=0, newmax=1):
    oldmin=oldvalue.min()
    oldmax=oldvalue.max()
    oldrange=oldmax-oldmin
    newrange=newmax-newmin
    #where in the old scale is this value (0...1)
    scale=(oldvalue-oldmin)/oldrange
    #place this scale in the new range
    newvalue=(newrange*scale)+newmin
    return newvalue

#The greyscale input image to detect edges on
#should be normalized to 0.0 to 1.0
im = imnorm(oilCropNwnr)

# Compute the Canny filter for two values of sigma
edges1 = filter.canny(im)
edges2 = filter.canny(im, sigma=3)

# display results
plt.figure()

plt.subplot(131)
plt.imshow(im, cmap=plt.cm.gray)
plt.axis('off')
plt.title('noisy image', fontsize=20)

plt.subplot(132)
plt.imshow(edges1, cmap=plt.cm.gray)
plt.axis('off')
plt.title('Canny filter, $\sigma=1$', fontsize=20)

plt.subplot(133)
plt.imshow(edges2, cmap=plt.cm.gray)
plt.axis('off')
plt.title('Canny filter, $\sigma=3$', fontsize=20)

plt.subplots_adjust(wspace=0.02, hspace=0.02, top=0.9,
                    bottom=0.02, left=0.02, right=0.98)


plt.show()