
# coding: utf-8

# In[25]:

# !TODO

# CMOD for HH and HV/VH - Hwang

# Ice from cross pol

# Svalbard/Barents Area Case Study

# Remove antenna pattern
# https://github.com/senbox-org/s1tbx/blob/master/s1tbx-op-calibration/src/main/java/org/esa/nest/gpf/ASARCalibrator.java

# Trimm the data by removing zero values from side rows and cols

# Read and calibrate SLC product

# Example from Fab (Chapron)
# s1a-s5-grd-hh-20140818t181248-20140818t181312-001998-001ef8-001-roughness
# https://mail.google.com/mail/u/0/#inbox/14961597bb0cef3c


# In[26]:

import zipfile

import datetime

from matplotlib.mlab import find
import xmltodict

from PIL import Image
import StringIO

from scipy.interpolate import RectBivariateSpline
from scipy.signal import wiener

import os

from numpy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2014, 10, 28)
__modified__ = datetime.datetime(2014, 11, 13)
__version__  = "1.0"
__status__   = "Development"


# In[27]:

from multiprocessing import cpu_count
numProcs = cpu_count()-2


# In[28]:

from numpy import linspace, arange

def imresize(image, size):
    """
    Resizes coefficient arrays using bivariate spline approximation.
    """
    m, n = image.shape
    X = linspace(0, m - 1, size[0])
    Y = linspace(0, n - 1, size[1])
    kx, ky = min([m - 1, 3]), min([n - 1, 3])
    interp = RectBivariateSpline(
        arange(m), arange(n), image, kx=kx, ky=ky)
    resized = interp(X, Y)
    return resized


# In[29]:

from math import atan

def windDirection(u, v):
    U = u.ravel()
    V = v.ravel()
    direction = zeros(size(U))
    for i in range(0, len(U)):
        if U[i] >= 0 and V[i] > 0: direction[i] = ((180 / pi) * atan(abs(U[i] / V[i])) + 180)
        if U[i] < 0 and V[i] > 0: direction[i] = (-(180 / pi) * atan(abs(U[i] / V[i])) + 180)
        if U[i] >= 0 and V[i] < 0: direction[i] = (-(180 / pi) * atan(abs(U[i] / V[i])) + 360)
        if U[i] < 0 and V[i] < 0: direction[i] = ((180 / pi) * atan(abs(U[i] / V[i])))
        if V[i] == 0 and U[i] > 0: direction[i] = 270
        if V[i] == 0 and U[i] < 0: direction[i] = 90
        if V[i] == 0 and U[i] == 0: direction[i] = 0
    return reshape(direction, v.shape)


# In[30]:

# Using NCEP
from createMapsEtopo1 import findSubsetIndices
import pygrib

def ncepGFSmodel(startTime, lats_2, lons_2):
    """
    NCEP GFS model wind for givven time, lat/lon crop 
    """
    ncepGFSmodel = {} # empty dict for ncepGFSmodel

    iPath_wind = '/media/SOLabNFS2/store/model/ncep/gfs/'

    # find the ncep gfs filename to open from ASAR filename
    baseHour = floor((startTime.hour+3/2)/6)*6
    baseHour = min(18, baseHour)
    if startTime.hour-baseHour>1.5:
        forecastHour = 3
    else:
        forecastHour = 0

    if startTime <= datetime.datetime(2014, 8, 19):
        ncepFileName = 'gfs' + startTime.strftime("%Y%m%d") + '/gfs.t' + '%.2d' %(baseHour) + 'z.master.grbf' + '%.2d' %(forecastHour)

        grbs = pygrib.open(iPath_wind + ncepFileName)

        u_wind = None
        v_wind = None

        # wind contains u=u_wind.values[:], Lats=u_wind.latlons()[0], Lons=u_wind.latlons()[1]
        for idx, msg_info in enumerate(grbs.select()):
            if msg_info['short_name'] == '10u':
                u_wind = grbs.message(idx + 1)
            elif msg_info['short_name'] == '10v':
                v_wind = grbs.message(idx + 1)

        u = u_wind.values[:]
        v = v_wind.values[:]
        lats_wind = u_wind.latlons()[0]
        lons_wind = u_wind.latlons()[1]
    else:
        ncepFileName = 'gfs.' + startTime.strftime("%Y%m%d") + '%.2d' %(baseHour) + '/gfs.t' + '%.2d' %(baseHour) + 'z.master.grbf' + '%.2d' %(forecastHour) + '.10m.uv.grib2'

        grbs = pygrib.open(iPath_wind + ncepFileName)

        u_wind = grbs.message(1)
        v_wind = grbs.message(2)
        u = u_wind['values']
        v = v_wind['values']
        lats_wind = u_wind['latitudes']
        lons_wind = u_wind['longitudes']
        lons_wind = reshape(lons_wind, (lons_wind.shape[0]/720, 720))
        lats_wind = reshape(lats_wind, (lats_wind.shape[0]/720, 720))


    #Make sure the longitude is between -180.00 .. 179.9
    lons_wind = map(lambda x : (lons_wind.ravel()[x]+180)-int((lons_wind.ravel()[x]+180)/360)*360-180, range(0,lons_wind.size))
    lons_wind = reshape(lons_wind, lats_wind.shape)
    # plt.close('all')
    # plt.imshow(lons_wind)
    # plt.colorbar()

#     #Make sure the latitudes is between -90.00 .. 89.9, starting from North - positive
#     lats_wind = map(lambda x : (lats_wind.ravel()[x]+90)-int((lats_wind.ravel()[x]+90)/180)*180-90, xrange(0,lats_wind.size))
#     lats_wind = reshape(lats_wind, lons_wind.shape)
#     if lats_wind[0,0] < lats_wind[-1,-1]:
#         lats_wind = flipud(lats_wind)
#         u = flipud(u)
#         v = flipud(v)
#     plt.close('all')
#     plt.imshow(lats_wind)
#     plt.colorbar()


    # find subset
    res = findSubsetIndices(lats_2.min(),lats_2.max(),lons_2.min(),lons_2.max(),lats_wind[:,0],lons_wind[0,:])
    # expand subset by 1 pixel for better further pyresample
    res[0]=res[0]-2
    res[1]=res[1]+2
    res[2]=res[2]-2
    res[3]=res[3]+2

    # crop the data
    u = u[int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    v = v[int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    ncepGFSmodel['lats_wind'] = lats_wind[int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    ncepGFSmodel['lons_wind'] = lons_wind[int(res[2]):int(res[3]),int(res[0]):int(res[1])]

    ncepGFSmodel['wind_dir'] = windDirection(u,v)
    ncepGFSmodel['wind_speed'] = sqrt(u**2 + v**2)

#     del u_wind, v_wind
    return ncepGFSmodel


# In[31]:

# READ THE RAW_COUNTS from GRD image

def read_raw_counts(fn, fileLocation, polarization):


    # Note that For SLC images BitsPerSample=32 and for GRD BitsPerSample=16

    # im = Image.open(inpath + fileLocation['s1aiwgrd' + polarization])

    im = zf.read(fn[:-4] + '.SAFE' + fileLocation[fn.lower().replace("_","")[0:8] + polarization][1:])
    im = StringIO.StringIO(im) #Encode the raw data to be used by Image.open()
    im = Image.open(im)        #Open the image

    return asarray(im)


# In[32]:

# READ the ANNOTATION

def read_anotation(fn, fileLocation, polarization):
    # open the fileLocation

    # annotation = open(inpath + fileLocation['products1aiwgrd' + polarization], "r") # Open a file in read-only mode
    # annotation = annotation.read() # read the file object

    annotation = zf.read(fn[:-4] + '.SAFE' + fileLocation['product' + fn.lower().replace("_","")[0:8] + polarization][1:])
    annotation = xmltodict.parse(annotation) # Parse the read document string

    # get geolocationGrid parameters from the Annotation Data Set Records (ADSR)
    # preallocate variables
    GEOgrid = {} # empty dict for GEOgrids

    geolocationGridPointList = annotation['product']['geolocationGrid']['geolocationGridPointList']
    GEOgrid['lats']  = zeros( ( int(geolocationGridPointList['@count']), 1) )
    GEOgrid['lons']  = zeros( GEOgrid['lats'].shape )
    GEOgrid['line']  = zeros( GEOgrid['lats'].shape, dtype=int )
    GEOgrid['pixel'] = zeros( GEOgrid['lats'].shape, dtype=int )
    GEOgrid['incidenceAngle'] = zeros( GEOgrid['lats'].shape )

    # read Geolocation grid points
    for n in range(int(geolocationGridPointList['@count'])):
        GEOgrid['lats'][n]  = float(geolocationGridPointList['geolocationGridPoint'][n]['latitude'])
        GEOgrid['lons'][n]  = float(geolocationGridPointList['geolocationGridPoint'][n]['longitude'])
        GEOgrid['line'][n]  = int(geolocationGridPointList['geolocationGridPoint'][n]['line'])
        GEOgrid['pixel'][n] = int(geolocationGridPointList['geolocationGridPoint'][n]['pixel'])
        GEOgrid['incidenceAngle'][n] = float(geolocationGridPointList['geolocationGridPoint'][n]['incidenceAngle'])


    # find zero pixel to rehape grid points to array
    ind = find(GEOgrid['pixel'] == 0)
    GEOgrid['pixel'] = reshape(GEOgrid['pixel'], (ind.size, GEOgrid['lats'].size/ind.size))
    GEOgrid['line']  = reshape(GEOgrid['line'], (ind.size, GEOgrid['lats'].size/ind.size))
    GEOgrid['lats']  = reshape(GEOgrid['lats'], (ind.size, GEOgrid['lats'].size/ind.size))
    GEOgrid['lons']  = reshape(GEOgrid['lons'], (ind.size, GEOgrid['lats'].size/ind.size))
    GEOgrid['incidenceAngle'] = reshape(GEOgrid['incidenceAngle'], (ind.size, GEOgrid['lats'].size/ind.size))
    
    return GEOgrid


# In[33]:

# READ the CALIBRATION LUTs

def read_clbrtn_luts(fn, fileLocation, polarization):
    # The calibration data set contains calibration information
    # and the beta nought, sigma nought, gamma and digital
    # number (DN) Look-up Tables (LUT)s that can be used for
    # absolute product calibration.

    # We take only the calibrationVector record
    # This record holds the calibration vectors and associated fields required to
    # derive radiometrically calibrated imagery from the image MDS.

    # open the fileLocation

    # calibration = open(inpath + fileLocation['calibrations1aiwgrd' + polarization], "r") # Open a file in read-only mode
    # calibration = calibration.read() # read the file object

    calibration = zf.read(fn[:-4] + '.SAFE' + fileLocation['calibration' + fn.lower().replace("_","")[0:8] + polarization][1:])
    calibration = xmltodict.parse(calibration) # Parse the read document string

    calibrationVectorList = calibration['calibration']['calibrationVectorList']

    cLUTs = {} # empty dict for CALIBRATION LUTs

    cLUTs['line']  = zeros( ( int(calibrationVectorList['@count']), 1), dtype=int )
    cLUTs['pixel'] = zeros( (cLUTs['line'] .shape[0], int(calibrationVectorList['calibrationVector'][0]['pixel']['@count'])), dtype=int )
    cLUTs['sigmaNought']  = zeros( cLUTs['pixel'].shape, dtype=float)

    # read Calibration Vector points
    for n in range(int(calibrationVectorList['@count'])):
        cLUTs['line'][n] = int(calibrationVectorList['calibrationVector'][n]['line'])
        pixel_     = calibrationVectorList['calibrationVector'][0]['pixel']['#text']
        cLUTs['pixel'][n,:] = asarray(pixel_.split(' '), dtype=int)
        sigmaNought_     = calibrationVectorList['calibrationVector'][0]['sigmaNought']['#text']
        cLUTs['sigmaNought'][n,:] = asarray(sigmaNought_.split(' '), dtype=float)

    return cLUTs


# In[34]:

# READ the NOISE LUTs

def read_noise_luts(fn, fileLocation, polarization):
    # The L1 Noise ADS provides a LUT – with values provided in linear power – 
    # that can be used to derive calibrated noise profiles which match the calibrated GRD data.

    # open the fileLocation

    noise = zf.read(fn[:-4] + '.SAFE' + fileLocation['noise' + fn.lower().replace("_","")[0:8] + polarization][1:])
    noise = xmltodict.parse(noise) # Parse the read document string

    noiseVectorList = noise['noise']['noiseVectorList']

    nLUTs = {} # empty dict for NOISE LUTs

    nLUTs['line']  = zeros( ( int(noiseVectorList['@count']), 1), dtype=int )
    nLUTs['pixel'] = zeros( (nLUTs['line'] .shape[0], int(noiseVectorList['noiseVector'][0]['pixel']['@count'])), dtype=int )
    nLUTs['noiseLut']  = zeros( nLUTs['pixel'].shape, dtype=float)

    # read Calibration Vector points
    for n in range(int(noiseVectorList['@count'])):
        nLUTs['line'][n] = int(noiseVectorList['noiseVector'][n]['line'])
        pixel_     = noiseVectorList['noiseVector'][0]['pixel']['#text']
        nLUTs['pixel'][n,:] = asarray(pixel_.split(' '), dtype=int)
        noiseLut_     = noiseVectorList['noiseVector'][0]['noiseLut']['#text']
        nLUTs['noiseLut'][n,:] = asarray(noiseLut_.split(' '), dtype=float)

    return nLUTs


# In[34]:




# In[35]:

# def xml2geo(self):
#         """
#         Reading geolocation grids from attributes
#         One may note that Ground Range, Multi-Look, Detected (GRD) products
#         lie in the ground range by azimuth surface,
#         with image coordinates oriented along ground range and flight direction.
#         Slant Range, Single-Look Complex (SLC) products
#         are images in the slant range by azimuth imaging plane,
#         in the image plane of satellite data acquisition.
#         """

# inpath = '/media/SOLabNFS2/tmp/sentinel-1/S1A_IW_GRDH_1SDV_20141004T155619_20141004T155644_002682_002FE5_BE58.SAFE/'
inpath = '/media/SOLabNFS2/tmp/sentinel-1/'

# READ THE MANIFEST - top level

# manifest = open(inpath + "manifest.safe", "r") # Open a file in read-only mode
# manifest = manifest.read()      # read the file object

# open zip file without extraction
# finnGulf
# fn = 'S1A_IW_GRDH_1SDV_20141004T155619_20141004T155644_002682_002FE5_BE58.zip'

# Svalbard-Barents
inpath = '/media/SOLabNFS2/tmp/sentinel-1/Svalbard-Barents/'
fn = 'S1A_EW_GRDH_1SDH_20141003T071117_20141003T071221_002662_002F6B_CAFE.zip'

fileNameList = []
for _dir, sub_dir, _files in os.walk(inpath):
    for fileName in _files:
        if fileName.startswith('S1A') and fileName.endswith('.zip') and fileName.find('RAW')==-1:
            fileNameList.append(fileName)

for fn in fileNameList:

	zf = zipfile.ZipFile(inpath+fn, 'r')

	manifest = zf.read(fn[:-4] + '.SAFE/manifest.safe')
	manifest = xmltodict.parse(manifest) # Parse the read document string

	# we have different XML/GeoTIFF fileLocations for vv-vh (or hh-hv) products, respectively
	# create new dict with file paths to the fileLocations
	fileLocation = {} # empty dict
	dataObject = manifest['xfdu:XFDU']['dataObjectSection']['dataObject']
	for n in range(len(dataObject)):
		if len(dataObject[n]['@ID']) > 45:
			k = dataObject[n]['@ID'][0:-45] # get the new key from @ID
		else:
			k = dataObject[n]['@ID'] # get the new key from @ID
		v = str(dataObject[n]['byteStream']['fileLocation']['@href']) # get the dict.value
		fileLocation[k] = v # assign to new dict
	#     locals()['fileLocation_'+k]=v # create local variable from 'fileLocation' dict.key


	# In[36]:

	# Get the productType/polarization

	metadataObject = manifest['xfdu:XFDU']['metadataSection']['metadataObject']
	for n in range(len(metadataObject)):
		if metadataObject[n]['@ID'] == 'generalProductInformation':
			productType = metadataObject[n]['metadataWrap']['xmlData']['s1sarl1:standAloneProductInformation']['s1sarl1:productType']
			transmitterReceiverPolarisation = metadataObject[n]['metadataWrap']['xmlData']['s1sarl1:standAloneProductInformation']['s1sarl1:transmitterReceiverPolarisation']

	polarization = []
	for p in range(0,len(transmitterReceiverPolarisation[0].lower())):
		polarization.append(transmitterReceiverPolarisation[p].lower())
	print "Available polarizations: \'%s\'" %polarization


	# In[37]:

	# %%timeit -n 1 -r 1

	raw_counts = {}
	lats_2 = {}
	lons_2 = {}
	incidenceAngle_2 = {}
	sigmaNought_2 = {}
	noiseLut_2 = {}
	sigma0 = {}

	for p in polarization:
		print "Reading raw_counts: \'%s\' polarization" %p

		# READ THE RAW_COUNTS from GRD image
		raw_counts[p] = read_raw_counts(fn, fileLocation, p)

	# READ the ANNOTATION
	GEOgrid = read_anotation(fn, fileLocation, polarization[0])

	# READ the CALIBRATION LUTs
	cLUTs = read_clbrtn_luts(fn, fileLocation, polarization[0])

	# READ the NOISE LUTs
	nLUTs = read_noise_luts(fn, fileLocation, polarization[0])

	#  ----------------------------------
	# INTERPOLATE DATA
	# Interpolate Geolocation grid points and calibration LUTs onto grid of raw_counts.shape

	# Serial Processing

	# Serial loop is faster than the Parallel Loop (see appropriate ipnb file)
	# 1 loops, best of 1: 31.5 s per loop
	# %%timeit -n 1 -r 1

	# create new grid of raw_counts shape
	line_2 = arange(raw_counts[p].shape[0])
	pixel_2 = arange(raw_counts[p].shape[1])

	# Interpolate onto a new grid
	lats_2 = RectBivariateSpline(GEOgrid['line'][:,0], GEOgrid['pixel'][0,:], GEOgrid['lats'], kx=2, ky=2)(line_2, pixel_2)
	lons_2 = RectBivariateSpline(GEOgrid['line'][:,0], GEOgrid['pixel'][0,:], GEOgrid['lons'], kx=2, ky=2)(line_2, pixel_2)


	for p in polarization:
		print "Interpolating LUTs: \'%s\' polarization" %p

		# interpolate incidenceAngle
		incidenceAngle_2[p] = RectBivariateSpline(GEOgrid['line'][:,0], GEOgrid['pixel'][0,:], GEOgrid['incidenceAngle'], kx=2, ky=2)(line_2, pixel_2)

		# interpolate sigmaNought
		sigmaNought_2[p] = RectBivariateSpline(cLUTs['line'][:,0], cLUTs['pixel'][0,:], cLUTs['sigmaNought'], kx=2, ky=2)(line_2, pixel_2)

		# interpolate noiseLut
		noiseLut_2[p] = RectBivariateSpline(nLUTs['line'][:,0], nLUTs['pixel'][0,:], nLUTs['noiseLut'], kx=2, ky=2)(line_2, pixel_2)

		# Apply Calibration, remove the thermal noise estimation and Convert to Intensity
		sigma0[p] = ( double(raw_counts[p])**2 - noiseLut_2[p] )/sigmaNought_2[p]**2


	# In[ ]:




	# In[ ]:




	# In[ ]:




	# In[38]:

	def swath_area_def(name='Temporal SWATH EPSG Projection 4326', proj='eqc', lonlim=(-180,180), latlim=(-90,90), ellps="WGS84", res=111.2e3, lat_ts=None, lat_0=None, lon_0=None):
		"""
		Convert given swath coordinates to pyresample area definition.
		The arguments are standard for Proj:
		name
		proj
		lonlim
		latlim
		ellipsoid
		resolution(meters)
		lat_ts (latitude of true scale)
		lat_0,lon_0 is central point
		EXAMPLE:
		
		epsg3426 is the default one
		for epsg3413:
		swath_area_def(name='Temporal SWATH EPSG Projection 3413', proj='stere', lonlim=(-180,180), latlim=(30,90), ellps="WGS84", res=111.2e3, lat_ts=70, lat_0=90, lon_0=-45)
			
		"""

		area_id = name.replace(" ", "_").lower()
		proj_id = area_id

		up    = min(latlim)
		down  = max(latlim)
		left  = min(lonlim)
		right = max(lonlim)
		
		if proj == 'eqc':
			p = Proj(proj=proj, llcrnrlat=up, urcrnrlat=down, llcrnrlon=left, urcrnrlon=right, ellps=ellps)
			proj4_args = '+proj=' + str(proj) + ' ' +              '+llcrnrlat=' + str(up) + ' ' +              '+urcrnrlat=' + str(down) + ' ' +              '+llcrnrlon=' + str(left) + ' ' +              '+urcrnrlon=' + str(right) + ' ' +              '+ellps=' + str(ellps)
		elif lat_ts!=None and lat_0!=None:
			# lat_ts is latitude of true scale.
			# lon_0,lat_0 is central point.
			p = Proj(proj=proj, lat_0=lat_0, lon_0=lon_0, lat_ts=lat_ts, ellps=ellps)
			proj4_args = '+proj=' + str(proj) + ' ' +              '+lat_0=' + str(lat_0) + ' ' +              '+lon_0=' + str(lon_0) + ' ' +              '+lat_ts=' + str(lat_ts) + ' ' +              '+ellps=' + str(ellps)
		elif lon_0!=None and lat_0!=None and lat_ts==None:
			# lon_0,lat_0 is central point.
			p = Proj(proj=proj, lat_0=lat_0, lon_0=lon_0, ellps=ellps)
			proj4_args = '+proj=' + str(proj) + ' ' +              '+lat_0=' + str(lat_0) + ' ' +              '+lon_0=' + str(lon_0) + ' ' +              '+ellps=' + str(ellps)
		elif lon_0==None and lat_0==None and lat_ts==None:
			# lon_0,lat_0 is central point.
			lat_0 = (min(latlim) + max(latlim)) / 2
			lon_0 = (min(lonlim) + max(lonlim)) / 2
			p = Proj(proj=proj, lat_0=lat_0, lon_0=lon_0, ellps=ellps)
			proj4_args = '+proj=' + str(proj) + ' ' +              '+lat_0=' + str(lat_0) + ' ' +              '+lon_0=' + str(lon_0) + ' ' +              '+ellps=' + str(ellps)

		# area_extent defined as (x_min, y_min, x_max, y_max)
		left_ex1, up_ex1 = p(left, up)
		right_ex1, up_ex2 = p(right, up)
		left_ex2, down_ex1 = p(left, down)
		right_ex2, down_ex2 = p(right, down)

		area_extent = (min(left_ex1, left_ex2),
					   min(up_ex1, up_ex2),
					   max(right_ex1, right_ex2),
					   max(down_ex1, down_ex2))

		# минимум из всех координат X, Y, максимум из всех координат X, Y
		# Такой результат даёт правильный area_extent для 3413
		# При этом для 4326 area_extent остаётся неизменным
		# area_def_3413 = swath_area_def(name='Temporal SWATH EPSG Projection 3413', proj='stere', \
		#                                lonlim=(-180,180), latlim=(30,90), ellps="WGS84", res=1500, \
		#                                lat_ts=70, lat_0=90, lon_0=-45)
		# Area extent: (-5050747.263141337, 0.0, 0.0, 5050747.263141336)
		area_extent = (min(left_ex1, left_ex2, right_ex1, right_ex2),
					   min(up_ex1, up_ex2, down_ex1, down_ex2),
					   max(left_ex1, left_ex2, right_ex1, right_ex2),
					   max(up_ex1, up_ex2, down_ex1, down_ex2))

		# Using abs() to avoid negative numbers of coloumns/rows as for epsg3413 for example
		xsize = abs(int((area_extent[2] - area_extent[0]) / res[0]))
		ysize = abs(int((area_extent[3] - area_extent[1]) / res[1]))

		swath_area_def = pr.utils.get_area_def(area_id, name, proj_id, proj4_args, xsize, ysize, area_extent)

		#~ print swath_area_def

		return swath_area_def



	# In[39]:

	scale = 5

	sigma0w = {}
	roughness = {}

	print "Scale set to: \'%s\' " %scale

	for p in polarization:
		print "Filtering Image: \'%s\' polarization" %p
		
		# filter the image
		sigma0w[p] = wiener(sigma0[p][::scale,::scale], mysize=(7,7), noise=None)
	#     sigma0w[p] = sigma0[p]

	# Close ZIP-file
	zf.close()


	# In[40]:

	# S1 Pixel resolution
	# we use pxlResSAR for further GSHHS rasterizing and reprojecting data with pyresample

	lonlim = (lons_2[::scale,::scale].min(),lons_2[::scale,::scale].max())
	latlim = (lats_2[::scale,::scale].min(),lats_2[::scale,::scale].max())

	# enlarge lonlims for cropping a bit larger area for masking
	lonlimGSHHS = (lonlim[0]-2, lonlim[1]+2)
	latlimGSHHS = (latlim[0]-2, latlim[1]+2)


	# Get first guess pixel resolution
	import distancelib
	pxlResSARm  = asarray(distancelib.getPixelResolution(lats_2[::scale,::scale], lons_2[::scale,::scale], lons_2[::scale,::scale].shape, 'km'))*1e3
	pxlResSARdeg  = asarray(distancelib.getPixelResolution(lats_2[::scale,::scale],   lons_2[::scale,::scale],   lons_2[::scale,::scale].shape, 'deg'))

	print "S1 cell resolution, %s deg"  % str(pxlResSARdeg)
	print "S1 cell resolution, %s m"  % str(pxlResSARm)


	# In[41]:

	import pyresample as pr
	from pyproj import Proj

	# Define areas with pyresample
	swath_def = pr.geometry.SwathDefinition(lons=lons_2[::scale,::scale], lats=lats_2[::scale,::scale])

	area_def_4326 = swath_area_def(name='Temporal SWATH EPSG Projection 4326', proj='eqc',
							  lonlim=lonlimGSHHS, latlim=latlimGSHHS, ellps="WGS84", res=pxlResSARm)


	# In[42]:

	# Get the SAR pixel resolution from the area_def for further identical shapes
	up    = min(latlimGSHHS)
	down  = max(latlimGSHHS)
	left  = min(lonlimGSHHS)
	right = max(lonlimGSHHS)
	area_extent_deg = (left, down, right, up)

	area_extent_deg_shape = area_def_4326.shape

	pxlResSARdeg = asarray( (abs(area_extent_deg[2] - area_extent_deg[0]) / float(area_extent_deg_shape[1]), abs(area_extent_deg[3] - area_extent_deg[1]) / float(area_extent_deg_shape[0])) )

	pxlResSARm = asarray( (area_def_4326.pixel_size_x, area_def_4326.pixel_size_y) )
	print "S1 cell resolution, %s deg"  % str(pxlResSARdeg)
	print "S1 cell resolution, %s m"  % str(pxlResSARm)


	# In[43]:

	# Apply Mask from GSHHS

	import gshhs_rasterize
	reload(gshhs_rasterize)

	# ESRI shapefile containing land polygons
	shapefile = '/media/SOLabNFS/store/auxdata/coastline/GSHHS_shp/f/GSHHS_f_L1.shp'

	# reproject GSHHS onto S1 grid before calculations
	print "Rasterizing Land Mask"
	mask_arr_4326 = gshhs_rasterize.gshhs_rasterize_4326(lonlimGSHHS, latlimGSHHS, pxlResSARdeg, area_def_4326.shape, True, shapefile)


	# In[44]:

	mask_arr_swath = pr.kd_tree.resample_nearest(area_def_4326, mask_arr_4326.ravel(), swath_def, radius_of_influence=4*pxlResSARm.max(), epsilon=0.5, fill_value=None, nprocs=numProcs)


	# In[ ]:

	# print mask_arr_4326.shape, mask_arr_swath.shape, sigma0w[p].shape


	# In[45]:

	# Nice Image (Roughness)

	sigma0wAvg = {}
	roughnessNrmlzd = {}

	for p in polarization:
		print "Nice Image: \'%s\' polarization" %p

		roughness[p] = ma.masked_where(mask_arr_swath, sigma0w[p])

		sigma0wAvg[p] = ma.median(roughness[p], axis=0)

		roughnessNrmlzd[p] = (roughness[p]-sigma0wAvg[p])/sigma0wAvg[p]


	# In[47]:

	# Reproject

	#~ roughness_4326 = {}

	#~ for p in polarization:
		#~ print "Reprojecting Image: \'%s\' polarization" %p
#~ 
		#~ roughness_4326[p] = pr.kd_tree.resample_nearest(swath_def, roughnessNrmlzd[p].ravel(), area_def_4326, radius_of_influence=4*pxlResSARm.max(), epsilon=0.5, fill_value=None, nprocs=numProcs)


	# In[53]:

	import os
	# Save images
	oPath = '/home/mag/tmp/'

	for p in polarization:
		print "Saving Image: \'%s\' polarization" %p
		plt.close('all')
		oFileName = os.path.join(oPath, fn[:-3]+p+'_bone_r.png')
		plt.imsave(oFileName, roughnessNrmlzd[p], vmin=-1, vmax=1, cmap=plt.cm.bone_r)
		oFileName = os.path.join(oPath, fn[:-3]+p+'_RdBu_r.png')
		plt.imsave(oFileName, roughnessNrmlzd[p], vmin=-1, vmax=1, cmap=plt.cm.RdBu_r)
		#~ oFileName = os.path.join(oPath, fn[:-3]+p+'_4326_bone_r.png')
		#~ plt.imsave(oFileName, roughness_4326[p], vmin=-1, vmax=1, cmap=plt.cm.bone_r)
		#~ oFileName = os.path.join(oPath, fn[:-3]+p+'_4326_RdBu_r.png')
		#~ plt.imsave(oFileName, roughness_4326[p], vmin=-1, vmax=1, cmap=plt.cm.RdBu_r)
       
