
# coding: utf-8

# In[1]:

# In[2]:

# In[3]:

from numpy import asarray, zeros, reshape, double, arange, ma, log10, diff, mean, flipud, floor, pi, sqrt, size
import pyresample as pr
from pyproj import Proj

# In[4]:

from IPython.html import widgets
# [widget for widget in dir(widgets) if widget.endswith('Widget')]


# In[5]:

import zipfile

import datetime

import matplotlib.pyplot as plt

from matplotlib.mlab import find
import xmltodict

from PIL import Image
import StringIO

from scipy.interpolate import RectBivariateSpline
from scipy.signal import wiener

import os

import distancelib

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2014, 10, 28)
__modified__ = datetime.datetime(2014, 5, 27)
__version__  = "1.0"
__status__   = "Development"


# In[6]:

from multiprocessing import cpu_count
numProcs = cpu_count()


# In[7]:

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


# In[8]:

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


# In[9]:

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
        try:
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
        except Exception:
            ncepFileName = 'gfs.' + startTime.strftime("%Y%m%d") + '%.2d' %(baseHour) + '/gfs.t' + '%.2d' %(baseHour) + 'z.pgrb2.0p25.f' + '%.3d' %(forecastHour)

            grbs = pygrib.open(iPath_wind + ncepFileName)

            u_wind = grbs.message(1)
            v_wind = grbs.message(2)
            u = u_wind['values']
            v = v_wind['values']
            lats_wind = u_wind['latitudes']
            lons_wind = u_wind['longitudes']
            lons_wind = reshape(lons_wind, (lons_wind.shape[0]/1440, 1440))
            lats_wind = reshape(lats_wind, (lats_wind.shape[0]/1440, 1440))

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
    ncepGFSmodel['u'] = u
    ncepGFSmodel['v'] = v
    ncepGFSmodel['baseHour'] = baseHour
    ncepGFSmodel['forecastHour'] = forecastHour
#     del u_wind, v_wind
    return ncepGFSmodel


# In[10]:

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

    up    = max(latlim)
    down  = min(latlim)
    left  = min(lonlim)
    right = max(lonlim)
    
    print 'up, down, left, right: ', round(up), round(down), round(left), round(right)

    area_id = name.replace(" ", "_").lower()
    proj_id = area_id

    if proj == 'eqc':
        p = Proj(proj=proj, llcrnrlat=up, urcrnrlat=down, llcrnrlon=left, urcrnrlon=right, ellps=ellps)
        proj4_args = '+proj=' + str(proj) + ' ' + \
             '+llcrnrlat=' + str(up) + ' ' + \
             '+urcrnrlat=' + str(down) + ' ' + \
             '+llcrnrlon=' + str(left) + ' ' + \
             '+urcrnrlon=' + str(right) + ' ' + \
             '+ellps=' + str(ellps)
    elif lat_ts!=None and lat_0!=None:
        # lat_ts is latitude of true scale.
        # lon_0,lat_0 is central point.
        p = Proj(proj=proj, lat_0=lat_0, lon_0=lon_0, lat_ts=lat_ts, ellps=ellps)
        proj4_args = '+proj=' + str(proj) + ' ' + \
             '+lat_0=' + str(lat_0) + ' ' + \
             '+lon_0=' + str(lon_0) + ' ' + \
             '+lat_ts=' + str(lat_ts) + ' ' + \
             '+ellps=' + str(ellps)
    elif lon_0!=None and lat_0!=None and lat_ts==None:
        # lon_0,lat_0 is central point.
        p = Proj(proj=proj, lat_0=lat_0, lon_0=lon_0, ellps=ellps)
        proj4_args = '+proj=' + str(proj) + ' ' + \
             '+lat_0=' + str(lat_0) + ' ' + \
             '+lon_0=' + str(lon_0) + ' ' + \
             '+ellps=' + str(ellps)
    elif lon_0==None and lat_0==None and lat_ts==None:
        # lon_0,lat_0 is central point.
        lat_0 = (up + down) / 2
        lon_0 = (right + left) / 2
        p = Proj(proj=proj, lat_0=lat_0, lon_0=lon_0, ellps=ellps)
        proj4_args = '+proj=' + str(proj) + ' ' + \
             '+lat_0=' + str(lat_0) + ' ' + \
             '+lon_0=' + str(lon_0) + ' ' + \
             '+ellps=' + str(ellps)

    left_ex1, up_ex1 = p(left, up)
    right_ex1, up_ex2 = p(right, up)
    left_ex2, down_ex1 = p(left, down)
    right_ex2, down_ex2 = p(right, down)

    if proj == 'stere':
        lon = (left+right)/2.0
        if (lon >=0 and lon <90) or (lon >=-360 and lon < -270):
            print 11111111111
            area_extent = (
                           min(left_ex1, left_ex2, right_ex1, right_ex2),
                           min(down_ex1, down_ex2, up_ex1, up_ex2),
                           max(left_ex1, left_ex2, right_ex1, right_ex2),
                           max(down_ex1, down_ex2, up_ex1, up_ex2)
                        )
        elif (lon >=90 and lon <180) or (lon >=-270 and lon < -180):
            print 2222222222222
            area_extent = (
                           max(left_ex1, left_ex2, right_ex1, right_ex2),
                           max(down_ex1, down_ex2, up_ex1, up_ex2),
                           min(left_ex1, left_ex2, right_ex1, right_ex2),
                           min(down_ex1, down_ex2, up_ex1, up_ex2)
                        )
        elif (lon >= 180 and lon < 270) or (lon >= -180 and lon < -90):
            print 333333333333
            area_extent = (
                           min(left_ex1, left_ex2, right_ex1, right_ex2),
                           min(down_ex1, down_ex2, up_ex1, up_ex2),
                           max(left_ex1, left_ex2, right_ex1, right_ex2),
                           max(down_ex1, down_ex2, up_ex1, up_ex2)
                        )
        else:
            print 44444444444444444
            area_extent = (
                           min(left_ex1, left_ex2, right_ex1, right_ex2),
                           min(down_ex1, down_ex2, up_ex1, up_ex2),
                           max(left_ex1, left_ex2, right_ex1, right_ex2),
                           max(down_ex1, down_ex2, up_ex1, up_ex2)
                        )
    else:
        # минимум из всех координат X, Y, максимум из всех координат X, Y
        # Такой результат даёт правильный area_extent для 3413
        # При этом для 4326 area_extent остаётся неизменным
        # area_def_3413 = swath_area_def(name='Temporal SWATH EPSG Projection 3413', proj='stere', \
        #                                lonlim=(-180,180), latlim=(30,90), ellps="WGS84", res=1500, \
        #                                lat_ts=70, lat_0=90, lon_0=-45)
        # Area extent: (-5050747.263141337, 0.0, 0.0, 5050747.263141336)
        area_extent = (
                        min(left_ex1, left_ex2, right_ex1, right_ex2),
                        min(up_ex1, up_ex2, down_ex1, down_ex2),
                        max(left_ex1, left_ex2, right_ex1, right_ex2),
                        max(up_ex1, up_ex2, down_ex1, down_ex2)
                    )
    
    #~ print 'left: ', left_ex1, left_ex2
    #~ print 'right: ', right_ex1, right_ex2
    #~ print 'up: ', up_ex1, up_ex2
    #~ print 'down: ', down_ex1, down_ex2

#     Using abs() to avoid negative numbers of coloumns/rows as for epsg3413 for example
    xsize = abs(int((area_extent[2] - area_extent[0]) / res[0]))
    ysize = abs(int((area_extent[3] - area_extent[1]) / res[1]))
    
    swath_area_def = pr.utils.get_area_def(area_id, name, proj_id, proj4_args, xsize, ysize, area_extent)

#     print swath_area_def

    return swath_area_def




# In[11]:

# READ THE RAW_COUNTS from GRD image

def read_raw_counts(fn, fileLocation, polarization, scale):


    # Note that For SLC images BitsPerSample=32 and for GRD BitsPerSample=16

    # im = Image.open(inpath + fileLocation['s1aiwgrd' + polarization])

    im = zf.read(fn[:-4] + '.SAFE' + fileLocation[fn.lower().replace("_","")[0:8] + polarization][1:])
    im = StringIO.StringIO(im) #Encode the raw data to be used by Image.open()
    im = Image.open(im)        #Open the image

    return asarray(im)[::scale,::scale]


# In[12]:

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


# In[13]:

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


# In[14]:

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

def readS1_raw_counts(inpath = '/media/SOLabNFS2/tmp/sentinel-1/Svalbard-Barents/', \
fn='S1A_EW_GRDH_1SDH_20141003T133957_20141003T134057_002666_002F83_0BE7.zip', resolution=None):
    """Reading raw_counts from Sentinel-1 images"""

    global zf
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


    # In[17]:

    # Get the productType/polarization

    metadataObject = manifest['xfdu:XFDU']['metadataSection']['metadataObject']
    for n in range(len(metadataObject)):
        if metadataObject[n]['@ID'] == 'generalProductInformation':
            productType = metadataObject[n]['metadataWrap']['xmlData']['s1sarl1:standAloneProductInformation']['s1sarl1:productType']
            transmitterReceiverPolarisation = metadataObject[n]['metadataWrap']['xmlData']['s1sarl1:standAloneProductInformation']['s1sarl1:transmitterReceiverPolarisation']

    polarization = []
    for p in range(0,len(transmitterReceiverPolarisation[0].lower())):
        if len(transmitterReceiverPolarisation[0].lower()) == 1:
            polarization = transmitterReceiverPolarisation.lower()
        else:
            polarization.append(transmitterReceiverPolarisation[p].lower())
    if isinstance(polarization, basestring):
        polarization = [polarization]
    print "Available polarizations: \'%s\'" %polarization

    # In[18]:

    # %%timeit -n 1 -r 1

    # set default resolution
    if resolution is None:
        resolution = 80

    raw_counts = {}

    GEOgrid = read_anotation(fn, fileLocation, polarization[0])

    # Find scale to reduce image to the specified resolution
    arrShape =  (GEOgrid['line'].max()+1, GEOgrid['pixel'].max()+1)
    scale = resolution/round(mean(asarray(distancelib.getPixelResolution(GEOgrid['lats'], \
                                                                         GEOgrid['lons'], \
                                                                         arrShape, 'km'))*1e3))
    for p in polarization:
        print "Reading raw_counts: \'%s\' polarization" %p
        # READ THE RAW_COUNTS from GRD image
        raw_counts[p] = read_raw_counts(fn, fileLocation, p, scale)

    # Close ZIP-file
    zf.close()

    return raw_counts, polarization, manifest, GEOgrid

class readS1:
    """\
    Initial class for Reading and Interpolating data from Sentinel-1 images
    Usage:
        s1 = readS1(inpath, fn)
            pn = pathname
            fn = filename
    """
    def __init__(self, inpath = '/media/SOLabNFS2/tmp/sentinel-1/Svalbard-Barents/', \
    fn='S1A_EW_GRDH_1SDH_20141003T133957_20141003T134057_002666_002F83_0BE7.zip', resolution = None):
        """Reading and Interpolating data from Sentinel-1 images"""

        # READ THE MANIFEST - top level

        # manifest = open(inpath + "manifest.safe", "r") # Open a file in read-only mode
        # manifest = manifest.read()      # read the file object

        # open zip file without extraction
        # finnGulf
        # fn = 'S1A_IW_GRDH_1SDV_20141004T155619_20141004T155644_002682_002FE5_BE58.zip'

        # Svalbard-Barents
        #~ inpath = '/media/SOLabNFS2/tmp/sentinel-1/Svalbard-Barents/'
        #~ fileNameList = ['S1A_EW_GRDH_1SDH_20141003T071321_20141003T071421_002662_002F6B_AFCF.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141003T071221_20141003T071321_002662_002F6B_0F56.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141003T071117_20141003T071221_002662_002F6B_CAFE.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141004T061601_20141004T061701_002676_002FBD_3300.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141004T061701_20141004T061800_002676_002FBD_EF64.zip',
                         #~ 'S1A_EW_GRDM_1SDH_20141004T061601_20141004T061701_002676_002FBD_B433.zip',
                         #~ 'S1A_EW_GRDM_1SDH_20141004T061357_20141004T061501_002676_002FBD_A7E4.zip',
                         #~ 'S1A_EW_GRDM_1SDH_20141004T061701_20141004T061800_002676_002FBD_89ED.zip',
                         #~ 'S1A_EW_GRDM_1SDH_20141004T061501_20141004T061601_002676_002FBD_A842.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141004T061501_20141004T061601_002676_002FBD_4B66.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141004T061357_20141004T061501_002676_002FBD_D530.zip',
                         #~ 'S1A_EW_GRDM_1SDH_20141003T151915_20141003T152019_002667_002F88_B77B.zip',
                         #~ 'S1A_EW_GRDM_1SDH_20141003T152119_20141003T152229_002667_002F88_1BCA.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141003T152119_20141003T152229_002667_002F88_8003.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141003T152019_20141003T152119_002667_002F88_9106.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141003T151915_20141003T152019_002667_002F88_C200.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141003T134257_20141003T134338_002666_002F83_5FD6.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141003T134057_20141003T134157_002666_002F83_7D9C.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141003T134157_20141003T134257_002666_002F83_30DE.zip',
                         #~ 'S1A_EW_GRDH_1SDH_20141003T133957_20141003T134057_002666_002F83_0BE7.zip']

        # inpath = '/media/SOLabNFS2/tmp/sentinel-1/'
        # fileNameList = ['S1A_IW_GRDH_1SDV_20141004T155619_20141004T155644_002682_002FE5_BE58.zip']


        #~ fn = fileNameList[-1]
        global zf
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


        # In[17]:

        # Get the productType/polarization

        metadataObject = manifest['xfdu:XFDU']['metadataSection']['metadataObject']
        for n in range(len(metadataObject)):
            if metadataObject[n]['@ID'] == 'generalProductInformation':
                productType = metadataObject[n]['metadataWrap']['xmlData']['s1sarl1:standAloneProductInformation']['s1sarl1:productType']
                transmitterReceiverPolarisation = metadataObject[n]['metadataWrap']['xmlData']['s1sarl1:standAloneProductInformation']['s1sarl1:transmitterReceiverPolarisation']

        polarization = []
        for p in range(0,len(transmitterReceiverPolarisation[0].lower())):
            if len(transmitterReceiverPolarisation[0].lower()) == 1:
                polarization = transmitterReceiverPolarisation.lower()
            else:
                polarization.append(transmitterReceiverPolarisation[p].lower())
        # If only one polarization, it is represented as basestring, we make it a list, to avoid mistakes in loops
        if isinstance(polarization, basestring):
            polarization = [polarization]
        print "Available polarizations: \'%s\'" %polarization

        # In[18]:

        # %%timeit -n 1 -r 1

        # set default resolution
        if resolution is None:
            resolution = 80

        raw_counts = {}
        lats_2 = {}
        lons_2 = {}
        incidenceAngle_2 = {}
        sigmaNought_2 = {}
        noiseLut_2 = {}
        sigma0 = {}

        # NB! if len(polarization[0]) == 1 then there is only one polarization, meaning polarization=='vv'
        # READ the ANNOTATION
        GEOgrid = read_anotation(fn, fileLocation, polarization[0])

        # READ the CALIBRATION LUTs
        cLUTs = read_clbrtn_luts(fn, fileLocation, polarization[0])

        # READ the NOISE LUTs
        nLUTs = read_noise_luts(fn, fileLocation, polarization[0])

        # Find scale to reduce image to the specified resolution
        arrShape =  (GEOgrid['line'].max()+1, GEOgrid['pixel'].max()+1)
        scale = resolution/round(mean(asarray(distancelib.getPixelResolution(GEOgrid['lats'], \
                                                                             GEOgrid['lons'], \
                                                                             arrShape, 'km'))*1e3))

        for p in polarization:
            print "Reading raw_counts: \'%s\' polarization" %p
            # READ THE RAW_COUNTS from GRD image
            raw_counts[p] = read_raw_counts(fn, fileLocation, p, scale)

        #  ----------------------------------
        # INTERPOLATE DATA
        # Interpolate Geolocation grid points and calibration LUTs onto grid of raw_counts.shape

        # Serial Processing

        # Serial loop is faster than the Parallel Loop (see appropriate ipnb file)
        # 1 loops, best of 1: 31.5 s per loop
        # %%timeit -n 1 -r 1

        # create new grid of raw_counts shape
        #~ line_2 = arange(raw_counts[p].shape[0])
        #~ pixel_2 = arange(raw_counts[p].shape[1])
        line_2 = arange(arrShape[0])[::scale]
        pixel_2 = arange(arrShape[1])[::scale]

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
            # for VH, HV - multiply S1 noiseLUTs by nLtCoeff=1e10
            # for VV, HH - multiply S1 noiseLUTs by nLtCoeff=sqrt(2)*1e10
            if p=='vv' or p=='hh':
                nLtCoeff=sqrt(2)*1e10
            elif p=='vh' or p=='hv':
                nLtCoeff=1e10

            sigma0[p] = ( double(raw_counts[p])**2 - noiseLut_2[p]*nLtCoeff )/sigmaNought_2[p]**2

        #~ Putting others vars to self
        self.polarization = polarization
        self.raw_counts = raw_counts
        self.incidenceAngle_2 = incidenceAngle_2
        self.sigmaNought_2 = sigmaNought_2
        self.noiseLut_2 = noiseLut_2
        self.sigma0 = sigma0
        self.lons_2 = lons_2
        self.lats_2 = lats_2
        self.GEOgrid =GEOgrid
        self.cLUTs = cLUTs
        self.nLUTs = nLUTs
        self.manifest = manifest

        # Close ZIP-file
        zf.close()
    #~ return incidenceAngle_2, sigmaNought_2, noiseLut_2, sigma0, lons_2, lats_2
