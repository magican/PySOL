# coding: utf-8

# In[1]:


# In[2]:

from numpy import asarray, zeros, reshape, double, arange, ma, log10, diff, mean, flipud, floor, pi, sqrt, size, fliplr, meshgrid, exp, cos, radians

# In[3]:


# In[4]:

import datetime

import matplotlib.pyplot as plt
import matplotlib
#~ matplotlib.use('Agg')
from pylab import *

from matplotlib.mlab import find

from PIL import Image

from scipy.signal import wiener

import os

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2014, 10, 28)
__modified__ = datetime.datetime(2015, 6, 19)
__version__  = "1.0"
__status__   = "Development"

def mkdirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)

# Svalbard-Barents
inpath = '/media/SOLabNFS2/tmp/different_SAR/sentinel-1/Svalbard-Barents/'
fileNameList = ['S1A_EW_GRDH_1SDH_20141003T071321_20141003T071421_002662_002F6B_AFCF.zip',
                'S1A_EW_GRDH_1SDH_20141003T071221_20141003T071321_002662_002F6B_0F56.zip',
                'S1A_EW_GRDH_1SDH_20141003T071117_20141003T071221_002662_002F6B_CAFE.zip',
                'S1A_EW_GRDH_1SDH_20141004T061601_20141004T061701_002676_002FBD_3300.zip',
                'S1A_EW_GRDH_1SDH_20141004T061701_20141004T061800_002676_002FBD_EF64.zip',
                'S1A_EW_GRDM_1SDH_20141004T061601_20141004T061701_002676_002FBD_B433.zip',
                'S1A_EW_GRDM_1SDH_20141004T061357_20141004T061501_002676_002FBD_A7E4.zip',
                'S1A_EW_GRDM_1SDH_20141004T061701_20141004T061800_002676_002FBD_89ED.zip',
                'S1A_EW_GRDM_1SDH_20141004T061501_20141004T061601_002676_002FBD_A842.zip',
                'S1A_EW_GRDH_1SDH_20141004T061501_20141004T061601_002676_002FBD_4B66.zip',
                'S1A_EW_GRDH_1SDH_20141004T061357_20141004T061501_002676_002FBD_D530.zip',
                'S1A_EW_GRDM_1SDH_20141003T151915_20141003T152019_002667_002F88_B77B.zip',
                'S1A_EW_GRDM_1SDH_20141003T152119_20141003T152229_002667_002F88_1BCA.zip',
                'S1A_EW_GRDH_1SDH_20141003T152119_20141003T152229_002667_002F88_8003.zip',
                'S1A_EW_GRDH_1SDH_20141003T152019_20141003T152119_002667_002F88_9106.zip',
                'S1A_EW_GRDH_1SDH_20141003T151915_20141003T152019_002667_002F88_C200.zip',
                'S1A_EW_GRDH_1SDH_20141003T134257_20141003T134338_002666_002F83_5FD6.zip',
                'S1A_EW_GRDH_1SDH_20141003T134057_20141003T134157_002666_002F83_7D9C.zip',
                'S1A_EW_GRDH_1SDH_20141003T134157_20141003T134257_002666_002F83_30DE.zip',
                'S1A_EW_GRDH_1SDH_20141003T133957_20141003T134057_002666_002F83_0BE7.zip',
                'S1A_IW_GRDH_1SDV_20150506T045446_20150506T045511_005796_00772F_623B.zip',
                'S1A_IW_GRDH_1SDV_20150506T161536_20150506T161601_005803_00775C_25DB.zip']

# Liza Polar Low 2015 June 01
#inpath = '/media/SOLabNFS2/tmp/different_SAR/sentinel-1/Liza_PL_01_June_2015/'
#fileNameList = ['S1A_EW_GRDM_1SDH_20150601T173700_20150601T173804_006183_0080D9_94AD.zip',
#                'S1A_EW_GRDM_1SDH_20150601T173804_20150601T173858_006183_0080D9_E1A6.zip',
#                'S1A_IW_GRDH_1SDV_20150601T173611_20150601T173644_006183_0080D8_B3F1.zip']

# Ania_Ladoga_29_May_2015/
#inpath = '/media/SOLabNFS2/tmp/different_SAR/sentinel-1/Ania_Ladoga_29_May_2015/'
#fileNameList = ['S1A_IW_GRDH_1SDV_20150603T154002_20150603T154027_006211_0081A9_5F10.zip',
#                'S1A_IW_GRDH_1SDV_20150529T041657_20150529T041722_006131_007F51_F751.zip',
#                'S1A_EW_GRDM_1SDH_20150517T153117_20150517T153221_005963_007AED_56B0.zip']

fn = fileNameList[4]

import re
prog = re.compile(r'(\d{8})')
file_date = prog.findall(fn)[0]

year = file_date[:4]
month = file_date[4:6]
day = file_date[6:]

# In[6]:

import readS1
from readS1 import *

# In[7]:

s1 = readS1(inpath=inpath, fn=fn)
# s1.__dict__['raw_counts']


# In[8]:

# get vars from s1 class
for k, v in s1.__dict__.iteritems():
    locals()[k]=v

# In[9]:

scale = 4

sigma0w = {}
roughness = {}

print "Scale set to: \'%s\' " %scale

for p in polarization:
    print "Filtering Image: \'%s\' polarization" %p
    
    # filter the image
    sigma0w[p] = wiener(sigma0[p][::scale,::scale], mysize=(7,7), noise=None)
#     sigma0w[p] = sigma0[p]


# In[10]:

# S1 Pixel resolution
# we use pxlResSAR for further GSHHS rasterizing and reprojecting data with pyresample

lonlim = (lons_2[::scale,::scale].min(),lons_2[::scale,::scale].max())
latlim = (lats_2[::scale,::scale].min(),lats_2[::scale,::scale].max())

# enlarge lonlims for cropping a bit larger area for masking
lonlimGSHHS = (lonlim[0]-1.0, lonlim[1]+1.0)
latlimGSHHS = (latlim[0]-1.0, latlim[1]+1.0)


# Get first guess pixel resolution
import distancelib
pxlResSARm  = asarray(distancelib.getPixelResolution(lats_2[::scale,::scale], lons_2[::scale,::scale], lons_2[::scale,::scale].shape, 'km'))*1e3
pxlResSARdeg  = asarray(distancelib.getPixelResolution(lats_2[::scale,::scale], lons_2[::scale,::scale],   lons_2[::scale,::scale].shape, 'deg'))


# In[11]:

import pyresample as pr
from pyproj import Proj

# Define areas with pyresample
swath_def = pr.geometry.SwathDefinition(lons=lons_2[::scale,::scale], lats=lats_2[::scale,::scale])

area_def_4326 = swath_area_def(name='Temporal SWATH EPSG Projection 4326', proj='eqc',
                          lonlim=lonlimGSHHS, latlim=latlimGSHHS, ellps="WGS84", res=pxlResSARm)


# In[12]:

# Get the SAR pixel resolution from the area_def for further identical shapes
up    = min(latlimGSHHS)
down  = max(latlimGSHHS)
left  = min(lonlimGSHHS)
right = max(lonlimGSHHS)
area_extent_deg = (left, down, right, up)

area_extent_deg_shape = area_def_4326.shape

pxlResSARdeg = asarray( (abs(area_extent_deg[2] - area_extent_deg[0]) / float(area_extent_deg_shape[1]),                 abs(area_extent_deg[3] - area_extent_deg[1]) / float(area_extent_deg_shape[0])) )

pxlResSARm = asarray( (area_def_4326.pixel_size_x, area_def_4326.pixel_size_y) )
print "S1 cell resolution, %s deg"  % str(pxlResSARdeg)
print "S1 cell resolution, %s m"  % str(pxlResSARm)


# In[13]:

# Apply Mask from GSHHS

import gshhs_rasterize
reload(gshhs_rasterize)

# ESRI shapefile containing land polygons
shapefile = '/media/SOLabNFS/store/auxdata/coastline/GSHHS_shp/f/GSHHS_f_L1.shp'

# reproject GSHHS onto S1 grid before calculations
print "Rasterizing Land Mask"
mask_arr_4326 = gshhs_rasterize.gshhs_rasterize_4326(lonlimGSHHS, latlimGSHHS, pxlResSARdeg, area_def_4326.shape, True,                                      shapefile)


# In[14]:

mask_arr_swath = pr.kd_tree.resample_nearest(area_def_4326, mask_arr_4326, swath_def, radius_of_influence=4*pxlResSARm.max(), epsilon=0.5, fill_value=None)


# In[15]:

print area_def_4326.shape, mask_arr_4326.shape
print mask_arr_swath.shape, sigma0w[p].shape, swath_def.shape


# In[16]:

# Nice Image (Roughness)

sigma0wAvg = {}
roughnessNrmlzd = {}

if len(polarization[0])>=2: # if 2 polarizations
    for p in polarization:
        print "Nice Image: \'%s\' polarization" %p
        roughness[p] = ma.masked_where(mask_arr_swath, sigma0w[p])
        sigma0wAvg[p] = ma.median(roughness[p], axis=0)
        roughnessNrmlzd[p] = (roughness[p]-sigma0wAvg[p])/sigma0wAvg[p]
elif len(polarization[0])==1: # if only 1 polarization
    p = polarization
    print "Nice Image: \'%s\' polarization" %p
    roughness[p] = ma.masked_where(mask_arr_swath, sigma0w[p])
    sigma0wAvg[p] = ma.median(roughness[p], axis=0)
    roughnessNrmlzd[p] = (roughness[p]-sigma0wAvg[p])/sigma0wAvg[p]  

del roughness, sigma0wAvg

# In[ ]:




# # Adding Model wind

# In[17]:

# Adding Model wind

# import xmltodict

# zf = zipfile.ZipFile(inpath+fn, 'r')
# manifest = zf.read(fn[:-4] + '.SAFE/manifest.safe')
# manifest = xmltodict.parse(manifest) # Parse the read document string
# zf.close()

startTime = datetime.datetime.strptime(manifest['xfdu:XFDU']['metadataSection']['metadataObject'][12]['metadataWrap']['xmlData']['safe:acquisitionPeriod']['safe:startTime'],                              "%Y-%m-%dT%H:%M:%S.%f")

ncepGFSmodelWind = ncepGFSmodel(startTime, lats_2, lons_2)


# In[18]:

# Reprojecting data

import distancelib

# Pixel resolution
# we use pxlResWind/pxlResSAR for further pyresample radius_of_influence and sigmas
pxlResWind = asarray(distancelib.getPixelResolution(ncepGFSmodelWind['lats_wind'],                                                     ncepGFSmodelWind['lons_wind'],                                                     ncepGFSmodelWind['lons_wind'].shape, 'km'))
# pxlResSAR  = asarray(distancelib.getPixelResolution(lats_2, lons_2, lons_2.shape, 'km'))*1e3

# Note pxlResWind is in KM, multiply by 1e3 for meters
print "S1 cell resolution, %s m"  % pxlResSARm
print "Wind cell resolution, %s km" % pxlResWind


# In[19]:

from scipy.interpolate import RectSphereBivariateSpline

def ncepGFSmodel2swath(lats, lons, data, lats_2, lons_2):

    func = RectSphereBivariateSpline(lats, lons, data)
    data_2 = func.ev(lats_2.ravel()*pi/180,                      lons_2.ravel()*pi/180)                     .reshape(lats_2.shape)
    return data_2


# In[20]:

# reproject NCEP onto S1 grid before calculations
# Using RectSphereBivariateSpline - Bivariate spline approximation over a rectangular mesh on a sphere
# as it is much more efficiant for full resolution
# as well as smoothes nicely the image

# We don't want to work with full res wind so scaling the image for about 100m resolution
# Adjust scale to get appropriate value

lts = flipud(ncepGFSmodelWind['lats_wind'])[:,0]*pi/180
lns = ncepGFSmodelWind['lons_wind'][0,:]*pi/180

lts_2 = lats_2[::scale,::scale]
lns_2 = lons_2[::scale,::scale]

ncepGFSmodelWindSwath = {}
ncepGFSmodelWindSwath['wind_speed'] = ncepGFSmodel2swath(lts, lns, flipud(ncepGFSmodelWind['wind_speed']), lts_2, lns_2)
ncepGFSmodelWindSwath['wind_dir']   = ncepGFSmodel2swath(lts, lns, flipud(ncepGFSmodelWind['wind_dir']),   lts_2, lns_2)
ncepGFSmodelWindSwath['u']   = ncepGFSmodel2swath(lts, lns, flipud(ncepGFSmodelWind['u']),   lts_2, lns_2)
ncepGFSmodelWindSwath['v']   = ncepGFSmodel2swath(lts, lns, flipud(ncepGFSmodelWind['v']),   lts_2, lns_2)
    
pxlResWindSwath = asarray(distancelib.getPixelResolution(lts_2,                                                     lns_2,                                                     lns_2.shape, 'km'))

print "Interpolated Wind cell resolution, %s km" % pxlResWindSwath


# In[21]:

# calculate bearing from initial lats/lons for further wind calculation
# Taking initial values as bearing is more accurate after interpolation than vice versa
bearing = zeros((GEOgrid['lons'].shape[0]-1,GEOgrid['lons'].shape[1]))

for n in range(0,GEOgrid['lons'].shape[1]):
    col = ([GEOgrid['lats'][:-1,n], GEOgrid['lons'][:-1,n]], [GEOgrid['lats'][1:,n], GEOgrid['lons'][1:,n]])
    for m in range(0,GEOgrid['lons'].shape[0]-1):
        bearing[m][n] = distancelib.bearing(asarray(col[0])[:,m], asarray(col[1])[:,m])

# interpolate to raw_counts.shape
bearing_2 = imresize(bearing, ncepGFSmodelWindSwath['wind_dir'].shape)


# In[22]:

def PR_Mouche(theta, phi):

    A_0 = 0.00650704
    B_0 = 0.128983
    C_0 = 0.992839
    A_HALF_PI = 0.00782194
    B_HALF_PI = 0.121405
    C_HALF_PI = 0.992839
    A_PI = 0.00598416
    B_PI = 0.140952
    C_PI = 0.992885

    P_0 = A_0 * exp(B_0* theta) + C_0
    P_HALF_PI = A_HALF_PI * exp(B_HALF_PI* theta) + C_HALF_PI
    P_PI = A_PI * exp(B_PI* theta) + C_PI
    
    C0 = (P_0 + P_PI + 2 * P_HALF_PI) / 4
    C1 = (P_0 - P_PI) / 2
    C2 = (P_0 + P_PI - 2 * P_HALF_PI) / 4
    
    P = C0 + C1 * cos(radians(phi)) + C2 * cos(radians(2 * phi))
    
    return P


# In[23]:

#NB! WINDDIR = 0 WHEN WIND BLOWS TOWARDS RADAR!

p = polarization[0]

wind_dir_model_swath_rel = 90 + bearing_2 - ncepGFSmodelWindSwath['wind_dir']

if p == 'hh':
    PR = PR_Mouche(incidenceAngle_2[p][::scale,::scale], wind_dir_model_swath_rel)
    try:
        from cmod_gpu import rcs2windOpenCl
        wind_speed_asar = rcs2windOpenCl(sar=sigma0w[p]*PR,                                          windir=wind_dir_model_swath_rel,                                          theta=incidenceAngle_2[p][::scale,::scale])
    except Exception:
        from cmod_vect import rcs2windPar
        wind_speed_asar = rcs2windPar(sigma0w[p]*PR, cmdv=5,                                       windir=wind_dir_model_swath_rel,                                       theta=incidenceAngle_2[p][::scale,::scale], nprocs=numProcs)
elif p == 'vv':
    try:
        from cmod_gpu import rcs2windOpenCl
        wind_speed_asar = rcs2windOpenCl(sar=sigma0w[p],                                          windir=wind_dir_model_swath_rel,                                          theta=incidenceAngle_2[p][::scale,::scale])
    except Exception:
        from cmod_vect import rcs2windPar
        wind_speed_asar = rcs2windPar(sigma0w[p], cmdv=5,                                       windir=wind_dir_model_swath_rel,                                       theta=incidenceAngle_2[p][::scale,::scale], nprocs=numProcs)


# In[24]:

# Add mask to initial NCEP data
area_def_ncep = pr.geometry.SwathDefinition(lons=ncepGFSmodelWind['lons_wind'], lats=ncepGFSmodelWind['lats_wind'])
mask_arr_ncep = pr.kd_tree.resample_nearest(area_def_4326, mask_arr_4326, area_def_ncep,                                              radius_of_influence=4*pxlResWind.max(), epsilon=0.5, fill_value=None)
ncepGFSmodelWind['wind_speed'] = ma.masked_where(mask_arr_ncep, ncepGFSmodelWind['wind_speed'])

# Add mask to ASAR wind and reprojected NCEP
wind_speed_asar = ma.masked_where(mask_arr_swath, wind_speed_asar)
ncepGFSmodelWindSwath['wind_speed'] = ma.masked_where(mask_arr_swath, ncepGFSmodelWindSwath['wind_speed'])

lats_2, lons_2, roughnessNrmlzd, wind_speed_asar

import simplekml
def create_KML_asar(area_extent, savepath):
    kml = simplekml.Kml()

    pol = kml.newpolygon(name='area_extent', visibility=1)
    pol.tessellate = 1

    pol.altitudemode = 'clampToGround'
    pol.outerboundaryis.coords = [(min(area_extent[0], area_extent[2]),
                                   min(area_extent[1], area_extent[3])),
                                  (max(area_extent[0], area_extent[2]),
                                   max(area_extent[1], area_extent[3]))]
    if type(savepath) == list:
        for _savepath in savepath:
            kml.save(_savepath)
    else:
        kml.save(savepath)

import gdal
import os
import sys
import gc
sys.path.append('/usr/bin')
from gdal2tiles import GDAL2Tiles

def create_asar_tiles(png_filename, tiles_output_dir, proj):
    local_argv = ['/usr/bin/gdal2tiles.py', '-p', 'raster', '-r', 'cubic',
                  '-s', proj, png_filename, tiles_output_dir]
    argv = gdal.GeneralCmdLineProcessor(local_argv)
    if argv:
        gdal2tiles = GDAL2Tiles(argv[1:])
        gdal2tiles.process()

from big_image import save_big_image

import multiprocessing
numProcs = multiprocessing.cpu_count() - 2
# for proj in ['EPSG:3413']:
proj = 'EPSG:3413'

pp_names = ['roughness', 'wind_speed']
out_path = '/media/SOLabNFS2/http/tiles/Sentinel-1'

#year = '2015'
#month = '06'
#day = '03'
fileName = fn[:-4]

area_def = swath_area_def(name='Temporal SWATH EPSG Projection 3413', proj='stere',
                          lonlim=(lons_2.min(),lons_2.max()),
                          latlim=(lats_2.min(),lats_2.max()),
                          ellps="WGS84", res=pxlResSARm, lat_ts=70, lat_0=90, lon_0=-45)
# Set the parameters for GSHHS masking
proj_ = '+units=m +ellps=WGS84 +lon_0=-45 +proj=stere +lat_ts=70 +lat_0=90'
proj_name = '3413'
units = 'm'

print area_def
#print "roughness.shape = ", roughnessNrmlzd.shape

if isinstance(polarization, basestring): 
    polarization = [polarization]

p_ws = polarization[0]

for p in polarization:
    oPath_3413_r = os.path.join(out_path, pp_names[0], p, 'epsg_3413', year, month, day, fileName)
   
    mkdirs(oPath_3413_r)
    oFileName = os.path.join(oPath_3413_r, fileName+'.png')
    if os.path.isfile(oFileName):
        continue

    roughness_res = pr.kd_tree.resample_nearest(swath_def, roughnessNrmlzd[p].ravel(), area_def,
                                            radius_of_influence=pxlResSARm.max(),
                                            epsilon=0.5, nprocs=numProcs, fill_value=None)
    print "resample_nearest done"

#        shapefile = '/media/SOLabNFS/store/auxdata/coastline/GSHHS_shp/f/GSHHS_f_L1.shp'
#        lonlim=(lons_2.min(),lons_2.max())
#        latlim=(lats_2.min(),lats_2.max())
#        lakes = True
#
#        mask_arr = gshhs_rasterize.gshhs_rasterize(lonlim, latlim, units, roughness_res.shape,
#                                                   proj_, proj_name, lakes, shapefile)
#        roughness_masked = ma.masked_where(mask_arr, roughness_res)

    gray()
    print 'save roughness_masked image, %s' % oFileName
#    imsave(oFileName, roughness_res, vmin=-1, vmax=1)
    save_big_image(roughness_res, oFileName, -1, 1, 4000)

    create_KML_asar(area_def.area_extent, os.path.join(oPath_3413_r, fileName+'.kml'))
    create_asar_tiles(oFileName, os.path.join(oPath_3413_r, 'tiles'), 'EPSG:3413')
    del roughness_res
    gc.collect()
#     del roughness_masked, roughness_res, roughness

del roughnessNrmlzd

oPath_3413_w = os.path.join(out_path, pp_names[1], p_ws, 'epsg_3413', year, month, day, fileName)
mkdirs(oPath_3413_w)
oFileNameWind = os.path.join(oPath_3413_w, fileName+'.png')
if not os.path.isfile(oFileNameWind):
    print "roughness_res done"
    wind_speed_asar_res = pr.kd_tree.resample_nearest(swath_def, wind_speed_asar.ravel(), area_def,
                                                      radius_of_influence=pxlResSARm.max(),
                                                      epsilon=0.5, nprocs=numProcs, fill_value=None)
    print "wind_speed_asar_res done"
    
    jet()
    print 'save wind_speed_asar_masked, image, %s' % oFileNameWind
    gc.collect()
    #imsave(oFileNameWind, wind_speed_asar_res, vmin=0, vmax=20)
    save_big_image(wind_speed_asar_res, oFileNameWind, 0, 20, 4000)
    gc.collect()
    
    #for _path in oPath:
    create_KML_asar(area_def.area_extent, os.path.join(oPath_3413_w, fileName+'.kml'))
    create_asar_tiles(oFileNameWind, os.path.join(oPath_3413_w, 'tiles'), 'EPSG:3413')
    del wind_speed_asar, wind_speed_asar_res
    gc.collect()

