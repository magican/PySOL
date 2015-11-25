# coding: utf-8
import datetime
import os
import re
import sys
import gc

# import gdal
from numpy import asarray, zeros, ma, flipud, pi, exp, cos, radians, log10
# import matplotlib
# import matplotlib.pyplot as plt
# from matplotlib.mlab import find
# matplotlib.use('Agg')
from pylab import *

# from PIL import Image
from scipy.signal import wiener
import pyresample as pr
# from pyproj import Proj
from scipy.interpolate import RectSphereBivariateSpline
import simplekml
from netCDF4 import Dataset as ncDataset

import readS1
from readS1 import *
import distancelib
import gshhs_rasterize

# sys.path.append('/usr/bin')
# from gdal2tiles import GDAL2Tiles

sys.path.append(
    os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )
)

from big_image import save_big_image
from redis_messages import redis_publish_metadata, redis_publish_newgranule
from Tiles.nctiles import create_nc_tiles

__author__ = 'Alexander Myasoedov'
__email__ = 'mag@rshu.ru'
__created__ = datetime.datetime(2014, 10, 28)
__modified__ = datetime.datetime(2015, 6, 19)
__version__ = "1.0"
__status__ = "Development"


def mkdirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)


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

    P_0 = A_0 * exp(B_0 * theta) + C_0
    P_HALF_PI = A_HALF_PI * exp(B_HALF_PI * theta) + C_HALF_PI
    P_PI = A_PI * exp(B_PI * theta) + C_PI

    C0 = (P_0 + P_PI + 2 * P_HALF_PI) / 4
    C1 = (P_0 - P_PI) / 2
    C2 = (P_0 + P_PI - 2 * P_HALF_PI) / 4

    P = C0 + C1 * cos(radians(phi)) + C2 * cos(radians(2 * phi))

    return P


def ncepGFSmodel2swath(lats, lons, data, lats_2, lons_2):
    func = RectSphereBivariateSpline(lats, lons, data)
    data_2 = func.ev(lats_2.ravel(),\
                     lons_2.ravel())\
                     .reshape(lats_2.shape)
    return data_2


def processingS1(inpath, fn, scale=1, resolution=80):
    #~ Setting default resolution to 80m, which in most cases sorresponds to scale=8
    pxlResSARm = resolution
    s1 = readS1(inpath=inpath, fn=fn, resolution=pxlResSARm, min_lat=35)
    if s1.lat_skip:
        return None
    # s1.__dict__['raw_counts']

    # get vars from s1 class
    # polarization
    # raw_counts
    # incidenceAngle_2
    # sigmaNought_2
    # noiseLut_2
    # sigma0
    # lons_2
    # lats_2
    # GEOgrid
    # cLUTs
    # nLUTs
    # manifest
    polarization = s1.polarization
    incidenceAngle_2 = s1.incidenceAngle_2
    sigma0 = s1.sigma0
    lons_2 = s1.lons_2
    lats_2 = s1.lats_2
    GEOgrid = s1.GEOgrid
    manifest = s1.manifest

    del s1
    gc.collect()

    sigma0w = {}
    roughness = {}

    print "Scale set to: \'%s\' " % scale

    for p in polarization:
        print "Filtering Image: \'%s\' polarization" % p

        # filter the image
        sigma0w[p] = wiener(
            sigma0[p][::scale, ::scale], mysize=(7, 7), noise=None
        )
        # sigma0w[p] = sigma0[p]

    del sigma0
    gc.collect()

    # S1 Pixel resolution
    # we use pxlResSAR for further GSHHS rasterizing and
    # reprojecting data with pyresample

    lonlim = (lons_2[::scale, ::scale].min(), lons_2[::scale, ::scale].max())
    latlim = (lats_2[::scale, ::scale].min(), lats_2[::scale, ::scale].max())

    # enlarge lonlims for cropping a bit larger area for masking
    lonlimGSHHS = (lonlim[0]-1.0, lonlim[1]+1.0)
    latlimGSHHS = (latlim[0]-1.0, latlim[1]+1.0)

    # Get first guess pixel resolution
    pxlResSARm = asarray(
        distancelib.getPixelResolution(
            lats_2[::scale, ::scale], lons_2[::scale, ::scale],
            lons_2[::scale, ::scale].shape, 'km'
        )
    )*1e3
    pxlResSARdeg = asarray(
        distancelib.getPixelResolution(
            lats_2[::scale, ::scale], lons_2[::scale, ::scale],
            lons_2[::scale, ::scale].shape, 'deg'
        )
    )

    # Define areas with pyresample
    swath_def = pr.geometry.SwathDefinition(
        lons=lons_2[::scale, ::scale], lats=lats_2[::scale, ::scale]
    )

    area_def_4326 = swath_area_def(name='Temporal SWATH EPSG Projection 4326',
                                   proj='eqc', lonlim=lonlimGSHHS,
                                   latlim=latlimGSHHS, ellps="WGS84",
                                   res=pxlResSARm)

    # Get the SAR pixel resolution from the area_def
    # for further identical shapes
    up = min(latlimGSHHS)
    down = max(latlimGSHHS)
    left = min(lonlimGSHHS)
    right = max(lonlimGSHHS)
    area_extent_deg = (left, down, right, up)

    area_extent_deg_shape = area_def_4326.shape

    pxlResSARdeg = asarray(
        (abs(area_extent_deg[2] - area_extent_deg[0]) /
            float(area_extent_deg_shape[1]),
         abs(area_extent_deg[3] - area_extent_deg[1]) /
            float(area_extent_deg_shape[0]))
    )

    pxlResSARm = asarray(
        (area_def_4326.pixel_size_x, area_def_4326.pixel_size_y)
    )
    print "S1 cell resolution, %s deg" % str(pxlResSARdeg)
    print "S1 cell resolution, %s m" % str(pxlResSARm)

    # Apply Mask from GSHHS
    reload(gshhs_rasterize)

    # ESRI shapefile containing land polygons
    shapefile = '/media/SOLabNFS/store/auxdata/coastline/GSHHS_shp/f/GSHHS_f_L1.shp'

    # reproject GSHHS onto S1 grid before calculations
    print "Rasterizing Land Mask"
    mask_arr_4326 = gshhs_rasterize.gshhs_rasterize_4326(
        lonlimGSHHS, latlimGSHHS, pxlResSARdeg, area_def_4326.shape,
        True, shapefile
    )

    del pxlResSARdeg

    mask_arr_swath = pr.kd_tree.resample_nearest(
        area_def_4326, mask_arr_4326, swath_def,
        radius_of_influence=4*pxlResSARm.max(), epsilon=0.5, fill_value=None
    )


    # Commented until roughness is needed
    # # Nice Image (Roughness)
    # sigma0wAvg = {}
    # roughnessNrmlzd = {}

    # if len(polarization[0]) >= 2:  # if 2 polarizations
    #     for p in polarization:
    #         print "Nice Image: \'%s\' polarization" % p

        # NB!!! DO NOT MASK LAND for Sigma0w and Roughness!!!!

    #         roughness[p] = ma.masked_where(mask_arr_swath, sigma0w[p])
    #         sigma0wAvg[p] = ma.median(roughness[p], axis=0)
    #         roughnessNrmlzd[p] = (roughness[p]-sigma0wAvg[p])/sigma0wAvg[p]
    # elif len(polarization[0]) == 1:  # if only 1 polarization
    #     p = polarization
    #     print "Nice Image: \'%s\' polarization" % p
    #     roughness[p] = ma.masked_where(mask_arr_swath, sigma0w[p])
    #     sigma0wAvg[p] = ma.median(roughness[p], axis=0)
    #     roughnessNrmlzd[p] = (roughness[p]-sigma0wAvg[p])/sigma0wAvg[p]
    # del roughness, sigma0wAvg


    # Adding Model wind

    # import xmltodict

    # zf = zipfile.ZipFile(inpath+fn, 'r')
    # manifest = zf.read(fn[:-4] + '.SAFE/manifest.safe')
    # manifest = xmltodict.parse(manifest) # Parse the read document string
    # zf.close()

    Objects = manifest['xfdu:XFDU']['metadataSection']['metadataObject']
    for Object in Objects:
        try:
            startTime = datetime.datetime.strptime(
                Object['metadataWrap']['xmlData']['safe:acquisitionPeriod']['safe:startTime'],
                "%Y-%m-%dT%H:%M:%S.%f"
            )
            break
        except:
            pass

    ncepGFSmodelWind = ncepGFSmodel(startTime, lats_2, lons_2)

    # Reprojecting data
    # Pixel resolution
    # we use pxlResWind/pxlResSAR for further pyresample
    # radius_of_influence and sigmas
    pxlResWind = asarray(
        distancelib.getPixelResolution(
            ncepGFSmodelWind['lats_wind'],
            ncepGFSmodelWind['lons_wind'],
            ncepGFSmodelWind['lons_wind'].shape, 'km'
        )
    )

    # Note pxlResWind is in KM, multiply by 1e3 for meters
    print "S1 cell resolution, %s m" % pxlResSARm
    print "Wind cell resolution, %s km" % pxlResWind

	# reproject NCEP onto S1 grid before calculations
	# Using RectSphereBivariateSpline - Bivariate spline approximation over a rectangular mesh on a sphere
	# as it is much more efficiant for full resolution
	# as well as smoothes nicely the image

	# We don't want to work with full res wind so scaling the image for about 100m resolution
	# Adjust scale to get appropriate value
	# scale = 5

	lts = ncepGFSmodelWind['lats_wind']
	lns = ncepGFSmodelWind['lons_wind']

	ncepGFSmodelWindSwath = {}
	# check that lns increasing, if not
	if ~all(diff(lns[0,:]) > 0):
	    lns_ = lns[0,:]
	    # we must start lns from -180 increasing to +180
	    # so we reconcatenate lns array and data array, by cropping part of array and putting in front
	    lns_ = concatenate((lns_[(lns_>=-180) & (lns_<0)], lns_[(lns_>=0)]),0)
	    ncepGFSmodelWindSwath['wind_speed'] = (ncepGFSmodelWind['wind_speed'])
	    ncepGFSmodelWindSwath['wind_speed'] = concatenate((ncepGFSmodelWind['wind_speed'][:,(lns_>=-180) & (lns_<0)],\
	                                                  ncepGFSmodelWind['wind_speed'][:,(lns_>=0)]),1)
	    ncepGFSmodelWindSwath['wind_dir'] = (ncepGFSmodelWind['wind_dir'])
	    ncepGFSmodelWindSwath['wind_dir'] = concatenate((ncepGFSmodelWind['wind_dir'][:,(lns_>=-180) & (lns_<0)],\
	                                                  ncepGFSmodelWind['wind_dir'][:,(lns_>=0)]),1)
	    # only then we do interpolate

	lts_2 = lats_2[::scale,::scale]
	lns_2 = lons_2[::scale,::scale]

	# RectSphereBivariateSpline uses lats and lons within the intervals (0, pi), (0, 2pi).
	# Make sure the latitude is between 0 .. 180
	if lts_2.min()<0:
	    lts_2 = lts_2 + 90
	    lts   = lts   + 90
	if lns_2.min()<0:
	    lns_2 = lns_2 + 180
	    lns   = lns   + 180

	try:
	    ncepGFSmodelWindSwath['wind_speed'] = griddata((lts.flatten(), lns.flatten()), (ncepGFSmodelWind['wind_speed']).flatten(), (lats_2, lons_2), method='cubic')
	    ncepGFSmodelWindSwath['wind_dir'] = griddata((lts.flatten(), lns.flatten()), (ncepGFSmodelWind['wind_dir']).flatten(), (lats_2, lons_2), method='cubic')
	    ncepGFSmodelWindSwath['u'] = griddata((lts.flatten(), lns.flatten()), (ncepGFSmodelWind['u']).flatten(), (lats_2, lons_2), method='cubic')
	    ncepGFSmodelWindSwath['v'] = griddata((lts.flatten(), lns.flatten()), (ncepGFSmodelWind['v']).flatten(), (lats_2, lons_2), method='cubic')
	except:
	#     ncepGFSmodelWindSwath['wind_speed']  = RectBivariateSpline(lts, lns, flipud(ncepGFSmodelWindSwath['wind_speed']), kx=2, ky=2)(lts_2, lns_2)
	    pass


	# FOR SOME REASON RectSphereBivariateSpline doesn't work very good, using griddata instead

    pxlResWindSwath = asarray(distancelib.getPixelResolution(lts_2, \
                                                        lns_2, \
                                                        lns_2.shape, 'km'))

    print "Interpolated Wind cell resolution, %s km" % pxlResWindSwath

    del lts, lns, lts_2, lns_2

    # calculate bearing from initial lats/lons for further wind calculation
    # Taking initial values as bearing is more accurate after
    # interpolation than vice versa
    bearing = zeros((GEOgrid['lons'].shape[0]-1, GEOgrid['lons'].shape[1]))

    for n in range(0, GEOgrid['lons'].shape[1]):
        col = ([GEOgrid['lats'][:-1, n], GEOgrid['lons'][:-1, n]],
               [GEOgrid['lats'][1:, n], GEOgrid['lons'][1:, n]])
        for m in range(0, GEOgrid['lons'].shape[0]-1):
            bearing[m][n] = distancelib.bearing(
                asarray(col[0])[:, m], asarray(col[1])[:, m]
            )

    # interpolate to raw_counts.shape
    bearing_2 = imresize(bearing, ncepGFSmodelWindSwath['wind_dir'].shape)

    # NB! WINDDIR = 0 WHEN WIND BLOWS TOWARDS RADAR!
    p = polarization[0]

    wind_dir_model_swath_rel = 90 + bearing_2 -\
        ncepGFSmodelWindSwath['wind_dir']

    del bearing, bearing_2

    if p == 'hh':
        PR = PR_Mouche(
            incidenceAngle_2[p][::scale, ::scale], wind_dir_model_swath_rel
        )
        try:
            from cmod.cmod_gpu import rcs2windOpenCl
            wind_speed_asar = rcs2windOpenCl(
                sar=sigma0w[p]*PR, windir=wind_dir_model_swath_rel,
                theta=incidenceAngle_2[p][::scale, ::scale]
            )
        except Exception:
            from cmod.cmod_vect import rcs2windPar
            wind_speed_asar = rcs2windPar(
                sigma0w[p]*PR, cmdv=5, windir=wind_dir_model_swath_rel,
                theta=incidenceAngle_2[p][::scale, ::scale], nprocs=numProcs
            )
        del PR
    elif p == 'vv':
        try:
            from cmod.cmod_gpu import rcs2windOpenCl
            wind_speed_asar = rcs2windOpenCl(
                sar=sigma0w[p], windir=wind_dir_model_swath_rel,
                theta=incidenceAngle_2[p][::scale, ::scale]
            )
        except Exception:
            from cmod.cmod_vect import rcs2windPar
            wind_speed_asar = rcs2windPar(
                sigma0w[p], cmdv=5, windir=wind_dir_model_swath_rel,
                theta=incidenceAngle_2[p][::scale, ::scale], nprocs=numProcs
            )

    gc.collect()

    # Add mask to initial NCEP data
    area_def_ncep = pr.geometry.SwathDefinition(
        lons=ncepGFSmodelWind['lons_wind'], lats=ncepGFSmodelWind['lats_wind']
    )
    mask_arr_ncep = pr.kd_tree.resample_nearest(
        area_def_4326, mask_arr_4326, area_def_ncep,
        radius_of_influence=4*pxlResWind.max(), epsilon=0.5, fill_value=None
    )

    del area_def_4326, mask_arr_4326, pxlResWind, area_def_ncep

    ncepGFSmodelWind['wind_speed'] = ma.masked_where(
        mask_arr_ncep, ncepGFSmodelWind['wind_speed']
    )

    # Add mask to ASAR wind and reprojected NCEP
    wind_speed_asar = ma.masked_where(mask_arr_swath, wind_speed_asar)
    ncepGFSmodelWindSwath['wind_speed'] = ma.masked_where(
        mask_arr_swath, ncepGFSmodelWindSwath['wind_speed']
    )

    del mask_arr_swath, ncepGFSmodelWind, mask_arr_ncep, ncepGFSmodelWindSwath
    gc.collect()

    return_values = {}
    return_values['lats'] = lats_2
    return_values['lons'] = lons_2
    # return_values['roughnessNrmlzd'] = roughnessNrmlzd
    return_values['sigma0w'] = sigma0w
    return_values['wind_speed_asar'] = wind_speed_asar
    return_values['pxlResSARm'] = pxlResSARm
    return_values['polarization'] = polarization
    return_values['swath_def'] = swath_def
    return_values['startTime'] = startTime
    return return_values


def create_KML(area_extent, savepath):
    kml = simplekml.Kml()

    pol = kml.newpolygon(name='area_extent', visibility=1)
    pol.tessellate = 1

    pol.altitudemode = 'clampToGround'
    # minx, miny, maxx, maxy
    pol.outerboundaryis.coords = [(min(area_extent[0], area_extent[2]),
                                   min(area_extent[1], area_extent[3])),
                                  (max(area_extent[0], area_extent[2]),
                                   max(area_extent[1], area_extent[3]))]
    if type(savepath) == list:
        for _savepath in savepath:
            kml.save(_savepath)
    else:
        kml.save(savepath)


# def create_asar_tiles(png_filename, tiles_output_dir, proj):
#     local_argv = ['/usr/bin/gdal2tiles.py', '-p', 'raster', '-r', 'cubic',
#                   '-s', proj, png_filename, tiles_output_dir]
#     argv = gdal.GeneralCmdLineProcessor(local_argv)
#     if argv:
#         gdal2tiles = GDAL2Tiles(argv[1:])
#         gdal2tiles.process()


def write_attrib_to_nc(nc_path, area_def, start_time, max_zoom_level,
                       resolution, polarizations):
    dataset = ncDataset(nc_path, 'a', format='NETCDF4')
    area_extent = area_def.area_extent
    # minx, miny, maxx, maxy
    dataset.setncattr('geospatial_lat_min',
                      min(area_extent[0], area_extent[2]))
    dataset.setncattr('geospatial_lon_min',
                      min(area_extent[1], area_extent[3]))
    dataset.setncattr('geospatial_lat_max',
                      max(area_extent[0], area_extent[2]))
    dataset.setncattr('geospatial_lon_max',
                      max(area_extent[1], area_extent[3]))

    dataset.setncattr('start_date', start_time.isoformat())

    dataset.setncattr('polarizations', ','.join(polarizations))
    dataset.setncattr('max_zoom_level', str(max_zoom_level))
    dataset.setncattr('resolution', str(resolution))


def publish_redis_messages(nc_path, area_def, max_zoom_level, resolution,
                           polarizations):
    redis_publish_newgranule(nc_path)

    area_extent = area_def.area_extent
    BBox_attrib = [
        min(area_extent[0], area_extent[2]),
        min(area_extent[1], area_extent[3]),
        max(area_extent[0], area_extent[2]),
        max(area_extent[1], area_extent[3]),
    ]
    resolution_list = [int(resolution * math.pow(2, i)) for i in range(max_zoom_level)]
    resolution_list.reverse()

    polarizations = map(lambda x: str(x), polarizations)

    message = {
        'PRODUCT_NAME': 'SOLAB-SENTINEL-1',
        'GRANULE_NAME': nc_path,
        'OUTPUT_DIRECTORY': os.path.basename(nc_path),
        'METADATA': {
            'bbox': BBox_attrib,
            'resolution': resolution_list,
            'variables': {
                # 'roughness': {
                #     'polarization': polarizations,
                #     'min': -1,
                #     'max': 1,
                # },
                'sigma0w': {
                    'polarization': polarizations,
                    'min': -35,
                    'max': 5,
                },
                'wind_speed': {
                    'polarization': [polarizations[0]],
                    'min': 0,
                    'max': 35,
                },
            }
        },
    }
    redis_publish_metadata(message)


def create_s1_nc(inpath, fn, out_dir, scale=1):
    if '_SLC_' in fn:
        print 'skipping SLC product for now'
        return
    prog = re.compile(r'(\d{8})')

    file_date = prog.findall(fn)[0]
    year = file_date[:4]
    month = file_date[4:6]
    day = file_date[6:]

    granule_name = fn.replace('.zip', '')
    nc_path = os.path.join(out_dir, year, month, day, granule_name+'.nc')
    if os.path.isfile(nc_path):
        print 'nc file exists: %s' % nc_path
        return

    print 'Start granule: ', fn

    resolution = 80
    return_values = processingS1(inpath, fn, scale, resolution)
    if return_values is None:
        return

    print 'Preprocessing done'
    # for key in return_values:
    #     locals()[key] = return_values[key]

    lats = return_values['lats']
    lons = return_values['lons']
    # roughnessNrmlzd = return_values['roughnessNrmlzd']
    sigma0w = return_values['sigma0w'] # convert to dB
    wind_speed_asar = return_values['wind_speed_asar']
    pxlResSARm = return_values['pxlResSARm']
    polarizations = return_values['polarization']
    swath_def = return_values['swath_def']
    start_time = return_values['startTime']

    # nc_path = os.path.join('/tmp', fn+'.nc')
    if not os.path.isdir(os.path.dirname(nc_path)):
        os.makedirs(os.path.dirname(nc_path))

    area_def = swath_area_def(
        name='Temporal SWATH EPSG Projection 3413', proj='stere',
        lonlim=(lons.min(), lons.max()), latlim=(lats.min(), lats.max()),
        ellps="WGS84", res=pxlResSARm, lat_ts=70, lat_0=90, lon_0=-45
    )
    pxlResSARm_max = pxlResSARm.max()

    # # roughness
    # for p in polarizations:
    #     print "Start polar: ", p
    #     roughness_res = pr.kd_tree.resample_nearest(
    #         swath_def, roughnessNrmlzd[p].ravel(), area_def,
    #         radius_of_influence=pxlResSARm_max,
    #         epsilon=0.5, nprocs=numProcs, fill_value=None
    #     )
    #     create_nc_tiles(roughness_res, 'roughness', nc_path, p)

    # sigma0w
    for p in polarizations:
        print "Start resampling polarisation: ", p
        sigma0w_res = pr.kd_tree.resample_nearest(
            swath_def, 10*log10(sigma0w[p].ravel()), area_def,
            radius_of_influence=pxlResSARm_max,
            epsilon=0.5, nprocs=numProcs, fill_value=None
        )
        create_nc_tiles(sigma0w_res, 'sigma0w', nc_path, p)

    # wind speed
    p_ws = polarizations[0]

    print "Start resampling wind speed, polarisation: ", p
    wind_speed_res = pr.kd_tree.resample_nearest(
        swath_def, wind_speed_asar.ravel(), area_def,
        radius_of_influence=pxlResSARm_max,
        epsilon=0.5, nprocs=numProcs, fill_value=None
    )

    max_zoom_level = create_nc_tiles(wind_speed_res, 'wind_speed', nc_path,
                                     p_ws)

    write_attrib_to_nc(nc_path, area_def, start_time,
                       max_zoom_level, resolution, polarizations)
    publish_redis_messages(nc_path, area_def, max_zoom_level, resolution,
                           polarizations)


if __name__ == '__main__':
    inpath = '/media/SOLabNFS2/tmp/different_SAR/sentinel-1/Ania_Ladoga_29_May_2015/'
    fileNameList = ['S1A_IW_GRDH_1SDV_20150603T154002_20150603T154027_006211_0081A9_5F10.zip',
                    'S1A_IW_GRDH_1SDV_20150529T041657_20150529T041722_006131_007F51_F751.zip',
                    'S1A_EW_GRDM_1SDH_20150517T153117_20150517T153221_005963_007AED_56B0.zip']

    fn = fileNameList[0]

    out_dir = '/media/SOLabNFS2/store/satellite/SOLAB-SENTINEL-1'

    create_s1_nc(inpath, fn, out_dir, scale=1)


# # Ania_Ladoga_29_May_2015/
# inpath = '/media/SOLabNFS2/tmp/different_SAR/sentinel-1/Ania_Ladoga_29_May_2015/'
# fileNameList = ['S1A_IW_GRDH_1SDV_20150603T154002_20150603T154027_006211_0081A9_5F10.zip',
#                 'S1A_IW_GRDH_1SDV_20150529T041657_20150529T041722_006131_007F51_F751.zip',
#                 'S1A_EW_GRDM_1SDH_20150517T153117_20150517T153221_005963_007AED_56B0.zip']

# fn = fileNameList[0]

# prog = re.compile(r'(\d{8})')
# file_date = prog.findall(fn)[0]

# year = file_date[:4]
# month = file_date[4:6]
# day = file_date[6:]
# # day = '04'

# print 'Start granule: ', fn

# return_values = processingS1(inpath, fn)

# print 'preprocessing done'
# for key in return_values:
#     locals()[key] = return_values[key]

# # lats_2 = return_values['lats_2']
# # lons_2 = return_values['lons_2']
# # roughnessNrmlzd = return_values['roughnessNrmlzd']
# # wind_speed_asar = return_values['wind_speed_asar']
# # pxlResSARm = return_values['pxlResSARm']
# # polarization = return_values['polarization']
# # swath_def = return_values['swath_def']

# del return_values

# numProcs = multiprocessing.cpu_count() - 2
# # for proj in ['EPSG:3413']:
# proj = 'EPSG:3413'

# pp_names = ['roughness', 'wind_speed']
# # out_path = '/media/SOLabNFS2/http/tiles/Sentinel-1'
# out_path = '/tmp/Sentinel-1'

# # year = '2015'
# # month = '06'
# # day = '03'
# fileName = fn[:-4]

# area_def = swath_area_def(
#     name='Temporal SWATH EPSG Projection 3413', proj='stere',
#     lonlim=(lons.min(), lons.max()), latlim=(lats.min(), lats.max()),
#     ellps="WGS84", res=pxlResSARm, lat_ts=70, lat_0=90, lon_0=-45
# )

# del lats, lons
# gc.collect()

# pxlResSARm_max = pxlResSARm.max()
# del pxlResSARm

# # Set the parameters for GSHHS masking
# proj_ = '+units=m +ellps=WGS84 +lon_0=-45 +proj=stere +lat_ts=70 +lat_0=90'
# proj_name = '3413'
# units = 'm'

# print area_def
# # print "roughness.shape = ", roughnessNrmlzd.shape

# if isinstance(polarization, basestring):
#     polarization = [polarization]

# p_ws = polarization[0]

# for p in polarization:
#     oPath_3413_r = os.path.join(out_path, pp_names[0], p, 'epsg_3413',
#                                 year, month, day, fileName)

#     mkdirs(oPath_3413_r)
#     oFileName = os.path.join(oPath_3413_r, fileName+'.png')
#     if os.path.isfile(oFileName):
#         continue

#     print 'Start roughness: %s polarization' % p
#     roughness_res = pr.kd_tree.resample_nearest(
#         swath_def, roughnessNrmlzd[p].ravel(), area_def,
#         radius_of_influence=pxlResSARm_max,
#         epsilon=0.5, nprocs=numProcs, fill_value=None
#     )
#     del roughnessNrmlzd[p]
#     print "resample_nearest done"

#     # shapefile = '/media/SOLabNFS/store/auxdata/coastline/GSHHS_shp/f/GSHHS_f_L1.shp'
#     # lonlim=(lons_2.min(),lons_2.max())
#     # latlim=(lats_2.min(),lats_2.max())
#     # lakes = True

#     # mask_arr = gshhs_rasterize.gshhs_rasterize(
#     #     lonlim, latlim, units, roughness_res.shape,
#     #     proj_, proj_name, lakes, shapefile
#     # )
#     # roughness_masked = ma.masked_where(mask_arr, roughness_res)

#     gray()
#     print 'save roughness_masked image, %s' % oFileName
# #    imsave(oFileName, roughness_res, vmin=-1, vmax=1)
#     save_big_image(oFileName, roughness_res, vmin=-1, vmax=1)

#     create_KML_asar(area_def.area_extent,
#                     os.path.join(oPath_3413_r, fileName+'.kml'))
#     create_asar_tiles(oFileName, os.path.join(oPath_3413_r, 'tiles'),
#                       'EPSG:3413')
#     del roughness_res
#     gc.collect()
# #     del roughness_masked, roughness_res, roughness

# del roughnessNrmlzd, polarization
# gc.collect()

# oPath_3413_w = os.path.join(out_path, pp_names[1], p_ws, 'epsg_3413',
#                             year, month, day, fileName)
# mkdirs(oPath_3413_w)
# oFileNameWind = os.path.join(oPath_3413_w, fileName+'.png')
# if not os.path.isfile(oFileNameWind):

#     print 'Start wind_speed: %s polarization' % p
#     wind_speed_asar_res = pr.kd_tree.resample_nearest(
#         swath_def, wind_speed_asar.ravel(), area_def,
#         radius_of_influence=pxlResSARm_max,
#         epsilon=0.5, nprocs=numProcs, fill_value=None
#     )
#     del swath_def, wind_speed_asar
#     print "wind_speed_asar_res done"

#     jet()
#     print 'save wind_speed_asar_masked, image, %s' % oFileNameWind
#     gc.collect()
#     # imsave(oFileNameWind, wind_speed_asar_res, vmin=0, vmax=20)
#     save_big_image(oFileNameWind, wind_speed_asar_res, vmin=0, vmax=20)
#     del wind_speed_asar_res
#     gc.collect()

#     # for _path in oPath:
#     create_KML_asar(area_def.area_extent,
#                     os.path.join(oPath_3413_w, fileName+'.kml'))
#     create_asar_tiles(oFileNameWind,
#                       os.path.join(oPath_3413_w, 'tiles'), 'EPSG:3413')
#     del area_def
#     gc.collect()
