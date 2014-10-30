__author__ = 'Alexander Myasoedov and Denis Spiridonov'

import simplekml
import os, datetime
import numpy
from numpy import *
from pylab import *
import pyresample as pr
from pyproj import Proj
import epr
import numpy.polynomial.polynomial as poly

def create_KML_asar(area_extent, savepath):
    kml = simplekml.Kml()

    pol = kml.newpolygon(name='area_extent', visibility=1)
    pol.tessellate = 1

    pol.altitudemode = 'clampToGround'
    pol.outerboundaryis.coords = [(area_extent[0], area_extent[1]), (area_extent[2], area_extent[3])]
    kml.save(savepath)

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

    up    = min(latlim)
    down  = max(latlim)
    left  = min(lonlim)
    right = max(lonlim)

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

    area_extent = (min(left_ex1, left_ex2),
                   min(up_ex1, up_ex2),
                   max(right_ex1, right_ex2),
                   max(down_ex1, down_ex2))

#     Using abs() to avoid negative numbers of coloumns/rows as for epsg3413 for example
    xsize = abs(int((area_extent[2] - area_extent[0]) / res))
    ysize = abs(int((area_extent[3] - area_extent[1]) / res))

    swath_area_def = pr.utils.get_area_def(area_id, name, proj_id, proj4_args, xsize, ysize, area_extent)

#     print swath_area_def

    return swath_area_def


def create_asar_image(iPath, oPath_4326, oPath_3413, fileName):
    #plt.close("all")
    product = epr.Product(os.path.join(iPath, fileName))

    band = product.get_band('proc_data')

    sc_w = product.get_scene_width()
    sc_h = product.get_scene_height()

    raw_counts = band.read_as_array(sc_w, sc_h)#, xoffset=100, yoffset=6500, xstep=2, ystep=2)
    lat = product.get_band('latitude').read_as_array(sc_w, sc_h)
    lon = product.get_band('longitude').read_as_array(sc_w, sc_h)
    incident_angle = product.get_band('incident_angle').read_as_array(sc_w, sc_h)

    raw_counts_trmmd = raw_counts
    # Trimming the array by removing zero values from rows and cols
    msk = []
    for m in range(raw_counts_trmmd.shape[0]):
        if raw_counts_trmmd[m, :].sum() == 0:
            msk.append(m)
    raw_counts_trmmd = numpy.delete(raw_counts_trmmd, msk, axis=0)
    lat = numpy.delete(lat, msk, axis=0)
    lon = numpy.delete(lon, msk, axis=0)
    incident_angle = numpy.delete(incident_angle, msk, axis=0)

    msk = []
    for n in range(raw_counts_trmmd.shape[1]):
        if raw_counts_trmmd[:, n].sum() == 0:
            msk.append(n)
    raw_counts_trmmd = numpy.delete(raw_counts_trmmd, msk, axis=1)
    lat = numpy.delete(lat, msk, axis=1)
    lon = numpy.delete(lon, msk, axis=1)
    incident_angle = numpy.delete(incident_angle, msk, axis=1)
    raw_counts = raw_counts_trmmd

    # Adding Sigma_0
    calibration_constant = \
    product.get_dataset('MAIN_PROCESSING_PARAMS_ADS').read_record(0).get_field('calibration_factors.1.ext_cal_fact').get_elems()
    # sigma0 = 10*log10( raw_counts**2*sin(incident_angle*pi/180)/calibration_constant )
    sigma0 = raw_counts**2*sin(incident_angle*pi/180)/calibration_constant

    print "    start filter"
    from scipy.signal import wiener
    sigma0w = wiener(sigma0, mysize=(7,7), noise=None)
    # sigma0w = sigma0

    pol = product.get_sph().get_field('MDS1_TX_RX_POLAR').get_elem()
    if pol == 'H/H':
        ph = (2.20495, -14.3561e-2, 11.28e-4)
        sigma0_hh_ref = exp((ph[0]+incident_angle*ph[1]+incident_angle**2*ph[2])*log(10))
        roughness = sigma0w/sigma0_hh_ref
    elif pol == 'V/V':
        pv = (2.29373, -15.393e-2, 15.1762e-4)
        sigma0_vv_ref = exp((pv[0]+incident_angle*pv[1]+incident_angle**2*pv[2])*log(10))
        roughness = sigma0w/sigma0_vv_ref

    # masking the arrays
    raw_counts = ma.masked_where(raw_counts == 0, raw_counts)
    roughness = ma.masked_where(raw_counts == 0, roughness)
    sigma0 = ma.masked_where(raw_counts == 0, sigma0)
    sigma0w = ma.masked_where(raw_counts == 0, sigma0w)
    default_fill_value = double(ma.default_fill_value(raw_counts))

    scale = 30

    pr.kd_tree.which_kdtree()
    pr.get_capabilities()

    scale = 1

    ln = lon[::scale,::scale]
    lt = lat[::scale,::scale]
    data = roughness[::scale,::scale]

    for proj in ['EPSG:4326', 'EPSG:3413']:
        print "    start projection %s" % proj
        if proj == 'EPSG:4326':
            oPath = oPath_4326
            area_def = swath_area_def(name='Temporal SWATH EPSG Projection 4326', proj='eqc',
                                      lonlim=(lon.min(),lon.max()), latlim=(lat.min(),lat.max()),
                                      ellps="WGS84", res=150)
        elif proj == 'EPSG:3413':
            oPath = oPath_3413
            area_def = swath_area_def(name='Temporal SWATH EPSG Projection 3413', proj='stere',
                                      lonlim=(lon.min(),lon.max()), latlim=(lat.min(),lat.max()),
                                      ellps="WGS84", res=150, lat_ts=70, lat_0=90, lon_0=-45)

        swath_def = pr.geometry.SwathDefinition(lons=ln, lats=lt)
        # result = pr.kd_tree.resample_nearest(swath_def, data.ravel(), area_def, radius_of_influence=30000, nprocs=2)
        # result = pr.kd_tree.resample_nearest(swath_def, data, area_def, radius_of_influence=50000, epsilon=0.5, fill_value=default_fill_value)
        result = pr.kd_tree.resample_nearest(swath_def, data.ravel(), area_def, radius_of_influence=300,
                                             epsilon=0.5, nprocs=2, fill_value=None)

        oFileName = os.path.join(oPath, fileName+'.png')
        gray()
        imsave(oFileName, result, vmin=0, vmax=2)

        create_KML_asar(area_def.area_extent, os.path.join(oPath, fileName+'.kml'))
    close()
    # gray()
    # pr.plot.show_quicklook(area_def, result, vmin=0, vmax=2, label='Test', num_meridians=45, num_parallels=10, coast_res='l')

import gdal
import sys
sys.path.append('/usr/bin')
from gdal2tiles import GDAL2Tiles

def create_asar_tiles(png_filename, tiles_output_dir, proj):
    local_argv = ['/usr/bin/gdal2tiles.py', '-p', 'raster', '-r', 'cubic',
                  '-s', proj, png_filename, tiles_output_dir]
    argv = gdal.GeneralCmdLineProcessor(local_argv)
    if argv:
        gdal2tiles = GDAL2Tiles(argv[1:])
        gdal2tiles.process()

from PIL import Image
def resaze_image(filepath, basewidth, savepath):
    img = Image.open(filepath)
    wpercent = (basewidth / float(img.size[0]))
    hsize = int(float(img.size[1]) * float(wpercent))
    img = img.resize((basewidth, hsize), Image.ANTIALIAS)
    img.save(savepath)

fileName = 'ASA_WSM_1PNPDE20110815_090644_000001903105_00309_49461_8731.N1'
asar_path = '/nfs1/store/satellite/asar'
pp_name = 'sigma0'

for _dir, sub_dir, _files in os.walk(asar_path):
    for fileName in _files:
        if fileName.startswith('ASA_') and fileName.endswith('.N1'):
            dt = datetime.datetime.strptime(fileName[14:22], '%Y%m%d').timetuple()
            year = str(dt.tm_year)
            day = str(dt.tm_yday)

            if len(day) == 1:
                day = '00' + day
            elif len(day) == 2:
                day = '0' + day
            print year, day

            out_path = '/tmp/ASAR'

            _dir = '/nfs1/store/satellite/asar/%s/%s' % (year, day)

            oPath_4326 = os.path.join(out_path, pp_name, '4326', year, day, fileName)
            oPath_3413 = os.path.join(out_path, pp_name, '3413', year, day, fileName)

            print "Start granule: %s" % fileName

            # create dirs
            if not os.path.isdir(oPath_4326):
                os.makedirs(oPath_4326)
            if not os.path.isdir(oPath_3413):
                os.makedirs(oPath_3413)
            # check png file
            png_4326_filename = os.path.join(oPath_4326, fileName+'.png')
            png_3413_filename = os.path.join(oPath_3413, fileName+'.png')
            if not os.path.isfile(png_4326_filename) or not os.path.isfile(png_3413_filename):
                print "Start create image for %s" % fileName
                create_asar_image(_dir, oPath_4326, oPath_3413, fileName)

            # create thumbs
            resaze_image(png_4326_filename, 1024, os.path.join(oPath_4326, fileName+'_1024.png'))
            resaze_image(png_4326_filename, 263, os.path.join(oPath_4326, fileName+'_263.png'))
            resaze_image(png_3413_filename, 1024, os.path.join(oPath_3413, fileName+'_1024.png'))
            resaze_image(png_3413_filename, 263, os.path.join(oPath_3413, fileName+'_263.png'))

            # check tiles
            tiles_4326_output_dir = os.path.join(oPath_4326, 'tiles')
            tiles_3413_output_dir = os.path.join(oPath_3413, 'tiles')

            if not os.path.isdir(tiles_4326_output_dir):
                create_asar_tiles(png_4326_filename, tiles_4326_output_dir, 'EPSG:4326')

            if not os.path.isdir(tiles_3413_output_dir):
                create_asar_tiles(png_3413_filename, tiles_3413_output_dir, 'EPSG:3413')
