# coding: utf-8

import sys
sys.path.append('/home/mag/Documents/repos/solab/posada/')
sys.path.append('/home/mag/Documents/repos/solab/posada/handlers')
sys.path.append('/home/mag/Documents/repos/solab/posada/handlers/Common')

from Common import distancelib

from numpy import floor, size, linspace, arange, row_stack, delete, round, \
                  double, asarray, array, fliplr, flipud, mean, pi, exp, \
                  sin, cos, radians, log10, ma, mean, repeat, zeros, reshape, \
                  float64, dot, sqrt, diff, concatenate

from math import atan

from multiprocessing import cpu_count

from scipy.interpolate import RectBivariateSpline, griddata

import pyresample as pr
from pyproj import Proj

from posada_logger import get_logger
logger = get_logger()

import os
import re
import datetime

import pygrib



def get_date_parameters(granule_name):
    ''' return year, month and day for current granule

    '''
    prog = re.compile(r'(\d{8})')
    file_date = prog.findall(granule_name)[0]
    prog = re.compile(r'(\d{6})')
    file_time = prog.findall(granule_name)[1]
    year = file_date[:4]
    month = file_date[4:6]
    day = file_date[6:]
    
    startTime = datetime.datetime.strptime(year+month+day+file_time,
    "%Y%m%d%H%M%S")

    return year, month, day, startTime

def mkdirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def format_extent_spacing(extent, spacing, extmax, midazimuth=False,
                           midrange=False):
    """Format (and check) extent and spacing."""
    # Check extent
    ext = round(extent).flatten()
    if ext.size != 4:
        raise Exception('extent must contain 4 elements')
    if (ext[0:2] < extmax[0:2]).any() or (ext[2:4] > extmax[2:4]).any():
        exttmp = array(ext)
        ext[0:2] = maximum(ext[0:2], extmax[0:2])
        ext[2:4] = minimum(ext[2:4], extmax[2:4])
        logger.warning('Warning : extent is outside SAR image, '+str(exttmp)+            ' becomes '+str(ext))
    if (ext[0:2] > ext[2:4]).any():
        raise Exception('extent[0:2] must be less or equal than '+                        'extent[2:4]')
    # Check spacing
    spa = round(spacing).flatten()
    if spa.size == 1:
        spa = repeat(spa[0], 2)
    elif spa.size == 2:
        pass
    else:
        raise Exception('spacing must contain 1 or 2 elements')
    if (spa < [1, 1]).any():
        spatmp = array(spa)
        spa = maximum(spa, [1, 1])
        logger.warning('Warning : spacing too small, '+str(spatmp)+' becomes '+            str(spa))
    if (spa > ext[2:4]-ext[0:2]).any():
        spatmp = array(spa)
        spa = minimum(spa, ext[2:4]-ext[0:2])
        logger.warning('Warning : spacing too large, '+str(spatmp)+' becomes '+            str(spa))
    # Make extent to be spacing modulo
    ext[2:4] -= (ext[2:4]-ext[0:2]) % spa
#     # 1D extent
#     if midazimuth == True:
#         dim = (ext[2]-ext[0]+1)/spa[0]
#         ext[0:3:2] = ext[0] + (dim-1)//2*spa[0] + [0, spa[0]-1]
#     if midrange == True:
#         dim = (ext[3]-ext[1]+1)/spa[1]
#         ext[1:4:2] = ext[1] + (dim-1)//2*spa[1] + [0, spa[1]-1]
    return (ext, spa)


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


def sigma02dB(data):
    """ Converting from linear units to dB """
    return 10*log10(data)


def swath_area_def(name='Temporal SWATH EPSG Projection 4326', proj='eqc', lonlim=(-180,180), latlim=(-90,90),
                   ellps="WGS84", res=111.2e3, lat_ts=None, lat_0=None, lon_0=None):
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

    if size(res) == 1:
        res = array((res, res))

    up    = max(latlim)
    down  = min(latlim)
    left  = min(lonlim)
    right = max(lonlim)

    # print 'up, down, left, right: ', round(up), round(down), round(left), round(right)

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
            logger.debug("11111111111") 
            area_extent = array((
                           min(left_ex1, left_ex2, right_ex1, right_ex2),
                           min(down_ex1, down_ex2, up_ex1, up_ex2),
                           max(left_ex1, left_ex2, right_ex1, right_ex2),
                           max(down_ex1, down_ex2, up_ex1, up_ex2)
                        ))
        elif (lon >=90 and lon <180) or (lon >=-270 and lon < -180):
            logger.debug("2222222222222")
            area_extent = array((
                           max(left_ex1, left_ex2, right_ex1, right_ex2),
                           max(down_ex1, down_ex2, up_ex1, up_ex2),
                           min(left_ex1, left_ex2, right_ex1, right_ex2),
                           min(down_ex1, down_ex2, up_ex1, up_ex2)
                        ))
        elif (lon >= 180 and lon < 270) or (lon >= -180 and lon < -90):
            logger.debug("333333333333")
            area_extent = array((
                           min(left_ex1, left_ex2, right_ex1, right_ex2),
                           min(down_ex1, down_ex2, up_ex1, up_ex2),
                           max(left_ex1, left_ex2, right_ex1, right_ex2),
                           max(down_ex1, down_ex2, up_ex1, up_ex2)
                        ))
        else:
            logger.debug("44444444444444444")
            area_extent = array((
                           min(left_ex1, left_ex2, right_ex1, right_ex2),
                           min(down_ex1, down_ex2, up_ex1, up_ex2),
                           max(left_ex1, left_ex2, right_ex1, right_ex2),
                           max(down_ex1, down_ex2, up_ex1, up_ex2)
                        ))
    else:
        # минимум из всех координат X, Y, максимум из всех координат X, Y
        # Такой результат даёт правильный area_extent для 3413
        # При этом для 4326 area_extent остаётся неизменным
        # area_def_3413 = swath_area_def(name='Temporal SWATH EPSG Projection 3413', proj='stere', \
        #                                lonlim=(-180,180), latlim=(30,90), ellps="WGS84", res=1500, \
        #                                lat_ts=70, lat_0=90, lon_0=-45)
        # Area extent: (-5050747.263141337, 0.0, 0.0, 5050747.263141336)
        area_extent = array((
                        min(left_ex1, left_ex2, right_ex1, right_ex2),
                        min(up_ex1, up_ex2, down_ex1, down_ex2),
                        max(left_ex1, left_ex2, right_ex1, right_ex2),
                        max(up_ex1, up_ex2, down_ex1, down_ex2)
                    ))

    modulox = (area_extent[2] - area_extent[0]) % res[0]
    moduloy = (area_extent[3] - area_extent[1]) % res[1]
    area_extent[0] = area_extent[0] + modulox
    area_extent[1] = area_extent[1] + moduloy
    area_extent[2] = area_extent[2]
    area_extent[3] = area_extent[3]

    # Using abs() to avoid negative numbers of coloumns/rows as for epsg3413 for example
    xsize = abs(int(round((area_extent[2] - area_extent[0]) / res[0])))
    ysize = abs(int(round((area_extent[3] - area_extent[1]) / res[1])))
    
    swath_area_def = pr.utils.get_area_def(area_id, name, proj_id, proj4_args, xsize, ysize, area_extent)

    return swath_area_def

def resample(data, lats_2, lons_2, latlim=None, lonlim=None, pxlRes=None, proj='EPSG:3413', numProcs=None):
    """
    Resampling data to specified projection
    """

    # #### Define areas
    logger.info("Defining areas.......")

    if latlim is None:
        latlim = (lats_2.min(), lats_2.max())
    if lonlim is None:
        lonlim = (lons_2.min(), lons_2.max())
    if pxlRes is None:
        pxlRes = 800
    if numProcs is None:
        numProcs=cpu_count()-1

    # Define areas with pyresample
    swath_def = pr.geometry.SwathDefinition(
        lons=lons_2, lats=lats_2
    )

    logger.info("Reprojecting to %s" % proj)
    if proj == 'EPSG:4326':
        area_def = swath_area_def(name='Temporal SWATH EPSG Projection 4326',
                                  proj='eqc',
                                  lonlim=lonlim,
                                  latlim=latlim, ellps="WGS84",
                                  res=pxlRes)
        # Set the parameters for GSHHS masking
        area_def.proj_ = '4326'
        area_def.proj_name = None
        area_def.units = 'deg'
    elif proj == 'EPSG:3413':
        area_def = swath_area_def(name='Temporal SWATH EPSG Projection 3413',
                                  proj='stere',
                                  lonlim=lonlim,
                                  latlim=latlim, ellps="WGS84",
                                  res=pxlRes,
                                  lat_ts=70, lat_0=90, lon_0=-45)
        # Set the parameters for GSHHS masking
        area_def.proj_ = '+units=m +ellps=WGS84 +lon_0=-45 +proj=stere +lat_ts=70 +lat_0=90'
        area_def.proj_name = '3413'
        area_def.units = 'm'


    # #### Recalculating Pixel resolution
    logger.info("Recalculating Pixel resolution.......")

    # Get the SAR pixel resolution from the area_def
    # for further identical shapes

    pxlRes = asarray(
        (abs(area_def.pixel_size_x), abs(area_def.pixel_size_y))
    )
    # logger.debug("       S1 cell resolution, %s m" % str(pxlRes))


    # #### Reprojecting non msked arrays
    logger.info("Reprojecting data to new projection(s).......")

    data_res = pr.kd_tree.resample_nearest(
        swath_def, data.ravel(), area_def,
        radius_of_influence=4*pxlRes.max(),
        epsilon=0.5, fill_value=None, nprocs=numProcs
    )

    lon = mean(lonlim)
    if proj == 'EPSG:3413' and (lon >=90 and lon <180) or (lon >=-270 and lon < -180):
        data_res = fliplr(flipud(data_res))

    return data_res, swath_def, area_def


def area_def_crop(area_def, new_shape, area_extent=None):
    """
    Trimming the area_def correspondingly to the normalized array
    (x_ll, y_ll, x_ur, y_ur)
    """
    area_id = area_def.area_id
    name = area_def.name
    proj_id = area_def.proj_id
    proj4_args = area_def.proj_dict
    xsize = new_shape[1]
    ysize = new_shape[0]

    area_def_crop = pr.utils.get_area_def(area_id, name, proj_id, proj4_args, xsize, ysize, area_extent)

    return area_def_crop

from Tiles.nctiles import create_nc_tiles, write_attrib_to_nc, array_size_normalize

def ncTilesMetainfo(data_res, lats_2, lons_2, area_def, latlim=None, lonlim=None):

    if latlim is None:
        latlim = (lats_2.min(), lats_2.max())
    if lonlim is None:
        lonlim = (lons_2.min(), lons_2.max())

    # #### Normalize array size to be a multyply of the tile size before timming the area
    old_shape = data_res.shape
    data_res = array_size_normalize(data_res, False)
    new_shape = data_res.shape

    # #### Trimming the area_def correspondingly to the normalized array
    lon = mean(lonlim)
    if (lon >=90 and lon <180) or (lon >=-270 and lon < -180):
        area_extent_ = (
            abs((new_shape[1])*area_def.pixel_size_y) + area_def.area_extent[2],
            abs((new_shape[0])*area_def.pixel_size_x) + area_def.area_extent[3],
            area_def.area_extent[2],
            area_def.area_extent[3]
            )
    else:
        area_extent_ = (
            abs((new_shape[1])*area_def.pixel_size_y) + area_def.area_extent[0],
            abs((new_shape[0])*area_def.pixel_size_x) + area_def.area_extent[1],
            area_def.area_extent[0],
            area_def.area_extent[1]
            )

    bbox_area_def = area_def_crop(area_def, new_shape, area_extent=area_extent_)

    from shapely.geometry import Polygon

    if lats_2[0,0] < lats_2[-1,0] and lats_2[0,-1] < lats_2[-1,-1]:
        geo_extent = (
        (lons_2[-1,-1],lats_2[-1,-1]),
        (lons_2[-1,0],lats_2[-1,0]),
        (lons_2[0,0],lats_2[0,0]),
        (lons_2[0,-1],lats_2[0,-1])
        )
    elif lats_2[0,0] > lats_2[-1,0] and lats_2[0,-1] > lats_2[-1,-1]:
        geo_extent = (
        (lons_2[0,0],lats_2[0,0]),
        (lons_2[0,-1],lats_2[0,-1]),
        (lons_2[-1,-1],lats_2[-1,-1]),
        (lons_2[-1,0],lats_2[-1,0])
        )

    geospatial_extent_ll = Polygon(geo_extent)

    # get list of available resolutions
    avail_resolution_list = [int(256*pow(2, i)) for i in range(15)]

    # get first guess size of area
    # size = int(256*math.pow(2, 14)/res)
    _size = abs(int(-5000000 - 5000000) / area_def.pixel_size_x)

    # get closest size of area from available resolutions list
    size = min(filter(lambda x: _size <= x,
                                  avail_resolution_list))

    # Extent is large enough, not to get errors about points outside area
    _area_def = pr.geometry.AreaDefinition(
               'epsg_3413_crude', 'NSIDC Polar Stereographic North EPSG:3413',
               'epsg_3413_crude',
               {'proj': 'stere', 'lat_0': '90',
                'lon_0': '-45', 'lat_ts': '70', 'ellps': 'WGS84',
                'datum': 'WGS84', 'units': 'm'}, size, size,
               [-6000000, -6000000, 6000000, 6000000]
           )

    try:
        geospatial_extent = (
            _area_def.get_proj_coords(data_slice=(_area_def.get_xy_from_lonlat(geo_extent[0][0], geo_extent[0][1])[::-1])),
            _area_def.get_proj_coords(data_slice=(_area_def.get_xy_from_lonlat(geo_extent[1][0], geo_extent[1][1])[::-1])),
            _area_def.get_proj_coords(data_slice=(_area_def.get_xy_from_lonlat(geo_extent[2][0], geo_extent[2][1])[::-1])),
            _area_def.get_proj_coords(data_slice=(_area_def.get_xy_from_lonlat(geo_extent[3][0], geo_extent[3][1])[::-1])),
            _area_def.get_proj_coords(data_slice=(_area_def.get_xy_from_lonlat(geo_extent[0][0], geo_extent[0][1])[::-1]))
            )
    except:
        # if most southern lat is less 45 degrees we take the most southern point available from the area_def    
        geospatial_extent = (
            _area_def.get_proj_coords(data_slice=(_area_def.get_xy_from_lonlat(geo_extent[0][0], geo_extent[0][1])[::-1])),
            _area_def.get_proj_coords(data_slice=(_area_def.get_xy_from_lonlat(geo_extent[1][0], geo_extent[1][1])[::-1])),
            _area_def.get_proj_coords(data_slice=((_area_def.shape[1]/2, _area_def.shape[0]))),
            _area_def.get_proj_coords(data_slice=((_area_def.shape[1]/2, _area_def.shape[0]))),
            _area_def.get_proj_coords(data_slice=(_area_def.get_xy_from_lonlat(geo_extent[0][0], geo_extent[0][1])[::-1]))
            )

    return data_res, bbox_area_def, geospatial_extent, geospatial_extent_ll










def findSubsetIndices(min_lat,max_lat,min_lon,max_lon,lats,lons):
    """Array to store the results returned from the function"""
    res=zeros((4),dtype=float64)
    minLon=min_lon; maxLon=max_lon

    distances1 = []; distances2 = []
    indices=[]; index=1

    for point in lats:
        s1 = max_lat-point # (vector subtract)
        s2 = min_lat-point # (vector subtract)
        distances1.append((dot(s1, s1), point, index))
        distances2.append((dot(s2, s2), point, index-1))
        index=index+1

    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])

    distances1 = []; distances2 = []; index=1

    for point in lons:
        s1 = maxLon-point # (vector subtract)
        s2 = minLon-point # (vector subtract)
        distances1.append((dot(s1, s1), point, index))
        distances2.append((dot(s2, s2), point, index-1))
        index=index+1

    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])

    # Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices
    minJ=indices[1][2]
    maxJ=indices[0][2]
    minI=indices[3][2]
    maxI=indices[2][2]

    res[0]=min(minI,maxI); res[1]=max(minI,maxI); res[2]=min(minJ,maxJ); res[3]=max(minJ,maxJ)
    return res


def windDirection(u, v):
    U = u.ravel()
    V = v.ravel()
    direction = zeros(size(U))
    for i in range(0, len(U)):
        if U[i] >= 0 and V[i] > 0:
            direction[i] = ((180 / pi) * atan(abs(U[i] / V[i])) + 180)
        if U[i] < 0 and V[i] > 0:
            direction[i] = (-(180 / pi) * atan(abs(U[i] / V[i])) + 180)
        if U[i] >= 0 and V[i] < 0:
            direction[i] = (-(180 / pi) * atan(abs(U[i] / V[i])) + 360)
        if U[i] < 0 and V[i] < 0:
            direction[i] = ((180 / pi) * atan(abs(U[i] / V[i])))
        if V[i] == 0 and U[i] > 0:
            direction[i] = 270
        if V[i] == 0 and U[i] < 0:
            direction[i] = 90
        if V[i] == 0 and U[i] == 0:
            direction[i] = 0
    return reshape(direction, v.shape)


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


from scipy.interpolate import RectSphereBivariateSpline

def ncepGFSmodel2swath(lats, lons, data, lats_2, lons_2):

    func = RectSphereBivariateSpline(lats, lons, data)
    data_2 = func.ev(lats_2.ravel(),\
                     lons_2.ravel())\
                     .reshape(lats_2.shape)
    return data_2


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


def addModelWind(startTime, lats_2, lons_2, numProcs=None):

    if numProcs is None:
        numProcs=cpu_count()-1

    ncepGFSmodelWind = ncepGFSmodel(startTime, lats_2, lons_2)
    ncep_def  = pr.geometry.GridDefinition (lons=ncepGFSmodelWind['lons_wind'], \
                                            lats=ncepGFSmodelWind['lats_wind'])
    swath_def = pr.geometry.SwathDefinition(lons=lons_2, lats=lats_2)
    pxlResWind = asarray(distancelib.getPixelResolution(ncepGFSmodelWind['lats_wind'], \
                                                        ncepGFSmodelWind['lons_wind'], \
                                                        ncepGFSmodelWind['lats_wind'].shape, 'km'))
    ncepGFSmodelWindSwath = {}
    ncepGFSmodelWindSwath['wind_speed'] = pr.kd_tree.resample_gauss(ncep_def, ncepGFSmodelWind['wind_speed'].ravel(), swath_def, \
                                         radius_of_influence=2*pxlResWind.min()*1e3, neighbours=12, \
                                         sigmas=pxlResWind.max()*1e3, fill_value=None, nprocs=numProcs)
    ncepGFSmodelWindSwath['wind_dir']   = pr.kd_tree.resample_gauss(ncep_def, ncepGFSmodelWind['wind_dir'].ravel(), swath_def, \
                                         radius_of_influence=2*pxlResWind.min()*1e3, neighbours=12, \
                                         sigmas=pxlResWind.max()*1e3, fill_value=None, nprocs=numProcs)
    ncepGFSmodelWindSwath['u'] = pr.kd_tree.resample_gauss(ncep_def, ncepGFSmodelWind['u'].ravel(), swath_def, \
                                         radius_of_influence=2*pxlResWind.min()*1e3, neighbours=12, \
                                         sigmas=pxlResWind.max()*1e3, fill_value=None, nprocs=numProcs)
    ncepGFSmodelWindSwath['v'] = pr.kd_tree.resample_gauss(ncep_def, ncepGFSmodelWind['v'].ravel(), swath_def, \
                                         radius_of_influence=2*pxlResWind.min()*1e3, neighbours=12, \
                                         sigmas=pxlResWind.max()*1e3, fill_value=None, nprocs=numProcs)
    return ncepGFSmodelWindSwath