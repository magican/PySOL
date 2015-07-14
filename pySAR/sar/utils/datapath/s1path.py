#!/usr/bin/env python
# coding=utf-8
"""
"""


import os
import glob
from datetime import datetime, timedelta, date
import re
import xml.etree.ElementTree as ET
import numpy as np
from shapely.geometry import Polygon, box
import warnings
#import pdb


DEFAULTROOT = '/home/cercache/project/mpc-sentinel1/data/esa'
MISSION = ['S1A', 'S1B']
PRODUCT = ['GRD', 'OCN', 'RAW', 'SLC']
MODE = ['EW', 'EW1', 'EW2', 'EW3', 'EW4', 'EW5',
        'IW', 'IW1', 'IW2', 'IW3',
        'SM', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6',
        'WV', 'WV1', 'WV2']
RESOLUTION = ['F', 'H', 'M', '_']
PCLASS = ['A', 'C', 'N', 'S']
POLARISATION = ['HH', 'HV', 'VH', 'VV']


def mission2missiondir(mission):
    """
    """
    missiondir = [mis.lower().replace('s', 'sentinel-') for mis in mission]
    return missiondir


def product2leveldir(product):
    """
    """
    prod2level = {'GRD':'L1', 'OCN':'L2', 'RAW':'L0', 'SLC':'L1'}
    leveldir = [prod2level[prod] for prod in product]
    return list(set(leveldir))


def mode2modedir(mode):
    """
    """
    modedir = []
    for mod in mode:
        if re.match(r'^EW[1-5]?$', mod) is not None:
            modedir.append('EW')
        elif re.match(r'^IW[1-3]?$', mod) is not None:
            modedir.append('IW')
        elif re.match(r'^S(M|[1-6])$', mod) is not None:
            modedir.append('SM')
        elif re.match(r'^WV[1-2]?$', mod) is not None:
            modedir.append('WV')
    return list(set(modedir))


def mode2modesafe(mode):
    """
    """
    modesafe = []
    for mod in mode:
        if re.match(r'^EW[1-5]?$', mod) is not None:
            modesafe.append('EW')
        elif re.match(r'^IW[1-3]?$', mod) is not None:
            modesafe.append('IW')
        elif mod == 'SM':
            modesafe += ['S'+str(num) for num in range(1, 7)]
        elif re.match(r'^S[1-6]$', mod) is not None:
            modesafe.append(mod)
        elif re.match(r'^WV[1-2]?$', mod) is not None:
            modesafe.append('WV')
    return list(set(modesafe))


def mode2modemsrt(mode):
    """
    """
    modemsrt = []
    for mod in mode:
        if mod == 'EW':
            modemsrt.append('ew')
            modemsrt += ['ew'+str(num) for num in range(1, 6)]
        elif re.match(r'^EW[1-5]$', mod) is not None:
            modemsrt.append(mod.lower())
        elif mod == 'IW':
            modemsrt.append('iw')
            modemsrt += ['iw'+str(num) for num in range(1, 4)]
        elif re.match(r'^IW[1-3]$', mod) is not None:
            modemsrt.append(mod.lower())
        elif mod == 'SM':
            modemsrt += ['s'+str(num) for num in range(1, 7)]
        elif re.match(r'^S[1-6]$', mod) is not None:
            modemsrt.append(mod.lower())
        elif mod == 'WV':
            modemsrt += ['wv'+str(num) for num in range(1, 3)]
        elif re.match(r'^WV[1-2]$', mod) is not None:
            modemsrt.append(mod.lower())
    return list(set(modemsrt))


def polarisation2polarisationsafe(polarisation):
    """
    """
    polarisationsafe = []
    for polar in polarisation:
        if polar == 'HH':
            polarisationsafe += ['SH', 'DH', 'HH']
        elif polar == 'VV':
            polarisationsafe += ['SV', 'DV', 'VV']
        elif polar == 'HV':
            polarisationsafe += ['DH', 'HV']
        elif polar == 'VH':
            polarisationsafe += ['DV', 'VH']
    return list(set(polarisationsafe))


def polarisation2polarisationmsrt(polarisation):
    """
    """
    polarisationmsrt = [polar.lower() for polar in polarisation]
    return polarisationmsrt


def split_path(path):
    """
    """
    spath = path.split(os.path.sep)
    while '' in spath:
        spath.remove('')
    return spath


def ydn2date(year, daynumber):
    """
    """
    return date(year, 01, 01) + timedelta(days=daynumber-1)


# def date2ydn(date):
#     """
#     """
#     timetup = date.timetuple()
#     return (timetup.tm_year, timetup.tm_yday)


def intersect_bbox(footprint, bounding_box):
    """
    """
    footprint[1:, 0] += np.round((footprint[0, 0]-footprint[1:, 0])/360.)*360
    meanbblon = (bounding_box[0] + bounding_box[2]) / 2.
    meanfplon = (footprint[:, 0].min() + footprint[:, 0].max()) / 2.
    footprint[:, 0] += np.round((meanbblon - meanfplon) / 360.) * 360
    bbox = box(bounding_box[0], bounding_box[1], bounding_box[2],
               bounding_box[3])
    footp = Polygon(footprint)
    return bbox.intersects(footp)


def intersect_bpoly(footprint, bounding_polygon):
    """
    """
    footprint[1:, 0] += np.round((footprint[0, 0]-footprint[1:, 0])/360.)*360
    bounding_lon = [point[0] for point in bounding_polygon]
    meanbplon = (min(bounding_lon) + max(bounding_lon)) / 2.
    meanfplon = (footprint[:, 0].min() + footprint[:, 0].max()) / 2.
    footprint[:, 0] += np.round((meanbplon - meanfplon) / 360.) * 360
    bpoly = Polygon(bounding_polygon)
    footp = Polygon(footprint)
    return bpoly.intersects(footp)


def safefootprint(safepath):
    """
    """
    footprint = []
    namespaces = {'safe': 'http://www.esa.int/safe/sentinel-1.0',
                  'gml': 'http://www.opengis.net/gml'}
    path = './metadataSection/metadataObject/metadataWrap/xmlData/' \
           'safe:frameSet/safe:frame/safe:footPrint/gml:coordinates'
    tree = ET.parse(os.path.join(safepath, 'manifest.safe'))
    root = tree.getroot()
    elems = root.findall(path, namespaces=namespaces)
    for elem in elems:
        corners = np.array(re.split(r',| ', elem.text), dtype='float32')
        footprint.append(corners.reshape((4, 2))[:, ::-1])
    return footprint


# def test_footprint_number():
#     for mod in ['WV', 'IW', 'EW', 'SM']:
#         for prod in ['SLC', 'GRD']:
#             print mod, prod
#             safepath = safepath_from_filters(mode=[mod], product=[prod],
#                                              pclass='S',
#                                              date_stop=datetime(2014, 04, 30))
#             print len(safefootprint(safepath[0]))


def safepath2msrtpath(safepath, mode=None, polarisation=None,
                      date_start=None, date_stop=None,
                      bounding_box=None, bounding_polygon=None,
                      **kwargs):
    """
    """
    # Make msrt test
    notest = lambda x: True
    if mode is None:
        modetest = notest
    else:
        modemsrt = mode2modemsrt(mode)
        modetest = lambda x: x[1] in modemsrt
    if polarisation is None:
        polartest = notest
    else:
        polarmsrt = polarisation2polarisationmsrt(polarisation)
        polartest = lambda x: x[3] in polarmsrt
    if date_start is None:
        date0test = notest
    else:
        date0test = lambda x: datetime.strptime(x[5], '%Y%m%dt%H%M%S') \
                    >= date_start
    if date_stop is None:
        date1test = notest
    else:
        date1test = lambda x: datetime.strptime(x[4], '%Y%m%dt%H%M%S') \
                    <= date_stop
    msrttest = lambda x: modetest(x) and polartest(x) and date0test(x) and \
               date1test(x)
    # Find msrt
    msrtpath = []
    for spath in safepath:
        safeproduct = os.path.basename(spath)[7:10]
        if safeproduct in ['GRD', 'SLC']:
            msrts = glob.glob(os.path.join(spath, 'measurement/*.tiff'))
        elif safeproduct == 'OCN':
            msrts = glob.glob(os.path.join(spath, 'measurement/*.nc'))
        elif safeproduct == 'RAW':
            warnings.warn('RAW *.dat are not filtered at all !')
            msrtpath += glob.glob(os.path.join(spath, '*.dat'))
            continue
        else:
            raise Exception('Unkown SAFE product : {}'.format(safeproduct))
        safemsrtpath = []
        for msrt in msrts:
            smsrt = os.path.basename(msrt).split('-')
            if msrttest(smsrt):
                safemsrtpath.append(msrt)
        if (bounding_box is not None or bounding_polygon is not None) and \
           len(safemsrtpath) != 0:
            footprint = safefootprint(spath)
            if os.path.basename(spath)[4:6] == 'WV':
                for msrt in safemsrtpath:
                    num = int(msrt[-8:-5])
                    if bounding_polygon is not None:
                        if intersect_bpoly(footprint[num-1], bounding_polygon):
                            msrtpath.append(msrt)
                    else:
                        if intersect_bbox(footprint[num-1], bounding_box):
                            msrtpath.append(msrt)
            else:
                if bounding_polygon is not None:
                    if intersect_bpoly(footprint[0], bounding_polygon):
                        msrtpath += safemsrtpath
                else:
                    if intersect_bbox(footprint[0], bounding_box):
                        msrtpath += safemsrtpath
        else:
            msrtpath += safemsrtpath
    return msrtpath


def safepath_from_safename(safename, root=DEFAULTROOT):
    """
    """
    safepath = []
    for safe in safename:
        missiondir = safe[0:3].lower().replace('s', 'sentinel-')
        leveldir = 'L'+safe[12]
        typedir = safe[0:14]
        if typedir[4] == 'S':
            typedir[5] = 'M'
        modedir = typedir[4:6]
        yeardir = safe[17:21]
        sdate = date(int(yeardir), int(safe[21:23]), int(safe[23:25]))
        ydaydir = '%03i' % sdate.timetuple().tm_yday
        spath = os.path.join(root, missiondir, leveldir, modedir, typedir,
                             yeardir, ydaydir, safe)
        if os.path.exists(spath):
            safepath.append(spath)
    return safepath


def safepath_from_filters(root=DEFAULTROOT, mission=None, product=None,
                          mode=None, resolution=None, pclass=None,
                          polarisation=None, orbit_start=None, orbit_stop=None,
                          date_start=None, date_stop=None, **kwargs):
    """
    """
    # Make walking tests according to inputs
    dowalk = []
    notest = lambda x: True
    ### root
    dowalk += [notest]
    ### {'sentinel-1a', 'sentinel-1b'}
    if mission is None:
        dowalk += [notest]
    else:
        missiondir = mission2missiondir(mission)
        dowalk += [lambda x: x[-1] in missiondir]
    ### {'L0', 'L1'}
    if product is None:
        dowalk += [notest]
    else:
        leveldir = product2leveldir(product)
        dowalk += [lambda x: x[-1] in leveldir]
    ### {'EW', 'IW', 'SM', 'WV'}
    if mode is None:
        dowalk += [notest]
    else:
        modedir = mode2modedir(mode)
        dowalk += [lambda x: x[-1] in modedir]
    ### {'S1A_WV_GRDM_1S', 'S1A_WV_SLC__1A', 'S1A_WV_SLC__1S', ...}
    if product is None:
        producttest = notest
    else:
        producttest = lambda x: x[-1][7:10] in product
    if resolution is None:
        resolutiontest = notest
    else:
        resolutiontest = lambda x: x[-1][10] in resolution
    if pclass is None:
        pclasstest = notest
    else:
        pclasstest = lambda x: x[-1][-1] in pclass
    dowalk += [lambda x: producttest(x) and resolutiontest(x) and pclasstest(x)]
    ### year
    if date_start is None:
        year0test = notest
    else:
        year0test = lambda x: int(x[-1]) >= date_start.year
    if date_stop is None:
        year1test = notest
    else:
        year1test = lambda x: int(x[-1]) <= date_stop.year
    dowalk += [lambda x: year0test(x) and year1test(x)]
    ### day number
    if date_start is None:
        day0test = notest
    else:
        day0test = lambda x: ydn2date(int(x[-2]), int(x[-1])) >= date_start.date()
    if date_stop is None:
        day1test = notest
    else:
        day1test = lambda x: ydn2date(int(x[-2]), int(x[-1])) <= date_stop.date()
    dowalk += [lambda x: day0test(x) and day1test(x)]
    ### safe
    if mode is None:
        modetest = notest
    else:
        modesafe = mode2modesafe(mode)
        modetest = lambda x: x[-1][4:6] in modesafe
    if polarisation is None:
        polartest = notest
    else:
        polarsafe = polarisation2polarisationsafe(polarisation)
        polartest = lambda x: x[-1][14:16] in polarsafe
    if date_start is None:
        date0test = notest
    else:
        date0test = lambda x: datetime.strptime(x[-1][33:48], '%Y%m%dT%H%M%S') \
                    >= date_start
    if date_stop is None:
        date1test = notest
    else:
        date1test = lambda x: datetime.strptime(x[-1][17:32], '%Y%m%dT%H%M%S') \
                    <= date_stop
    if orbit_start is None:
        orbit0test = notest
    else:
        orbit0test = lambda x: int(x[-1][49:55]) >= orbit_start
    if orbit_stop is None:
        orbit1test = notest
    else:
        orbit1test = lambda x: int(x[-1][49:55]) <= orbit_stop
    dowalk += [lambda x: modetest(x) and polartest(x) and date0test(x) and \
               date1test(x) and orbit0test(x) and orbit1test(x)]
    # Walk
    safepath = []
    rootdepth = len(split_path(root))
    levelstop = len(dowalk)-1
    for path, dirs, files in os.walk(root, topdown=True):
        spath = split_path(path)
        level = len(spath) - rootdepth
        if level < levelstop:
            if not dowalk[level](spath):
                del dirs[:]
        else:
            del dirs[:]
            if dowalk[level](spath):
                safepath.append(path)
    return safepath


def s1path(root=DEFAULTROOT, safename=None, onlysafe=False, **kwargs):
    """
    """
    if safename is not None:
        safepath = safepath_from_safename(safename, root=root)
    else:
        safepath = safepath_from_filters(root=root, **kwargs)
    if onlysafe == True:
        return safepath
    else:
        msrtpath = safepath2msrtpath(safepath, **kwargs)
        return msrtpath
