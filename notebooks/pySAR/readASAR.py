# coding: utf-8

import datetime

import os, sys

import epr

from CommonSAR import *


def readASAR(iPath, fileName, pxlRes=800.0, image_size_skip=40*1e6, area_bbox=45):

    logger.info(os.path.join(iPath, fileName))
    try:
        product = epr.Product(os.path.join(iPath, fileName))
    except:
        logger.error('unable to read file')
        return False

    try:
        band = product.get_band('proc_data')
    except epr.EPRValueError:
        logger.error('unable to get band "proc_data": epr_get_band_id: band not found')
        return False

    sc_w = double(product.get_scene_width())
    sc_h = double(product.get_scene_height())

    # Skipping if image too large
    logger.debug('sc_w*sc_h = %s MPs', str(round(sc_w * sc_h / 1e6)))
    if sc_w*sc_h > image_size_skip:
        logger.debug("ASAR Image too large, skipping...")
        return False

    # Get lat/lon from geolocation grid
    dataset = product.get_dataset('GEOLOCATION_GRID_ADS')
    fltp_lats = map(
        lambda x:
        dataset.read_record(x).get_field('first_line_tie_points.lats').get_elems(),
        range(dataset.get_num_records())
    )
    lltp_lats = map(
        lambda x:
        dataset.read_record(x).get_field('last_line_tie_points.lats').get_elems(),
        range(dataset.get_num_records())
    )
    fltp_lons = map(
        lambda x:
        dataset.read_record(x).get_field('first_line_tie_points.longs').get_elems(),
        range(dataset.get_num_records())
    )
    lltp_lons = map(
        lambda x:
        dataset.read_record(x).get_field('last_line_tie_points.longs').get_elems(),
        range(dataset.get_num_records())
    )

    fltp_lats = asarray(double(fltp_lats))/1e6
    lltp_lats = asarray(double(lltp_lats))/1e6
    fltp_lons = asarray(double(fltp_lons))/1e6
    lltp_lons = asarray(double(lltp_lons))/1e6

    lats = row_stack((fltp_lats, lltp_lats[-1, :]))
    lons = row_stack((fltp_lons, lltp_lons[-1, :]))

    # Skipping if no area overlap
    if mean(lats[:]) <= area_bbox:
        logger.debug("No area overlap. Skipping...")
        return False

    lats = fliplr(lats)
    lons = fliplr(lons)

    # Find scale to reduce image to the specified resolution
    arrShape =  asarray([sc_w, sc_h])
    _lats = asarray([lats[0,0], lats[-1,-1], lats[0,-1], lats[-1,0]])
    _lons = asarray([lons[0,0], lons[-1,-1], lons[0,-1], lons[-1,0]])
    imageRes = round(mean(asarray(distancelib.getPixelResolution(_lats, _lons, arrShape, 'km'))*1e3))
    scale = pxlRes/imageRes

    extMax = (0.,0.,arrShape[0]-1,arrShape[1]-1)
    ext    = (0.,0.,arrShape[0]-1,arrShape[1]-1)

    # Format extent/spacing
    ext, spa = format_extent_spacing(extent=ext, spacing = scale, extmax=extMax)
    
    # Read data with stepping=spacing
    try:
        raw_counts = band.read_as_array(sc_w, sc_h, xstep=spa[0], ystep=spa[1])
        incident_angle = product.get_band('incident_angle').read_as_array(sc_w, sc_h, xstep=spa[0], ystep=spa[1])
    except epr.EPRValueError:
        logger.error("EPRValueError")
        return False

    
    lats_2 = imresize(lats, raw_counts.shape)
    lons_2 = imresize(lons, raw_counts.shape)

#     if lats.max() <= 35:
#         logger.debug("skipping no area overlap")
#         return False

    # Trimming the array by removing zero values from rows and cols
    msk = []
    for m in range(raw_counts.shape[0]):
        if raw_counts[m, :].sum() == 0:
            msk.append(m)
    raw_counts = delete(raw_counts, msk, axis=0)
    lats_2 = delete(lats_2, msk, axis=0)
    lons_2 = delete(lons_2, msk, axis=0)
    incident_angle = delete(incident_angle, msk, axis=0)
    polarization = product.get_sph().get_field('MDS1_TX_RX_POLAR').get_elem()

    msk = []
    for n in range(raw_counts.shape[1]):
        if raw_counts[:, n].sum() == 0:
            msk.append(n)
    raw_counts = delete(raw_counts, msk, axis=1)
    lats_2 = delete(lats_2, msk, axis=1)
    lons_2 = delete(lons_2, msk, axis=1)
    incident_angle = delete(incident_angle, msk, axis=1)

    # Adding Sigma_0
    calibration_constant = product.get_dataset('MAIN_PROCESSING_PARAMS_ADS').read_record(0).get_field('calibration_factors.1.ext_cal_fact').get_elems()
    # sigma0 = 10*log10( raw_counts**2*sin(incident_angle*pi/180)/calibration_constant )
    sigma0 = raw_counts**2*sin(incident_angle*pi/180)/calibration_constant
    
    return sigma0, lats_2, lons_2, incident_angle, polarization, lats, lons