# -*- coding: utf-8 -*-

"""
Description :
    Give the geolocation of a MODIS 250m pixel specified by its indices [i_250m_along_track, i_250m_cross_track] for products where the (latitude, longitude) datasets are given at 1km resolution
    For the interpolation of the 250m positions, it takes care of borders of scan. In this cases, the position will be extrapolated.
    It also deals with the change date meridian crossing

Usage :
    python modis_250m_subpixel_geolocation.py <infile> <i_250m_along_track> <i_250m_cross_track>

    It will print out the ( lat, lon ) geolocation of the given pixel

    Optionally, this script can also display a figure of the pixels around the requested one.

Prerequisites :
    [REQUIRED]
        - python >= 2.5
        - numpy
        - pyhdf

    [OPTIONAL]
        When enabling the plot of the grid ( by setting __DEBUG__ to True ), you need to have installed the matplotlib library wirth a fonctionnal graphical toolkit support

Author :
    CGTD-ICARE/UDEV Nicolas PASCAL ( nicolas.pascal-at-icare.univ-lille1.fr  )

License :
    This file must be used under the terms of the CeCILL.
    This source file is licensed as described in the file COPYING, which
    you should have received as part of this distribution.  The terms
    are also available at
    http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

History :
--------
v0.1.0 : 2009/11/06
 - creation
"""

import math
import sys

import warnings
warnings.filterwarnings("ignore")

import numpy as np
from pyhdf import SD

__DEBUG__ = False
if __DEBUG__ is True :
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

# --- 1km -> 250m ---
sz_sc_1km  = 1354
sz_sc_250m = 5416

def itk_1km_to_250m ( i_tk_1km ) :
    """
    return the 250m grid index along track of a 1km pixel
    """
    return 1.5 + 4. * i_tk_1km

def itk_250m_to_1km ( i_tk_250m ) :
    """
    return the 1km grid index along track of a 250m pixel
    """
    return ( i_tk_250m - 1.5 ) / 4.

def isc_1km_to_250m ( i_sc_1km ) :
    """
    return the 250m grid index cross track of a 1km pixel
    """
    return 4. * i_sc_1km

def isc_250m_to_1km ( i_sc_250m ) :
    """
    return the 1km grid index cross track of a 250m pixel
    """
    return i_sc_250m / 4.

def get_250m_pix_pos ( itk_250m, isc_250m, lat_1km, lon_1km ) :
    """
    return the (lat,lon) of a 250m pixel specified with its indexes, when the geolocation datasets are given at 1km resolution
    @param itk_250m grid index of the 250m pixel along-track
    @param isc_250m grid index of the 250m pixel cross-track
    @param lat_1km latitudes dataset at 1km resolution
    @param lon_1km longitudes dataset at 1km resolution
    @return the ( lat, lon ) of the 250m pixel  [ itk_250m, isc_250m ]
    """
    # change date meridian case sentinel
    is_crossing_change_date = False

    # check 250m indexes validity
    sz_tk_1km = lat_1km.shape [0]
    sz_tk_250m = 4 * sz_tk_1km
    if ( isc_250m < 0 ) or ( isc_250m > sz_sc_250m - 1 ) :
        raise ValueError ( "Invalid scan index %d. Must be in range [%d,%d]"%(isc_250m,0,sz_sc_250m-1) )
    if ( itk_250m < 0 ) or ( itk_250m > sz_tk_250m - 1 ) :
        raise ValueError ( "Invalid track index %d. Must be in range [%d,%d]"%(itk_250m,0,sz_tk_250m-1) )

    # --- set the bounding 1km pixels to take for interpolation
    # set the (track,scan) indexes of the 1km pixel in the 250m grid
    itk_1km = itk_250m_to_1km ( itk_250m )
    isc_1km = isc_250m_to_1km ( isc_250m )

    #print "i_250m=[%.2f, %.2f] -> i_1km=[%.2f, %.2f]"%(itk_250m, isc_250m, itk_1km,isc_1km)

    # the width of one scan, in number of pixels, at 250m resolution
    w_scan_250m = 40
    # - extra/interpolation 1km pixels along track
    if (   itk_250m % w_scan_250m ) <= 1 : # extrapolate start of scan
        itk_top_1km    = math.ceil ( itk_1km )
        itk_bottom_1km = itk_top_1km + 1
    elif ( itk_250m % w_scan_250m ) >= 38 : # extrapolate end of scan
        itk_top_1km    = math.floor ( itk_1km )
        itk_bottom_1km = itk_top_1km - 1
    else : # general case : middle of scan
        itk_top_1km    = math.floor ( itk_1km )
        itk_bottom_1km = itk_top_1km + 1
    # - extra/interpolation 1km pixels along track
    if ( isc_1km >= 1353. ) : # extrapolate end of scan line
        isc_left_1km  = math.floor ( isc_1km ) - 1
        isc_right_1km = math.floor ( isc_1km )
    else : # general case : interpolation
        isc_left_1km  = math.floor ( isc_1km )
        isc_right_1km = isc_left_1km + 1

    #print "itk_top_1km=%d itk_bottom_1km=%d isc_left_1km=%d isc_right_1km=%d"%(itk_top_1km, itk_bottom_1km, isc_left_1km, isc_right_1km)

    # --- set the 1km track lines position ; left border ---
    lat_left_250m, lon_left_250m   = get_y_pos_1km_to_250m (
                isc_left_1km, itk_top_1km, itk_bottom_1km,
                lat_1km, lon_1km,
                itk_250m )
    # --- set the 1km track lines position ; right border ---
    lat_right_250m, lon_right_250m = get_y_pos_1km_to_250m (
                isc_right_1km, itk_top_1km, itk_bottom_1km,
                lat_1km, lon_1km,
                itk_250m )

    #print "left_250m=[%f,%f] right_250m=[%f,%f]"%(lat_left_250m, lon_left_250m, lat_right_250m, lon_right_250m)

    # check for change date meridian case
    if  abs ( lon_right_250m  - lon_left_250m  ) > 180. :
        is_crossing_change_date = True
        # all negative longitudes will be incremented of 360 before interpolation
        if lon_left_250m < 0. :
            lon_left_250m += 360.
        if lon_right_250m < 0. :
            lon_right_250m += 360.

    # for each track line position, interpolate along scan to retrieve the 250m geolocation
    lat, lon = get_x_pos_1km_to_250m ( lat_left_250m,  lon_left_250m,  isc_left_1km,
                                     lat_right_250m, lon_right_250m, isc_right_1km,
                                     isc_250m )
    #print "geolocation = [%f, %f]"%(lat,lon)

    # in case of change date crossing, turn values > 180. to negative
    if lon > 180.:
        lon -= 360.
    elif lon < -180.:
        lon += 360.
    return lat, lon

def get_x_pos_1km_to_250m ( lat_left_250m,  lon_left_250m,  isc_left_1km,
                            lat_right_250m, lon_right_250m, isc_right_1km,
                            isc_250m
                        ) :
    """
    retrieve the position of a 250m pixel set by its cross-track index and the position of 2 along-track lines set by 2 sucessive 1km
    pixels
    PRECONDITION : change date meridian case must have already been treated, and so lon_left_250m or lon_right_250m can have values > 180.
    @param lat_left_250m latitude of the 250m along track pixel that is the intersection between cross-track line at index isc_250m and
    the along-track line that joins the 2 sucessive bounding 1km pixels on the left of the point to interpolate
    @param lon_left_250m longitude of the 250m along track pixel that is the intersection between cross-track line at index isc_250m and
    the along-track line that joins the 2 sucessive bounding 1km pixels on the left of the point to interpolate
    @param isc_left_1km cross-track index at 1km of the left border
    @param lat_right_250m latitude of the 250m along track pixel that is the intersection between cross-track line at index isc_250m and
    the along-track line that joins the 2 sucessive bounding 1km pixels on the right of the point to interpolate
    @param lon_right_250m longitude of the 250m along track pixel that is the intersection between cross-track line at index isc_250m and
    the along-track line that joins the 2 sucessive bounding 1km pixels on the right of the point to interpolate
    @param isc_right_1km cross-track index at 1km of the right border
    @return (lat,lon) of the pixel at [itk_250m, isc_1km], interpolated between the 2 successive 1km pixels p1 and p2
    """
    # make sure P1 and P1 are 2 successive pixels along-track
    if abs ( isc_left_1km - isc_right_1km ) != 1 :
        raise ValueError ( "The 2 borders are on the same along-track line and must be successive" )

    # coordinates of left and right border in the 250m grid
    isc_left_250m  = isc_1km_to_250m ( isc_left_1km  )
    isc_right_250m = isc_1km_to_250m ( isc_right_1km )

    #print "isc_left_1km=%d isc_right_1km=%d"%(isc_left_1km, isc_right_1km)
    #print "isc_left_250m=%d isc_right_250m=%d"%(isc_left_250m, isc_right_250m)
    #print "isc_250m_min=%d isc_250m_max=%d"%(isc_250m_min,isc_250m_max)

    # linear interpolation on the position
    alpha_lon = ( lon_right_250m - lon_left_250m ) / ( isc_right_250m - isc_left_250m )
    alpha_lat = ( lat_right_250m - lat_left_250m ) / ( isc_right_250m - isc_left_250m )

    lat = lat_left_250m + alpha_lat * ( isc_250m - isc_left_250m )
    lon = lon_left_250m + alpha_lon * ( isc_250m - isc_left_250m )

    return lat, lon

def get_y_pos_1km_to_250m ( isc_1km, itk_p1_1km, itk_p2_1km,
                            lat_1km, lon_1km,
                            itk_250m ) :
    """
    return the position of 250m pixel defined by its along-track index between 2 successive 1km pixels on a same cross-track line
    @warning If the 2 pixels are crossing the changing date meridian, the negative longitude will be returned
    with an increment of 360 to manage correctly the interpolation
    @param isc_1km index of the cross-track line in the 1km grid
    @param itk_p1_1km index of the first along-track pixel in the 1km grid
    @param itk_p2_1km index of the second along-track pixel in the 1km grid
    @param lat_1km latitudes dataset at 1km resolution
    @param lon_1km longitudes dataset at 1km resolution
    @param itk_250m index of the along-track 250m along-track pixel to interpolate
    @return (lat,lon) of the pixel at [itk_250m, isc_1km], interpolated between the 2 successive 1km pixels p1 and p2
    """
    # make sure P1 and P1 are 2 successive pixels along-track
    if abs ( itk_p1_1km - itk_p2_1km ) != 1 :
        raise ValueError ( "P1 and P2 are on the same cross-track line and must be successive" )

    # lat, lon of the 1km bounding pixels
    p1_1km_lat = lat_1km [ itk_p1_1km, isc_1km ]
    p1_1km_lon = lon_1km [ itk_p1_1km, isc_1km ]
    p2_1km_lat = lat_1km [ itk_p2_1km, isc_1km ]
    p2_1km_lon = lon_1km [ itk_p2_1km, isc_1km ]

    # check for change date meridian particular case

    # change date meridian case sentinel
    if abs ( p1_1km_lon - p2_1km_lon ) > 180. :
        if p1_1km_lon < 0. :
            p1_1km_lon += 360.
        elif p2_1km_lon < 0. :
            p2_1km_lon += 360.

    # coordinates of p1, p2 in the 250m grid
    itk_p1_250m = itk_1km_to_250m ( itk_p1_1km )
    itk_p2_250m = itk_1km_to_250m ( itk_p2_1km )

    #print "itk_p1_250m=%f itk_p2_250m=%f"%( itk_p1_250m, itk_p2_250m )

    # linear interpolation on the position
    alpha_lon = ( p2_1km_lon - p1_1km_lon ) / ( itk_p2_250m - itk_p1_250m )
    alpha_lat = ( p2_1km_lat - p1_1km_lat ) / ( itk_p2_250m - itk_p1_250m )
    lon = p1_1km_lon + alpha_lon * ( itk_250m - itk_p1_250m )
    lat = p1_1km_lat + alpha_lat * ( itk_250m - itk_p1_250m )

    # in case of change date crossing, turn values > 180. to negative
    if lon > 180. :
        lon -= 360.
    elif lon < -180. :
        lon += 360.

    return lat, lon

def get_1km_pix_to_plot ( itk_250m, isc_250m, lat_1km, lon_1km ) :
    """
    return the arrays of (lat,lon) and of [i,j] labels of the 1km pixelsto draw arounf the (itk_250m, isc_250m) pixel
    """
    v_itk_1km = [ ]
    v_isc_1km = [ ]
    sc_250m_width = 40
    itk_1km = round ( itk_250m_to_1km ( itk_250m ) )
    isc_1km = max ( 0, min ( sz_sc_1km - 1, round ( isc_250m_to_1km ( isc_250m ) ) ) )

    if ( itk_250m % sc_250m_width ) <= 2 : # start of scan
        #itk_1km = round ( itk_1km )
        if isc_250m <= 2 :
            v_itk_1km = [   itk_1km,     itk_1km,     itk_1km,
                            itk_1km + 1, itk_1km + 1, itk_1km + 1,
                            itk_1km + 2, itk_1km + 2, itk_1km + 2 ]
            v_isc_1km = [  isc_1km,     isc_1km + 1, isc_1km + 2,
                            isc_1km,     isc_1km + 1, isc_1km + 2,
                            isc_1km,     isc_1km + 1, isc_1km + 2 ]
        elif isc_250m >= 5410 :
            v_itk_1km = [  itk_1km,     itk_1km    , itk_1km,
                            itk_1km + 1, itk_1km + 1, itk_1km + 1,
                            itk_1km + 2, itk_1km + 2, itk_1km + 2 ]
            v_isc_1km = [  isc_1km,     isc_1km - 1, isc_1km - 2,
                            isc_1km,     isc_1km - 1, isc_1km - 2,
                            isc_1km,     isc_1km - 1, isc_1km - 2 ]
        else :
            v_itk_1km = [  itk_1km,     itk_1km    , itk_1km,
                            itk_1km + 1, itk_1km + 1, itk_1km + 1,
                            itk_1km + 2, itk_1km + 2, itk_1km + 2 ]
            v_isc_1km = [  isc_1km,     isc_1km - 1, isc_1km + 1,
                            isc_1km,     isc_1km - 1, isc_1km + 1,
                            isc_1km,     isc_1km - 1, isc_1km + 1 ]
    elif ( itk_250m % sc_250m_width ) >= 37 : # end of scan
        #itk_1km = round ( itk_1km )
        if isc_250m <= 2 :
            v_itk_1km = [  itk_1km,     itk_1km    , itk_1km,
                           itk_1km - 1, itk_1km - 1, itk_1km - 1,
                           itk_1km - 2, itk_1km - 2, itk_1km - 2 ]
            v_isc_1km = [  isc_1km,     isc_1km + 1, isc_1km + 2,
                           isc_1km,     isc_1km + 1, isc_1km + 2,
                           isc_1km,     isc_1km + 1, isc_1km + 2 ]
        elif isc_250m >= 5410 :
            v_itk_1km = [  itk_1km,     itk_1km    , itk_1km,
                            itk_1km - 1, itk_1km - 1, itk_1km - 1,
                            itk_1km - 2, itk_1km - 2, itk_1km - 2 ]
            v_isc_1km = [  isc_1km,     isc_1km - 1, isc_1km - 2,
                           isc_1km,     isc_1km - 1, isc_1km - 2,
                           isc_1km,     isc_1km - 1, isc_1km - 2 ]
        else :
            v_itk_1km = [  itk_1km,     itk_1km    , itk_1km,
                            itk_1km - 1, itk_1km - 1, itk_1km - 1,
                            itk_1km - 2, itk_1km - 2, itk_1km - 2 ]
            v_isc_1km = [  isc_1km,     isc_1km - 1, isc_1km + 1,
                            isc_1km,     isc_1km - 1, isc_1km + 1,
                            isc_1km,     isc_1km - 1, isc_1km + 1 ]
    else : # middle of scan
        #itk_1km = round ( itk_1km )
        if isc_1km == 0 :
            v_itk_1km = [  itk_1km,     itk_1km    , itk_1km,
                            itk_1km - 1, itk_1km - 1, itk_1km - 1,
                            itk_1km + 1, itk_1km + 1, itk_1km + 1 ]
            v_isc_1km = [  isc_1km,     isc_1km + 1, isc_1km + 2,
                            isc_1km,     isc_1km + 1, isc_1km + 2,
                            isc_1km,     isc_1km + 1, isc_1km + 2 ]
        elif isc_1km == 1353 :
            v_itk_1km = [  itk_1km,     itk_1km    , itk_1km,
                            itk_1km - 1, itk_1km - 1, itk_1km - 1,
                            itk_1km + 1, itk_1km + 1, itk_1km + 1 ]
            v_isc_1km = [  isc_1km,     isc_1km - 1, isc_1km - 2,
                            isc_1km,     isc_1km - 1, isc_1km - 2,
                            isc_1km,     isc_1km - 1, isc_1km - 2 ]
        else :
            v_itk_1km = [  itk_1km,     itk_1km    , itk_1km,
                            itk_1km - 1, itk_1km - 1, itk_1km - 1,
                            itk_1km + 1, itk_1km + 1, itk_1km + 1 ]
            v_isc_1km = [  isc_1km,     isc_1km - 1, isc_1km + 1,
                            isc_1km,     isc_1km - 1, isc_1km + 1,
                            isc_1km,     isc_1km - 1, isc_1km + 1 ]
    npix = len(v_itk_1km)
    # build the labels array
    v_label = [ "[%d,%d]"%(v_itk_1km[i],v_isc_1km[i]) for i in xrange ( npix ) ]

    return v_itk_1km, v_isc_1km, v_label

def get_bounds_1km_to_250m( itk_1km, isc_1km ) :
    """
    return the 250m pixel indexes limits in the 1km pixel [ itk_1km, isc_1km ] footprint
    """

    # set the (track,scan) indexes of the 1km pixel in the 1km grid
    itk_250m = itk_1km_to_250m ( itk_1km )
    isc_250m = isc_1km_to_250m ( isc_1km )

    # set the 250m indexes of pixels to interpolate along track
    itk_250m_min    =  int ( itk_250m - 1.5 )
    itk_250m_max    =  int ( itk_250m + 1.5 )

    # general case : 2 interpolations done along scan : [isc-1, isc] then [isc, isc+1]
    isc_250m_min  = isc_250m - 2
    isc_250m_max  = isc_250m + 2
    if   ( isc_1km == 0 ) :
        isc_250m_min = 0
    elif ( isc_1km >= sz_sc_1km - 1 ) :
        isc_250m_max  = isc_250m + 3

    #print itk_1km, itk_250m_min, itk_250m_max
    #print isc_1km, isc_250m_min, isc_250m_max

    return itk_250m_min, itk_250m_max, isc_250m_min, isc_250m_max

##############  MAIN  ###############
def main():
    """
    Main script
    """
    if len ( sys.argv ) != 4 :
        raise ValueError ( "Invalid number of arguments.\nUsage : python modis_1km_to_250m_geolocation.py <infile> <itrack_250m> <iscan_250m>" )

    infile     =       sys.argv [ 1 ]
    itk_250m   = int ( sys.argv [ 2 ] )
    isc_250m   = int ( sys.argv [ 3 ] )
    #print itk_250m, isc_250m
    #print itk_250m_to_1km ( itk_250m )
    #itk_1km = max ( 0, math.floor ( itk_250m_to_1km ( itk_250m ) ) + 1 )
    #isc_1km = math.floor ( isc_250m_to_1km ( isc_250m ) )
    #print itk_1km, isc_1km

    # -----------------------------------
    # --- Load the 1km (lat,lon) data ---
    # -----------------------------------
    hdf  = SD.SD ( infile )

    # read the latitudes of the points
    sds  = hdf.select ( "Latitude" )
    lat_1km = sds.get ( )
    sds.endaccess()

    # read the longitudes of the points
    sds  = hdf.select ( "Longitude" )
    lon_1km = sds.get ( )
    sds.endaccess()

    hdf.end()

    # ---------------------------------------------
    # --- compute the (lat,lon) of a 250m pixel ---
    # ---------------------------------------------
    lat, lon = get_250m_pix_pos ( itk_250m, isc_250m, lat_1km, lon_1km )
    print "%f\t%f"%(lat,lon)

    if __DEBUG__ :
        # ----------------------------------
        # --- set the 1km pixels to plot ---
        # ----------------------------------
        sz_tk_1km, sz_sc_1km = lat_1km.shape
        #if itk_1km >= sz_tk_1km :
            #itk_1km = sz_tk_1km - 1
        #if isc_1km >= sz_sc_1km :
            #isc_1km = sz_sc_1km - 1
        v_itk_1km, v_isc_1km, v_label = get_1km_pix_to_plot ( itk_250m, isc_250m, lat_1km, lon_1km )

        # set the lat,lon of those points
        v_lat_1km = lat_1km [ v_itk_1km, v_isc_1km ]
        v_lon_1km = lon_1km [ v_itk_1km, v_isc_1km ]

        # ----------------------------------
        # --- plot 1km pixels
        # ----------------------------------
        # - init the plot where will be drawn the pixels
        fig = plt.figure ()
        ax  = plt.subplot (111)

        plt.xlabel ('Longitude' )
        plt.ylabel ('Latitude' )
        plt.axis   ( 'equal' )
        ax.get_xaxis().set_major_formatter( ticker.FormatStrFormatter('%.2f') )
        ax.get_yaxis().set_major_formatter( ticker.FormatStrFormatter('%.2f') )

        plt.axis   ( [ v_lon_1km.min(), v_lon_1km.max(), v_lat_1km.min(), v_lat_1km.max() ] )

        #print "v_lon_1km "+str(v_lon_1km)
        #print "v_lat_1km "+str(v_lat_1km)

        # - plot the 250m pixel
        plt.scatter ( [ lon ], [ lat ], color='m', marker='x', linewidth = 1 )
        ax.text     ( lon, lat, "(%d,%d)"%( itk_250m, isc_250m ), color = 'm' )

        # - plot the 1km pixels
        npix = len ( v_lat_1km )
        for i in xrange ( npix ) :
            _lat   = v_lat_1km [ i ]
            _lon   = v_lon_1km [ i ]
            label = v_label   [ i ]
            sz_marker = 1
            if i == 0 :
                # set a biggest marker for pixel to interpolate
                sz_marker = 5
            plt.scatter ( [ _lon ], [ _lat ], color='r', marker='o', linewidth = sz_marker, label = "_" )
            ax.text     ( _lon, _lat, label )

        # ------------------------------------------
        # --- Compute and plot inner 250m pixels ---
        # ------------------------------------------
        itk_1km = round ( itk_250m_to_1km ( itk_250m ) )
        isc_1km = min ( sz_sc_1km - 1, round ( isc_250m_to_1km ( isc_250m ) ) )

        itk_250m_min, itk_250m_max, isc_250m_min, isc_250m_max = get_bounds_1km_to_250m( itk_1km, isc_1km )
        v_lat_250m = []
        v_lon_250m = []
        v_itk_250m = []
        v_isc_250m = []

        for itk_250m in xrange ( itk_250m_min, itk_250m_max + 1 ) :
            for isc_250m in xrange ( isc_250m_min, isc_250m_max + 1 ) :
                v_itk_250m.append ( itk_250m )
                v_isc_250m.append ( isc_250m )
                # intert/extrapolate the lat, lon of this 250m pixel using the 1km grid
                lat, lon = get_250m_pix_pos ( itk_250m, isc_250m, lat_1km, lon_1km )
                v_lat_250m.append ( lat )
                v_lon_250m.append ( lon )
        v_itk_250m = np.array ( v_itk_250m )
        v_isc_250m = np.array ( v_isc_250m )
        v_lat_250m = np.array ( v_lat_250m )
        v_lon_250m = np.array ( v_lon_250m )

        plt.scatter ( v_lon_250m [ v_lon_250m != -9999. ].flat, v_lat_250m [ v_lat_250m != -9999. ].flat, \
                    color='b', marker='+', linewidth = 1, label = "_" )

        plt.show ()

if __name__ == "__main__":
  main()
