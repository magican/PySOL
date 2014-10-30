# -*- coding: utf-8 -*-

"""
Description :
    Give the geolocation of a MODIS 1km pixel specified by its indices [i_1km_along_track, i_1km_cross_track] for products where the (latitude, longitude) datasets are given at 5km resolution
    For the interpolation of the 1km positions, it takes care of borders of scan. In this cases, the position will be extrapolated.
    It also deals with the change date meridian crossing.

Limitation :
    This script works well with MYD05 and MYD06 products but needs a minor adaptations for fully running on MY021KM products.
    Due to the difference of along-scan dimension size ( 271 pixels for MYD021KM and 270 for others ), the pixels with a 1km scan index
    in range [1349,1352] are always extrapolated using the 5km scan pixels [268,269].whereas they can be interpolated between [269,270] in
    the MYD021KM case.

Usage :
    python modis_5km_to_1km_geolocation.py <infile> <i_1km_along_track> <i_1km_cross_track>

    It will print out the ( lat, lon ) geolocation of the given pixel.

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
v0.1.0 : 2009/11/25
 - creation
"""

import math
import sys

import warnings
warnings.filterwarnings("ignore")

import numpy as np

__DEBUG__ = False
if __DEBUG__ is True :
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

# --- 5km -> 1km ---
sz_sc_5km = 271 # number of 5km pixels along scan
# sz_sc_5km = 270

def itk_5km_to_1km ( i_tk_5km ) :
    """
    return the 1km grid index along track of a 5km pixel
    """
    return 2 + 5 * i_tk_5km

def itk_1km_to_5km ( i_tk_1km ) :
    """
    return the 5km grid index along track of a 1km pixel
    """
    return ( i_tk_1km - 2. ) / 5.

def isc_5km_to_1km ( i_sc_5km ) :
    """
    return the 1km grid index cross track of a 5km pixel
    """
    return 2. + 5. * i_sc_5km

def isc_1km_to_5km ( i_sc_1km ) :
    """
    return the 5km grid index cross track of a 1km pixel
    """
    return  ( i_sc_1km - 2. ) / 5.

def get_1km_pix_pos ( itk_1km, isc_1km, lat_5km, lon_5km ) :
    """
    return the (lat,lon) of a 1km pixel specified with its indexes, when the geolocation datasets are given at 5km resolution
    @param itk_1km grid index of the 1km pixel along-track
    @param isc_1km grid index of the 1km pixel cross-track
    @param lat_5km latitudes dataset at 5km resolution
    @param lon_5km longitudes dataset at 5km resolution
    @return the ( lat, lon ) of the 1km pixel  [ itk_1km, isc_1km ]
    """
    # check 1km indexes validity
    sz_tk_5km = lat_5km.shape[0] + 1
    sz_tk_1km = 5 * sz_tk_5km
    sz_sc_1km = ( 5 * sz_sc_5km ) + 6
    if ( isc_1km < 0 ) or ( isc_1km > sz_sc_1km - 1 ) :
        raise ValueError ( "Invalid scan index %d. Must be in range [%d,%d]"%(isc_1km,0,sz_sc_1km-1) )
    if ( itk_1km < 0 ) or ( itk_1km > sz_tk_1km - 1 ) :
        raise ValueError ( "Invalid track index %d. Must be in range [%d,%d]"%(itk_1km,0,sz_tk_1km-1) )

    # --- set the bounding 5km pixels to take for interpolation
    # set the (track,scan) indexes of the 5km pixel in the 1km grid
    itk_5km = itk_1km_to_5km ( itk_1km )
    isc_5km = isc_1km_to_5km ( isc_1km )

    #print "i_1km=[%.2f, %.2f] -> i_5km=[%.2f, %.2f]"%(itk_1km, isc_1km, itk_5km, isc_5km)

    # the width of one scan, in number of pixels, at 1km resolution
    w_scan_1km = 10.
    # - extra/interpolation 5km pixels along track
    if ( itk_1km % w_scan_1km ) <= 2 : # extrapolate start of scan
        itk_top_5km    = math.ceil ( itk_5km )
        itk_bottom_5km = itk_top_5km + 1
    elif ( itk_1km % w_scan_1km ) >= 7 : # extrapolate end of scan
        itk_top_5km    = math.floor ( itk_5km )
        itk_bottom_5km = itk_top_5km - 1
    else : # general case : middle of scan
        itk_top_5km    = math.floor ( itk_5km )
        itk_bottom_5km = itk_top_5km + 1
    # - extra/interpolation 5km pixels cross track
    if ( isc_1km <= 2 ) : # extrapolate start of scan line
        isc_left_5km  = 0
        isc_right_5km = 1
    elif ( isc_5km >= ( sz_sc_5km - 1 )  ) : # extrapolate end of scan line
        isc_left_5km  = sz_sc_5km - 2
        isc_right_5km = sz_sc_5km - 1
    else : # general case : interpolation
        isc_left_5km  = math.floor ( isc_5km )
        isc_right_5km = isc_left_5km + 1

    #print "itk_top_5km=%d itk_bottom_5km=%d isc_left_5km=%d isc_right_5km=%d"%(itk_top_5km, itk_bottom_5km, isc_left_5km, isc_right_5km)

    # --- set the 5km track lines position ; left border ---
    lat_left_1km, lon_left_1km   = get_y_pos_5km_to_1km (
                isc_left_5km, itk_top_5km, itk_bottom_5km,
                lat_5km, lon_5km,
                itk_1km )
    # --- set the 5km track lines position ; right border ---
    lat_right_1km, lon_right_1km = get_y_pos_5km_to_1km (
                isc_right_5km, itk_top_5km, itk_bottom_5km,
                lat_5km, lon_5km,
                itk_1km )

    #print "left_1km=[%f,%f] right_1km=[%f,%f]"%(lat_left_1km, lon_left_1km, lat_right_1km, lon_right_1km)

    # check for change date meridian case
    if  abs ( lon_right_1km  - lon_left_1km  ) > 180. :
        # all negative longitudes will be incremented of 360 before interpolation
        if lon_left_1km < 0. :
            lon_left_1km += 360.
        if lon_right_1km < 0. :
            lon_right_1km += 360.

    # for each track line position, interpolate along scan to retrieve the 1km geolocation
    lat, lon = get_x_pos_5km_to_1km ( lat_left_1km,  lon_left_1km,  isc_left_5km,
                                     lat_right_1km, lon_right_1km, isc_right_5km,
                                     isc_1km )
    #print "geolocation = [%f, %f]"%(lat,lon)

    # in case of change date crossing, turn values > 180. to negative < -180 to positive
    if lon > 180. :
        lon = lon - 360.
    elif lon < -180. :
        lon = lon + 360.
    return lat, lon

def get_x_pos_5km_to_1km ( lat_left_1km,  lon_left_1km,  isc_left_5km,
                           lat_right_1km, lon_right_1km, isc_right_5km,
                           isc_1km
                        ) :
    """
    retrieve the position of a 1km pixel set by its cross-track index and the position of 2 along-track lines set by 2 sucessive 5km
    pixels
    PRECONDITION : change date meridian case must have already been treated, and so lon_left_1km or lon_right_1km can have values > 180.
    @param lat_left_1km latitude of the 1km along track pixel that is the intersection between cross-track line at index isc_1km and
    the along-track line that joins the 2 sucessive bounding 5km pixels on the left of the point to interpolate
    @param lon_left_1km longitude of the 1km along track pixel that is the intersection between cross-track line at index isc_1km and
    the along-track line that joins the 2 sucessive bounding 5km pixels on the left of the point to interpolate
    @param isc_left_5km cross-track index at 5km of the left border
    @param lat_right_1km latitude of the 1km along track pixel that is the intersection between cross-track line at index isc_1km and
    the along-track line that joins the 2 sucessive bounding 5km pixels on the right of the point to interpolate
    @param lon_right_1km longitude of the 1km along track pixel that is the intersection between cross-track line at index isc_1km and
    the along-track line that joins the 2 sucessive bounding 5km pixels on the right of the point to interpolate
    @param isc_right_5km cross-track index at 5km of the right border
    @return (lat,lon) of the pixel at [itk_1km, isc_5km], interpolated between the 2 successive 5km pixels p1 and p2
    """
    # make sure P1 and P1 are 2 successive pixels along-track
    if abs ( isc_left_5km - isc_right_5km ) != 1 :
        raise ValueError ( "The 2 borders are on the same along-track line and must be successive" )

    # coordinates of left and right border in the 1km grid
    isc_left_1km  = isc_5km_to_1km ( isc_left_5km )
    isc_right_1km = isc_5km_to_1km ( isc_right_5km )

    #print "isc_left_5km=%d isc_right_5km=%d"%(isc_left_5km, isc_right_5km)
    #print "isc_left_1km=%d isc_right_1km=%d"%(isc_left_1km, isc_right_1km)
    #print "isc_1km_min=%d isc_1km_max=%d"%(isc_1km_min,isc_1km_max)

    # linear interpolation on the position
    alpha_lon = ( lon_right_1km - lon_left_1km ) / ( isc_right_1km - isc_left_1km )
    alpha_lat = ( lat_right_1km - lat_left_1km ) / ( isc_right_1km - isc_left_1km )

    lat = lat_left_1km + alpha_lat * ( isc_1km - isc_left_1km )
    lon = lon_left_1km + alpha_lon * ( isc_1km - isc_left_1km )

    return lat, lon

def get_y_pos_5km_to_1km (  isc_5km, itk_p1_5km, itk_p2_5km,
                            lat_5km, lon_5km,
                            itk_1km ) :
    """
    return the position of 1km pixel defined by its along-track index between 2 successive 5km pixels on a same cross-track line
    @warning If the 2 pixels are crossing the changing date meridian, the negative longitude will be returned
    with an increment of 360 to manage correctly the interpolation
    @param isc_5km index of the cross-track line in the 5km grid
    @param itk_p1_5km index of the first along-track pixel in the 5km grid
    @param itk_p2_5km index of the second along-track pixel in the 5km grid
    @param lat_5km latitudes dataset at 5km resolution
    @param lon_5km longitudes dataset at 5km resolution
    @param itk_1km index of the along-track 1km along-track pixel to interpolate
    @return (lat,lon) of the pixel at [itk_1km, isc_5km], interpolated between the 2 successive 5km pixels p1 and p2
    """
    # make sure P1 and P1 are 2 successive pixels along-track
    if abs ( itk_p1_5km - itk_p2_5km ) != 1 :
        raise ValueError ( "P1 and P2 are on the same cross-track line and must be successive" )

    # lat, lon of the 5km bounding pixels
    p1_5km_lat = lat_5km [ itk_p1_5km, isc_5km ]
    p1_5km_lon = lon_5km [ itk_p1_5km, isc_5km ]
    p2_5km_lat = lat_5km [ itk_p2_5km, isc_5km ]
    p2_5km_lon = lon_5km [ itk_p2_5km, isc_5km ]

    # check for change date meridian case
    if  abs ( p1_5km_lon  - p2_5km_lon  ) > 180. :
        # all negative longitudes will be incremented of 360 before interpolation
        if p1_5km_lon < 0. :
            p1_5km_lon += 360.
        if p2_5km_lon < 0. :
            p2_5km_lon += 360.

    # coordinates of p1, p2 in the 1km grid
    itk_p1_1km = itk_5km_to_1km ( itk_p1_5km )
    itk_p2_1km = itk_5km_to_1km ( itk_p2_5km )

    #print "itk_p1_1km=%f itk_p2_1km=%f"%( itk_p1_1km, itk_p2_1km )

    # linear interpolation on the position
    alpha_lon = ( p2_5km_lon - p1_5km_lon ) / ( itk_p2_1km - itk_p1_1km )
    alpha_lat = ( p2_5km_lat - p1_5km_lat ) / ( itk_p2_1km - itk_p1_1km )
    lon = p1_5km_lon + alpha_lon * ( itk_1km - itk_p1_1km )
    lat = p1_5km_lat + alpha_lat * ( itk_1km - itk_p1_1km )

    if lon > 180. :
        lon = lon - 360.
    elif lon < -180. :
        lon = lon + 360.

    return lat, lon

def get_5km_pix_to_plot ( itk_5km, isc_5km, lat_5km, lon_5km ) :
    """
    return the arrays of (lat,lon) and of [i,j] labels of the 5km pixels that bound the one at (itk_5km, isc_5km).
    The first element of these arrays represent the (itk_5km, isc_5km) pixel
    """
    v_itk_5km = [ ]
    v_isc_5km = [ ]
    sc_width = 2 # width of the scan along track
    if ( itk_5km % sc_width ) == 0 : # start of scan
        v_itk_5km = [  itk_5km,     itk_5km    , itk_5km,
                       itk_5km + 1, itk_5km + 1, itk_5km + 1 ]
        if isc_5km == 0 :
            v_isc_5km = [  isc_5km,     isc_5km + 1, isc_5km + 2,
                           isc_5km,     isc_5km + 1, isc_5km + 2 ]
        elif isc_5km == 269 :
            v_isc_5km = [  isc_5km,     isc_5km - 1, isc_5km - 2,
                            isc_5km,     isc_5km - 1, isc_5km - 2 ]
        else :
            v_isc_5km = [   isc_5km,     isc_5km - 1, isc_5km + 1,
                            isc_5km,     isc_5km - 1, isc_5km + 1 ]
    elif ( itk_5km % sc_width ) == ( sc_width - 1 ) : # end of scan
        v_itk_5km = [  itk_5km,     itk_5km    , itk_5km,
                       itk_5km - 1, itk_5km - 1, itk_5km - 1 ]
        if isc_5km == 0 :
            v_isc_5km = [  isc_5km,     isc_5km + 1, isc_5km + 2,
                           isc_5km,     isc_5km + 1, isc_5km + 2 ]
        elif isc_5km == 269 :
            v_isc_5km = [  isc_5km,     isc_5km - 1, isc_5km - 2,
                           isc_5km,     isc_5km - 1, isc_5km - 2 ]
        else :
            v_isc_5km = [  isc_5km,     isc_5km - 1, isc_5km + 1,
                           isc_5km,     isc_5km - 1, isc_5km + 1 ]

    npix = len( v_itk_5km )
    # build the labels array
    v_label = [ "[%d,%d]"%(v_itk_5km[i],v_isc_5km[i]) for i in xrange ( npix ) ]

    return v_itk_5km, v_isc_5km, v_label

def get_bounds_5km_to_1km( itk_5km, isc_5km ) :
    """
    return the 1km pixel indexes limits in the 5km pixel [ itk_5km, isc_5km ] footprint
    """

    # set the (track,scan) indexes of the 5km pixel in the 5km grid
    itk_1km = itk_5km_to_1km ( itk_5km )
    isc_1km = isc_5km_to_1km ( isc_5km )

    # set the 1km indexes of pixels to interpolate along track
    itk_1km_min    =  itk_1km - 2
    itk_1km_max    =  itk_1km + 2

    # general case : 2 interpolations done along scan : [isc-1, isc] then [isc, isc+1]
    isc_1km_min  = isc_1km - 2
    isc_1km_max  = isc_1km + 2
    # if last 5km pixel along scan, only 4 1km pixels in the 5km footprint in this direction
    if ( isc_5km == sz_sc_5km - 1 ) :
        isc_1km_max  = isc_1km + 6

    return itk_1km_min, itk_1km_max, isc_1km_min, isc_1km_max

##############  MAIN  ###############
def main():
    """
    Main script
    """
    if len ( sys.argv ) != 4 :
        raise ValueError ( "Invalid number of arguments.\nUsage : python modis_5km_to_1km_geolocation.py <infile> <i_track_1km> <iscan_1km>" )

    infile  =       sys.argv [ 1 ]
    itk_1km = int ( sys.argv [ 2 ] )
    isc_1km = int ( sys.argv [ 3 ] )

    # -----------------------------------
    # --- Load the 5km (lat,lon) data ---
    # -----------------------------------
    hdf  = SD.SD ( infile )

    # read the latitudes of the points
    sds  = hdf.select ( "Latitude" )
    lat_5km = sds.get ( )
    sds.endaccess()

    # read the longitudes of the points
    sds  = hdf.select ( "Longitude" )
    lon_5km = sds.get ( )
    sds.endaccess()

    hdf.end()

    # ---------------------------------------------
    # --- compute the (lat,lon) of a 1km pixel ---
    # ---------------------------------------------
    lat, lon = get_1km_pix_pos ( itk_1km, isc_1km, lat_5km, lon_5km )
    print "%f\t%f"%( lat, lon )

    if __DEBUG__ :
        # ----------------------------------
        # --- set the 5km pixels to plot ---
        # ----------------------------------
        itk_5km = itk_1km_to_5km ( itk_1km )
        isc_5km = isc_1km_to_5km ( isc_1km )

        itk_5km = max ( 0, round ( itk_1km_to_5km ( itk_1km ) ) )
        isc_5km = max ( 0, min ( round ( isc_1km_to_5km ( isc_1km ) ), sz_sc_5km - 1 ) )

        sz_tk_5km = lat_5km[0].shape
        if itk_5km >= sz_tk_5km :
            itk_5km = sz_tk_5km - 1
        if isc_5km >= sz_sc_5km :
            isc_5km = sz_sc_5km - 1
        v_itk_5km, v_isc_5km, v_label = get_5km_pix_to_plot ( itk_5km, isc_5km, lat_5km, lon_5km )

        # set the lat,lon of those points
        v_lat_5km = lat_5km [ v_itk_5km, v_isc_5km ]
        v_lon_5km = lon_5km [ v_itk_5km, v_isc_5km ]

        # ----------------------------------
        # --- plot 5km pixels
        # ----------------------------------
        # - init the plot where will be drawn the pixels
        fig = plt.figure ()
        ax  = plt.subplot (111)

        plt.xlabel ( 'Longitude' )
        plt.ylabel ( 'Latitude' )
        plt.axis   ( 'equal' )
        ax.get_xaxis().set_major_formatter( ticker.FormatStrFormatter ( '%.2f' ) )
        ax.get_yaxis().set_major_formatter( ticker.FormatStrFormatter ( '%.2f' ) )

        plt.axis   ( [ v_lon_5km.min(), v_lon_5km.max(), v_lat_5km.min(), v_lat_5km.max() ] )

        #print "v_lon_5km "+str(v_lon_5km)
        #print "v_lat_5km "+str(v_lat_5km)

        # - plot the 1km pixel
        plt.scatter ( [ lon ], [ lat ], color='m', marker='x', linewidth = 1 )
        ax.text     ( lon, lat, "(%d,%d)"%( itk_1km, isc_1km ), color = 'm' )

        # - plot the 5km pixels
        npix = len ( v_lat_5km )
        for i in xrange ( npix ) :
            _lat   = v_lat_5km [ i ]
            _lon   = v_lon_5km [ i ]
            label = v_label   [ i ]
            sz_marker = 1
            if i == 0 :
                # set a biggest marker for pixel to interpolate
                sz_marker = 5
            plt.scatter ( [ _lon ], [ _lat ], color='r', marker='o', linewidth = sz_marker, label = "_" )
            ax.text     ( _lon, _lat, label )

        # --------------------------------------------------------------------
        # --- Compute and plot the in 1km pixels inside the 5km footprint  ---
        # --------------------------------------------------------------------
        itk_1km_min, itk_1km_max, isc_1km_min, isc_1km_max = get_bounds_5km_to_1km( itk_5km, isc_5km )
        v_lat_1km = []
        v_lon_1km = []
        v_itk_1km = []
        v_isc_1km = []

        for itk_1km in xrange ( itk_1km_min, itk_1km_max + 1 ) :
            for isc_1km in xrange ( isc_1km_min, isc_1km_max + 1 ) :
                v_itk_1km.append ( itk_1km )
                v_isc_1km.append ( isc_1km )
                # intert/extrapolate the lat, lon of this 1km pixel using the 5km grid
                lat, lon = get_1km_pix_pos ( itk_1km, isc_1km, lat_5km, lon_5km )
                v_lat_1km.append ( lat )
                v_lon_1km.append ( lon )
        v_itk_1km = np.array ( v_itk_1km )
        v_isc_1km = np.array ( v_isc_1km )
        v_lat_1km = np.array ( v_lat_1km )
        v_lon_1km = np.array ( v_lon_1km )

        #print "v_lon_1km "  + str ( v_lon_1km )
        #print "v_lat_1km "  + str ( v_lat_1km )
        #print "v_itk_1km " + str ( v_itk_1km )
        #print "v_isc_1km " + str ( v_isc_1km )

        plt.scatter ( v_lon_1km [ v_lon_1km != -9999. ].flat, v_lat_1km [ v_lat_1km != -9999. ].flat, \
                    color='b', marker='+', linewidth = 1, label = "_" )

        plt.show ()

if __name__ == "__main__":
  main()
