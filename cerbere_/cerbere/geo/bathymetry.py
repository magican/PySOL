# -*- coding: utf-8 -*-
"""
.. module::cerbere.geo.bathymetry

Class to handle bathymetry information and distance to coast.

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import os
import math

from netCDF4 import Dataset

RES1MIN = "GridOne.nc"
RES1MIN_DIST = "gridone_dist2shore.nc"
RES30S = "GEBCO_08.nc"


class Bathymetry:
    """A class to handle bathymetry information and distance to coast."""

    def __init__(self, resolution=RES1MIN, filename=None, cache=False):
        if filename is None:
            filename = os.path.join(os.path.dirname(__file__),
                                    'resources',
                                    resolution)
        self._bathyfile = None
        self._bathyfile = Dataset(filename,
                                  'r',
                                  format='NETCDF3_CLASSIC')
        if cache:
            self.bathymetry = self._bathyfile.variables['z'][:]
        else:
            self.bathymetry = self._bathyfile.variables['z']
#             print self.bathymetry.shape
        self._resolution = resolution
        return

    def __del__(self):
        if self._bathyfile is not None:
            self._bathyfile.close()

    def get_distance_to_shore(self, lon, lat):
        """Return the distance to the coast, in meters"""
        if (self._resolution == RES1MIN_DIST):
            lon_in_min = math.floor(lon) * 60 + (lon - math.floor(lon)) * 60.
            lat_in_min = math.floor(lat) * 60 + (lat - math.floor(lat)) * 60.
            i = int(180. * 60 + math.floor(lon_in_min))\
                + 21601 * int(90. * 60. - math.floor(lat_in_min))
        else:
            raise Exception("Class must be initialized with a distance file"
                            "like RES1MIN_DIST")
        distance = self.bathymetry[i]
        return distance

    def get_depth(self, lon, lat):
        """Return the depth, in meters"""
        # open Land Mask
        if (self._resolution == RES1MIN):
            lon_in_min = math.floor(lon) * 60 + (lon - math.floor(lon)) * 60.
            lat_in_min = math.floor(lat) * 60 + (lat - math.floor(lat)) * 60.
            i = int(180. * 60 + math.floor(lon_in_min))\
                + 21601 * int(90. * 60. - math.floor(lat_in_min))
        elif (self._resolution == RES30S):
            lon_in_min = math.floor(lon) * 60 + (lon - math.floor(lon)) * 60.
            lat_in_min = math.floor(lat) * 60 + (lat - math.floor(lat)) * 60.
            i = int(180. * 60 * 2. + math.floor(lon_in_min * 2.))\
                + 43200 * int(90. * 60. * 2. - math.floor(lat_in_min * 2.))
        else:
            raise Exception("Class must be initialized with a bathymetry file"
                            "like RES1MIN or RES30S")
        bathy = self.bathymetry[i]
        return bathy


# =====================================================
# test
# =====================================================
if __name__ == "__main__":

    bathy1 = Bathymetry(RES1MIN)
    bathy2 = Bathymetry(RES30S)
    bathy3 = Bathymetry(RES1MIN_DIST)

    # test depth
    print bathy1.get_depth(-180., 90),\
        bathy2.get_depth(-180., 90), bathy3.get_distance_to_shore(-180., 90)
    print bathy1.get_depth(0., 90.),\
        bathy2.get_depth(0., 90.), bathy3.get_distance_to_shore(0., 90.)
    print bathy1.get_depth(90., 90.),\
        bathy2.get_depth(90., 90.), bathy3.get_distance_to_shore(90., 90.)
    print bathy1.get_depth(180., 90.),\
        bathy2.get_depth(180., 90.), bathy3.get_distance_to_shore(180., 90.)
    print bathy1.get_depth(-180., 89.9,),\
        bathy2.get_depth(180., 89.9,), bathy3.get_distance_to_shore(180., 89.9,)
    print bathy1.get_depth(45., 20.),\
        bathy2.get_depth(45., 20.), bathy3.get_distance_to_shore(45., 20.)
    print bathy1.get_depth(37.65167, -0.315),\
        bathy2.get_depth(37.65167, -0.315),\
        bathy3.get_distance_to_shore(37.65167, -0.315)
    print bathy1.get_depth(40.745, 1.456667),\
        bathy2.get_depth(40.745, 1.456667),\
        bathy3.get_distance_to_shore(40.745, 1.456667)
    print bathy1.get_depth(40.745, 2.456667),\
        bathy2.get_depth(40.745, 2.456667),\
        bathy3.get_distance_to_shore(40.745, 2.456667)
    print bathy1.get_depth(-6.197, 48.455),\
        bathy2.get_depth(-6.197, 48.455),\
        bathy3.get_distance_to_shore(-6.197, 48.455)