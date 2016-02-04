# encoding: utf-8
"""
cerbere.geo.spatialshape
=========================

Routines for handling the shape of a feature

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
import numpy

try:
    from shapely.wkt import loads as OGRGeometry
except ImportError:
    # Fallback to GeoDjango
    try:
        from django.contrib.gis.gdal import OGRGeometry
    except ImportError:
        from django.contrib.gis.geos import GEOSGeometry as OGRGeometry

from ..datamodel.swath import Swath
from ..datamodel.grid import Grid
from ..datamodel.trajectory import Trajectory


def get_swath_spatialshape(self):
    """Return the spatial shape of a swath feature"""
    lats = self.get_lat()
    lons = self.get_lon()
    lon0 = lons.min()
    lon1 = lons.max()
    lon2 = lon1
    lon3 = lon0
    lat0 = lats.max()
    lat1 = lat0
    lat2 = lats.min()
    lat3 = lat2
#     if isinstance(lats, numpy.ma.MaskedArray):
#         valid_rows = lats.any(1)
#         valid_cols = lats.any(0)
#     row0 = numpy.where(valid_rows)[0][0]
#     row1 = numpy.where(valid_rows)[0][-1]
#     cell0 = numpy.where(valid_cols)[0][0]
#     cell1 = numpy.where(valid_cols)[0][-1]
#     lat0, lat1 = numpy.ma.flatnotmasked_edges(lats[row0, cell0:cell1])
#     lon0, lon1 = numpy.ma.flatnotmasked_edges(lons[row0, cell0:cell1])
#     lat3, lat2 = numpy.ma.flatnotmasked_edges(lats[row1, cell0:cell1])
#     lon3, lon2 = numpy.ma.flatnotmasked_edges(lons[row1, cell0:cell1])
    polystr = 'POLYGON((%0.2f %0.2f, %0.2f %0.2f, %0.2f %0.2f, %0.2f %0.2f,%0.2f %0.2f))'\
            % (lon0, lat0,
               lon1, lat1,
               lon2, lat2,
               lon3, lat3,
               lon0, lat0)
    print polystr
    return OGRGeometry(polystr)


def get_spatialshape(feature):
    """Return the spatial shape of a feature"""
    if isinstance(feature, Swath) or isinstance(feature, Trajectory):
        return get_swath_spatialshape(feature)
    elif isinstance(feature, Grid):
        lonmin, latmin, lonmax, latmax = feature.get_bbox()
        polystr = 'POLYGON((%0.2f %0.2f, %0.2f %0.2f, %0.2f %0.2f, %0.2f %0.2f,%0.2f %0.2f))' %\
            (lonmin, latmin,
             lonmin, latmax,
             lonmax, latmax,
             lonmax, latmin,
             lonmin, latmin)
        return OGRGeometry(polystr)
    else:
        return None

# 
#     def get_spatialshape_trajectory(self):
#         """
#         returns the polygonal spatial shape of the feature as a geometry object. 
#         """
#         lats = self.getLats()
#         lons = self.getLons()
#         minLat = lats.min()
#         maxLat = lats.max()
#         minLon = lons.min()
#         maxLon = lons.max()
#         # crossing of dateline
# #        for i,lon in enumerate(lons[1:]):
# #            if abs(lon - lons[i-1]) > 180.:
# #                    minLon += 360.
# #                    break    
#         polystr = 'POLYGON((%0.2f %0.2f, %0.2f %0.2f, %0.2f %0.2f, %0.2f %0.2f,%0.2f %0.2f))' % \
#             (minLon,minLat,minLon,maxLat,maxLon,maxLat,maxLon,minLat,minLon,minLat)                
#         return OGRGeometry( polystr )
# 

