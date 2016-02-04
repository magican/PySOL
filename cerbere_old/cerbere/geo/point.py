'''
Created on 29 mars 2012

@author: jfpiolle
'''
import numpy

EARTH_MEAN_RADIUS = 6.37123E6  # meters

class Point(object):
    '''
    classdocs
    '''


    def __init__(self, lon, lat):
        '''
        Constructor
        '''
        self.lon = lon
        self.lat = lat
        
    def __str__(self):
        return "Lat : %02f   Lon : %02f" % (self.lat, self.lon)

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    @classmethod
    def get_distance(cls, point1, point2):
        """
        return the distance in meters between two points
        """
        return Point.get_distance_latlon(point1.lat,
                                         point1.lon,
                                         point2.lat,
                                         point2.lon)

    @classmethod
    def get_distance_latlon(cls, lats, lons, lat, lon):
        # convert to radian
        lats = lats * numpy.pi / 180.
        lons = lons * numpy.pi / 180.
        lat = numpy.radians(lat)
        lon = numpy.radians(lon)

        dlat = lats - lat
        dlon = lons - lon

        a = (numpy.sin(dlat / 2)) ** 2 \
                + numpy.cos(lat) * numpy.cos(lats) * (numpy.sin(dlon / 2)) ** 2
        c = 2 * numpy.arcsin(numpy.minimum(1, numpy.sqrt(a)))
        return c * EARTH_MEAN_RADIUS
