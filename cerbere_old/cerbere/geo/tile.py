'''
Created on 28 mars 2012

@author: jfpiolle
'''

try:
    from shapely.wkt import loads as OGRGeometry
except ImportError:
    # Fallback to GeoDjango
    try:
        from django.contrib.gis.gdal import OGRGeometry
    except ImportError:
        from django.contrib.gis.geos import GEOSGeometry as OGRGeometry

class Tile(object):
    '''
    classdocs
    '''


    def __init__(self, offsets=None,timerange=None,bbox=None,meridian180=False):
        '''
        Constructor
        @param offsets: min and max offsets of the data in product (up to 3 dimensions)
        @param timerange: tuple(start date, end date) of the time coverage of the data
        @param meridian180: flag for meridian 180 crossing (True if the tile is crossing)
        '''
        if offsets != None:
            self.xMinOffset = None
            self.xMaxOffset = None
            self.yMinOffset = None
            self.yMaxOffset = None
            self.zMinOffset = None
            self.zMaxOffset = None
            if len(offsets) == 2:
                x0,x1 = offsets
                self.xMinOffset = x0
                self.xMaxOffset = x1
            elif len(offsets) == 4:
                x0,x1,y0,y1 = offsets
                self.xMinOffset = x0
                self.xMaxOffset = x1
                self.yMinOffset = y0
                self.yMaxOffset = y1
        if timerange != None:
            self.rangeBeginningDateTime = timerange[0]
            self.rangeEndingDateTime =  timerange[1]
        if bbox != None:
            coordstr = ""
            for i,pt in enumerate(bbox):
                if i==0:
                    coordstr = "%0.2f %0.2f" % (pt.lon,pt.lat)
                else:
                    coordstr += ", %0.2f %0.2f" % (pt.lon,pt.lat)
            polystr = 'POLYGON((%s))' % coordstr
            self.bbox = OGRGeometry( polystr )
        self.meridian180 = meridian180
