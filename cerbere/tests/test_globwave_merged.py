'''
Created on 2 aout 2012

@author: jfpiolle
'''

import logging, os

import datamodel.pointcollection
import mapper.ncfile


import matplotlib.colors

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s', datefmt='%d/%m/%Y %I:%M:%S')


def getColorMap( rgbFile = "medspiration.rgb" ):
    colors = []
    rgbFile = os.path.join(os.path.dirname(__file__), rgbFile )
    nbCol=0
    for line in open( rgbFile ):
        r,g,b = [int(c) for c in line.split()]
        colors.append( [r/255.,g/255.,b/255.] )
        nbCol += 1
    return( matplotlib.colors.ListedColormap(colors, name="custom", N=nbCol) )


f = mapper.ncfile.NCFile( URL = "/home/obistroz/data/public/swath/altimeters/waves/data/2000/01/wm_20000101.nc" )
feature = datamodel.pointcollection.PointCollection()
feature.load( f )
val = feature.get_lon()
print val
val = feature.get_values('wind_speed')[:]

feature.display_map(data=val, palette=getColorMap('wind.pal'), pretty=True, range=[0,20,0.25], output='tata.png')
