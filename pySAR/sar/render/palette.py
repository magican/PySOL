'''
Created on 16 august 2012

@author: jfpiolle
'''

import os

import matplotlib


def getColorMap( rgbFile = "medspiration.rgb" ):
    '''
    Load a RGB palette provided in ascii file
    '''
    colors = []
    rgbFile = os.path.join(os.path.dirname(__file__), 'palette', rgbFile )
    nbCol=0
    for line in open( rgbFile ):
        r,g,b = [int(c) for c in line.split()]
        colors.append( [r/255.,g/255.,b/255.] )
        nbCol += 1
    return( matplotlib.colors.ListedColormap(colors, name="custom", N=nbCol) )
