'''
Created on 18 juin 2012

Test reading/displaying swath file

@author: jfpiolle
'''
import logging

import datamodel.swath
import mapper.qscathdffile
import plot.palette

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s', datefmt='%d/%m/%Y %I:%M:%S')


ncf = mapper.qscathdffile.QSCATHDFFile( URL = "/home/cerdata/provider/podaac/satellite/l2b/quikscat/seawinds/25km/200727/QS_S2B41907.20071891450" )

swath = datamodel.swath.Swath()
swath.load( ncf )

swath.display_map('wind_speed_selection', range=[0,25,0.25], contouring="contour", gridStep=0.50, \
                palette=plot.palette.getColorMap('wind.pal'), output='swath_qscatl2b.png')
