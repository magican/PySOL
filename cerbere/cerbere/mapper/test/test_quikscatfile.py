'''
Created on 9 aoit 2013

@author: jfpiolle
'''
import logging

from cerbere.datamodel.swath import Swath
from cerbere.mapper.qscathdffile import QSCATHDFFile
from cerbere.plot.mapping import Mapping

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s', datefmt='%d/%m/%Y %I:%M:%S')

ncf = QSCATHDFFile(url='/home/cerdata/provider/podaac/satellite/l2b/quikscat/seawinds/25km/2009/271/QS_S2B53526.20092722205')

swath = Swath()
swath.load(ncf)

#print swath.get_bbox()




lons = swath.get_lon()
for i in range(lons.shape[1]):
    print lons[:,i].min(), lons[:,i].max(), lons[:,i]

lats = swath.get_lat()
for i in range(lats.shape[1]):
    print lats[:,i].min(), lats[:,i].max(), lats[:,i]

#qmap = Mapping(swath, fieldname='lon', maskland=False, legendposition=None)
#qmap.save('qscat2.png')

from matplotlib import pyplot

toto = pyplot.imshow(lons)
pyplot.colorbar()
pyplot.savefig('toto.png')