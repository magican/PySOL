'''
Created on 9 aoit 2013

@author: jfpiolle
'''
import logging

from cerbere.datamodel.grid import Grid
from cerbere.mapper.ghrsstncfile import GHRSSTNCFile


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s', datefmt='%d/%m/%Y %I:%M:%S')


ncf = GHRSSTNCFile(url='/home/cerdata/provider/ghrsst/satellite/l4/glob/odyssea-nrt/data/2015/075/20150316-IFR-L4_GHRSST-SSTfnd-ODYSSEA-GLOB_010-v2.0-fv1.0.nc')
print ncf.get_fieldnames()

feature = Grid()
feature.load(ncf)

print feature.get_bbox()
print feature.get_fieldnames()
print feature.get_times()
