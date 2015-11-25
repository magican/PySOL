'''
Created on 9 aoit 2013

@author: jfpiolle
'''
import logging

from cerbere.datamodel.swath import Swath
from cerbere.mapper.ascatifrncfile import AscatIFRNCFile


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-5s %(message)s', datefmt='%d/%m/%Y %I:%M:%S')


ncf = AscatIFRNCFile(url='/home3/bagan/private/swath/ascat/bufr/Cdf_orbit_250//2010/07/06/19266.nc')
swath=Swath()
swath.load(ncf)

print swath.get_bbox()

print swath.get_fieldnames()

print swath.get_times()