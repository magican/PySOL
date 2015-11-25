'''
Created on 30 juil. 2013

@author: jfpiolle
'''

import mapper.gribfile

f = '/home/cerdata/provider/ncep/model/cfsr/monthly/1999/03/pgbhnl.gdas.199903.grb2'

gribm = mapper.gribfile.GribFile(url=f)

gribm.open()

attrs = gribm.read_field_attributes('msl')
for k in attrs:
    print k, ' : ', attrs[k]

print "LAT"
field = gribm.read_field('lat')
print field
data = field.get_values()
print data

print "LON"
field = gribm.read_field('lon')
print field
data = field.get_values()
print data

print "TIME"
field = gribm.read_field('time')
print field
data = field.get_values()
print data

print "MSL"
field = gribm.read_field('msl')
print field
data = field.get_values()
print data
