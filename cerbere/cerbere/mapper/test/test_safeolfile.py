'''
Created on 15 janv. 2015

@author: jfpiolle
'''
import sys

from cerbere.mapper.safeolfile import SAFEOLFile
from cerbere.datamodel.swath import Swath
from cerplot.mapping import CerMap


fname = sys.argv[1]
print fname
if 'S3A_OL_1_ERR' in fname or 'S3A_OL_1_EFR' in fname:
    testfield = "Oa01_radiance"
elif 'S3A_OL_2_LRR' in fname or 'S3A_OL_2_LFR' in fname:
    testfield = "IWV"
elif 'S3A_OL_2_WRR' in fname or 'S3A_OL_2_WFR' in fname:
    testfield = "Oa01_reflectance"

ncf = SAFEOLFile(url=fname)

print 'OPEN'
ncf.open()

print '\n\nFIELDS :'
fields = ncf.get_fieldnames()
for f in fields:
    print f

print '\n\nDIMENSIONS :'
dims = ncf.get_dimensions()
for f in dims:
    print f

print '\n\nATTRIBUTES :'
attrs = ncf.read_global_attributes()
for attr in attrs:
    print attr

print '\n\nGEOLOCATION DIMENSIONS'
print "Lat: ", ncf.get_dimensions('lat')
print "Lon: ", ncf.get_dimensions('lon')
print "Time: ", ncf.get_dimensions('time')

print '\n\nREAD GEOLOCATION'
lats = ncf.read_values('lat')
lons = ncf.read_values('lon')
times = ncf.read_values('time')

print times.shape, times[0, :]
print times.shape, times[1, :]


print '\n\nREAD START/END TIMES'
print ncf.get_start_time()
print ncf.get_end_time()

print "\n\nREAD SST FIELD"
print ncf.read_field(testfield)
print ncf.get_dimensions(testfield)
field = ncf.read_field(testfield)
print field.get_dimnames()
print ncf.get_full_dimensions(fieldname=testfield)
data = ncf.read_values(testfield)
print data.min(), data.max(), data.count()

ncf.close()

# LOAD THROUGH MODEL
print "\n\nLOAD SWATH MODEL"
swath = Swath()
ncf2 = ncf.__class__(url=fname)
swath.load(ncf2)

print "\n\nGET ALL FIELDS"
print swath.get_fieldnames()
for fieldname in swath.get_fieldnames():
    print swath.get_field(fieldname)


# DISPLAY
m = CerMap(swath, fieldname=testfield)
m.show()
