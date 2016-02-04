'''
Created on 15 janv. 2015

@author: jfpiolle
'''
import sys

from cerbere.mapper.safeolfile import SAFEOLFile

test = int(sys.argv[1])


if test == 0:
    fname = "/home1/ananda/data/v1/S3A_OL_1_ERR____20080108T065652_20080108T074025_20150520T171028_2613_065_006______LN2_D_NT____.SEN3"
    ncf = SAFEOLFile(url=fname)
    testfield = "Oa01_radiance"
elif test == 1:
    fname = "/home1/ananda/data/v1/S3A_OL_2_LRR____20080108T065652_20080108T074025_20150521T125214_2613_084_049______LN1_D_NT_001.SEN3"
    ncf = SAFEOLFile(url=fname)
    testfield = "IWV"
elif test == 2:
    fname = "/home1/ananda/data/v1/S3A_OL_2_WRR____20080108T065652_20080108T074025_20150521T130759_2613_084_049______MAR_D_NT_001.SEN3"
    ncf = SAFEOLFile(url=fname)
    testfield = "Oa01_reflectance"


from cerbere.datamodel.swath import Swath


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
from cerbereutils.plot.mapping import CerMap

m = CerMap(swath, testfield, contouring='scatter')
m.save('toto.png')

