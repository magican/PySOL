'''
Created on 15 janv. 2015

@author: jfpiolle
'''
import sys

test = int(sys.argv[1])

if test == 0:
    from cerbere.mapper.safeslfile import SAFESLIRFile
    fname = "/export/home/sumba/project/felyx/data/S3A_SL_2_WCT____20080802T025305_20080802T044153_20150513T151745_6528_091_303______MAR_D_NT_001.SEN3"
    #fname = "/home1/ananda/data/v1/S3A_SL_2_WCT____20080802T025305_20080802T044153_20150513T151745_6528_091_303______MAR_D_NT_001.SEN3"
    ncf = SAFESLIRFile(url=fname)
    testfield = "SST"
elif test == 1:
    from cerbere.mapper.safeslfile import SAFESL500AFile
    fname = "/home1/ananda/data/S3A_SL_1_RBT____20130621T100932_20130621T101146_20141201T130554_0133_001_002______LN1_D_NR____.SEN3"
    ncf = SAFESL500AFile(url=fname)
    testfield = "S2_radiance_an"
elif test == 2:
    from cerbere.mapper.safeslfile import SAFESL500BFile
    fname = "/home1/ananda/data/S3A_SL_1_RBT____20130621T100932_20130621T101146_20141201T130554_0133_001_002______LN1_D_NR____.SEN3"
    ncf = SAFESL500BFile(url=fname)
    testfield = "S4_radiance_bo"
elif test == 3:
    from cerbere.mapper.safeslfile import SAFESL500TDIFile
    fname = "/home1/ananda/data/S3A_SL_1_RBT____20130621T100932_20130621T101146_20141201T130554_0133_001_002______LN1_D_NR____.SEN3"
    ncf = SAFESL500TDIFile(url=fname)
    testfield = "S4_radiance_cn"


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

print '\n\nREAD START/END TIMES'
print ncf.get_start_time()
print ncf.get_end_time()

print "\n\nREAD SST FIELD"
print ncf.read_field(testfield)
print ncf.get_dimensions(testfield)
field = ncf.read_field(testfield)
print field.get_dimnames()
print ncf.get_full_dimensions(fieldname=testfield)

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

print "\nREAD VALUES"
values_i = ncf2.read_values('lat')
print values_i.shape
values_o = ncf2.read_values('latitude_io')
print values_o.shape
print "Values at row 1000 :"
print values_i[1000, 500:550]
print values_o[1000, 500:550]
diff = (values_o - values_i)
print diff.min(), diff.max()

print "\nREAD SUBSET 1"
values_i = swath.get_lat(slices={'row':slice(1000,1001), 'cell': slice(500, 550)}, cache=False)
values_o = swath.get_values('latitude_io', slices={'row':slice(1000,1001), 'cell': slice(500, 550)}, cache=False)

print values_i, values_i.shape
print values_o, values_o.shape

print "\nREAD SUBSET 2"
values_i = swath.get_lat(slices={'row':slice(1000,1001), 'cell': slice(1200, 1250)}, cache=False)
values_o = swath.get_values('latitude_io', slices={'row':slice(1000,1001), 'cell': slice(1200, 1250)}, cache=False)

print "Nadir :", values_i, values_i.shape
print "Oblique : ", values_o, values_o.shape

print "\nREAD SUBSET 3"
values_i = swath.get_lat(slices={'row':slice(1000,1001), 'cell': slice(0, 50)}, cache=False)
values_o = swath.get_values('latitude_io', slices={'row':slice(1000,1001), 'cell': slice(0, 50)}, cache=False)

print "Nadir :", values_i, values_i.shape
print "Oblique : ", values_o, values_o.shape

print "\nREAD SUBSET 4"
values_i = swath.get_lat(slices={'row':slice(1000,1001), 'cell': slice(1000, 1300)}, cache=False)
values_o = swath.get_values('latitude_io', slices={'row':slice(1000,1001), 'cell': slice(1000, 1300)}, cache=False)

print "Nadir :", values_i, values_i.shape
print "Oblique : ", values_o, values_o.shape


# DISPLAY
from cerbereutils.plot.mapping import CerMap

#m = CerMap(swath, testfield, area=[-8., 48., 24., 60.])
m = CerMap(swath, testfield, contouring='scatter')
m.save('toto.png')

