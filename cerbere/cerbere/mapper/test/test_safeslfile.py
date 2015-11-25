'''
Created on 15 janv. 2015

@author: jfpiolle
'''
import sys
import os
import argparse
import importlib

from cerbere.datamodel.swath import Swath
from cerbere.mapper.ncfile import NCFile, WRITE_NEW
from cerplot.mapping import CerMap


SLSTR_MAPPERS = ['SAFESLIRFile', 'SAFESL500AFile',
                 'SAFESL500BFile', 'SAFESL500TDIFile']
TEST_FIELD = {
    'SAFESLIRFile': 'N3_sea_surface_temperature',
    'SAFESL500AFile': 'S2_radiance_an',
    'SAFESL500BFile': 'S4_radiance_bo',
    'SAFESL500TDIFile': 'S4_radiance_cn'
    }
OBLIQUE_LAT_FIELD = {
    'SAFESLIRFile': 'latitude_io',
    'SAFESL500AFile': 'latitude_ao',
    'SAFESL500BFile': 'latitude_bo',
    'SAFESL500TDIFile': 'latitude_co'
    }

parser = argparse.ArgumentParser(
    description="Test reading a SLSTR product")

# Set need argument
parser.add_argument("input",
                    help="full path to the input file (file to be read)")
parser.add_argument("mapper",
                    help="name of the mapper class to use")
args = parser.parse_args()

fname = args.input
if not args.mapper or args.mapper not in SLSTR_MAPPERS:
    raise Exception("No valid mapper")
else:
    mapper = args.mapper
    source_mapper = getattr(
        importlib.import_module(
            'cerbere.mapper.safeslfile'
            ),
        mapper,
        )

ncf = source_mapper(url=fname)
testfield = TEST_FIELD[mapper]


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
values_o = ncf2.read_values(OBLIQUE_LAT_FIELD[mapper])
print values_o.shape
print "Values at row 1000 :"
print values_i[1000, 550:600]
print values_o[1000, 550:600]
diff = (values_o - values_i)
print diff.min(), diff.max()

print "\nREAD SUBSET 1"
values_i = swath.get_lat(slices={'row': slice(1000, 1001),
                                 'cell': slice(550, 600)},
                         cache=False)
values_o = swath.get_values(OBLIQUE_LAT_FIELD[mapper],
                            slices={'row': slice(1000, 1001),
                                    'cell': slice(550, 600)},
                            cache=False)

print values_i, values_i.shape
print values_o, values_o.shape

print "\nREAD SUBSET 2"
values_i = swath.get_lat(slices={'row': slice(1000, 1001),
                                 'cell': slice(1200, 1250)},
                         cache=False)
values_o = swath.get_values(OBLIQUE_LAT_FIELD[mapper],
                            slices={'row': slice(1000, 1001),
                                    'cell': slice(1200, 1250)},
                            cache=False)

print "Nadir :", values_i, values_i.shape
print "Oblique : ", values_o, values_o.shape

print "\nREAD SUBSET 3"
values_i = swath.get_lat(slices={'row': slice(1000, 1001),
                                 'cell': slice(0, 50)},
                         cache=False)
values_o = swath.get_values(OBLIQUE_LAT_FIELD[mapper],
                            slices={'row': slice(1000, 1001),
                                    'cell': slice(0, 50)},
                            cache=False)

print "Nadir :", values_i, values_i.shape
print "Oblique : ", values_o, values_o.shape

print "\nREAD SUBSET 4"
values_i = swath.get_lat(slices={'row': slice(1000, 1001),
                                 'cell': slice(1000, 1300)},
                         cache=False)
values_o = swath.get_values(OBLIQUE_LAT_FIELD[mapper],
                            slices={'row': slice(1000, 1001),
                                    'cell': slice(1000, 1300)},
                            cache=False)

print "Nadir :", values_i, values_i.shape
print "Oblique : ", values_o, values_o.shape

# extract subset and save
subset = swath.extract_subset(slices={'row': slice(1000, 1001),
                                      'cell': slice(500, 550)})
# save subset
print("Save subset")
subsetfname = 'subset.nc'
if os.path.exists(subsetfname):
    os.remove(subsetfname)
oncf = NCFile(url=subsetfname, mode=WRITE_NEW, ncformat='NETCDF4')
subset.save(oncf)

# DISPLAY
#m = CerMap(swath, fieldname=testfield, area=[-8., 48., 24., 60.])
m = CerMap(swath, fieldname=testfield)
m.show()
