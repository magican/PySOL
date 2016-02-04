"""
checker
=======

Tool for checking a new mapper

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
import sys
import os
import logging
import importlib
import traceback

import netCDF4

from cerbere.mapper.ncfile import NCFile
from cerbere.mapper.abstractmapper import WRITE_NEW

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)-5s %(message)s',
                    datefmt='%d/%m/%Y %I:%M:%S')

model = sys.argv[1]
mapper = sys.argv[2]
infile = sys.argv[3]

logging.info("1. Checking mapper class %s", mapper)
logging.info("========================")

logging.info("Load mapper class")
if '.' in mapper:
    mapperpath = 'cerbere.mapper.' + mapper.split('.')[0].lower()
    mappername = mapper.split('.')[1]
else:
    mapperpath = 'cerbere.mapper.' + mapper.lower()
    mappername = mapper
classreader = getattr(
    importlib.import_module(mapperpath),
    mappername
    )
logging.info("OK")
logging.info("")

logging.info("Open file")
f = classreader(url=infile)
f.open()
logging.info("OK")
logging.info("")

logging.info("Read fields")
fields = f.get_fieldnames()
for fieldname in fields:
    logging.info('-------- %s --------', fieldname)
    try:
        field = f.read_field(fieldname)
        logging.info(field)
    except:
        logging.error("Could not read field %s", fieldname)
        logging.error(traceback.print_exc())
logging.info("")

logging.info("Read latitudes")
data = f.read_values('lat')
logging.info(data)
logging.info("min : %s", data.min())
logging.info("max : %s", data.max())
logging.info("")

logging.info("Read longitudes")
data = f.read_values('lon')
logging.info(data)
logging.info("min : %s", data.min())
logging.info("max : %s", data.max())
logging.info("")


logging.info("Read times")
data = f.read_values('time')
logging.info(data)
logging.info("min : %s", netCDF4.num2date(data.min(), f.read_field(
                f.get_geolocation_field('time')).units))
logging.info("max : %s", netCDF4.num2date(data.max(), f.read_field(
                f.get_geolocation_field('time')).units))
logging.info("")

logging.info("Read metadata")
logging.info('start time : %s', f.get_start_time())
logging.info('end time : %s', f.get_end_time())
attributes = f.read_global_attributes()
for attr in attributes:
    logging.info("%s : %s", attr, f.read_global_attribute(attr))


f.close()

logging.info("2. Testing datamodel")
logging.info("====================")
modelreader = getattr(
            importlib.import_module(
                    'cerbere.datamodel.' + model.lower()
                    ),
            model
            )

logging.info("Instantiate the datamodel")
f = classreader(url=infile)
modelobj = modelreader()
modelobj.load(f)


logging.info("3. Testing subsetting")
logging.info("=====================")

logging.info("Create subset")

geodims = modelobj.get_geolocation_dimsizes().values()

width = 5

if model == 'Swath':
    rows, cells = geodims
    r0, r1 = rows / 2 - width, rows / 2 + width
    c0, c1 = cells / 2 - width, cells / 2 + width
    print "Subset "
    print "row : ", r0, r1
    print "cell: ", c0, c1
    subset = modelobj.extract_subset(slices={'row': slice(r0, r1, 1),
                                             'cell': slice(c0, c1, 1)})
elif model == 'Image':
    rows, cells = geodims
    r0, r1 = rows / 2 - width, rows / 2 + width
    c0, c1 = cells / 2 - width, cells / 2 + width
    print "Subset "
    print "row : ", r0, r1
    print "cell: ", c0, c1
    subset = modelobj.extract_subset(slices={'row': slice(r0, r1, 1),
                                             'cell': slice(c0, c1, 1)})
elif model == 'Grid':
    nj, ni = geodims
    j0, j1 = nj / 2 - width, nj / 2 + width
    i0, i1 = ni / 2 - width, ni / 2 + width
    subset = modelobj.extract_subset(slices={'y': slice(j0, j1, 1),
                                             'x': slice(i0, i1, 1)})

# save subset
logging.info("Save subset")
subsetfname = 'hrdds21.nc'
if os.path.exists(subsetfname):
    os.remove(subsetfname)
oncf = NCFile(url=subsetfname, mode=WRITE_NEW, ncformat='NETCDF4')
subset.save(oncf)
oncf.close()

# read subset
logging.info("Read subset")
f = NCFile(url=subsetfname)
modelobj = modelreader()
modelobj.load(f)
