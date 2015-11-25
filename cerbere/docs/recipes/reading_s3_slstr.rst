=============================
Reading Sentinel-3 SLSTR data
=============================

Important foreword
==================

As usual, you have to instantiate a mapper of the correct class for this type
of files. However, there is not a single mapper for all SLSTR product types
and, worst, there is not even a single mapper for the same SAFE product.

This is due to the fact that the same SAFE product mixes different views and
grids, each with different dimensions.

.. important::

  From **cerbere** perspective, each grid is considered as a different product,
  with its own mapper class. However it merges oblique and nadir views as if
  they were contained in the same file, so you don't need to instanciate a
  mapper for each view.

Note that all mapper classes inherit from the same parent class, `SAFESLFile`.

Refer to the following table to decide which mapper you need to use:

+------------------------------------------------+------------------------+---------------+
| Dataset                                        | Mapper class           | Datamodel     |
+================================================+========================+===============+
| S3A_SL_2_WCT Nadir (in)                        | SAFESLIRFile           | Swath         |
| S3A_SL_1_RBT 1km Nadir (in)                    |                        |               |
| S3A_SL_2_WCT Oblique (io)                      |                        |               |
| S3A_SL_1_RBT 1km Oblique (io)                  |                        |               |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT 500m & SWIR A Stripe Nadir (an)   | SAFESL500AFile         | Swath         |
| S3A_SL_1_RBT 500m & SWIR A Stripe Oblique (ao) |                        |               |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT 500m & SWIR B Stripe Nadir (bn)   | SAFESL500BFile         | Swath         |
| S3A_SL_1_RBT 500m & SWIR B Stripe Oblique (bo) |                        |               |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT TDI Nadir (cn)                    | SAFESL500TDIFile       | Swath         |
| S3A_SL_1_RBT TDI Oblique (co)                  |                        |               |
+------------------------------------------------+------------------------+---------------+

Let's for instance work with the SAFE L2 WCT product::

  S3A_SL_2_WCT____20130621T101013_20130621T101053_20141201T092032_0039_009_022______MAR_O_NR_001.SEN3

We will open first the Nadir view:

.. code-block:: python

  from cerbere.mapper.safeslfile import SAFESLIRFile
  fname = 'S3A_SL_2_WCT____20130621T101013_20130621T101053_20141201T092032_0039_009_022______MAR_O_NR_001.SEN3'
  # this creates the mapper object (equivalent to opening a file)
  fd = SAFESLIRFile(url=fname)


Discovery of the file content
=============================

Similarly to the NetCDF API, **cerbere** provides methods to inspect the
content of a file as if it was self described (which is not the case of all
formats).

Getting the list of fields::

  fieldnames = fd.get_fieldnames()
  print fieldnames


.. note::
  Note that only the geophysical fields are returned. Geolocation (spatial and
  temporal) fields are filtered out.

Getting a specific field and printing its properties::

  field = fd.read_field('SST')
  print field

Getting the names of the file dimensions::

  dims = fd.get_dimensions()
  print dims

Getting the dimensions of a given field::

  dims = fd.get_dimensions(fieldname='SST')
  
  # or, using the field object :
  field = fd.read_field('SST')
  dims = field.get_dimnames()

You get the dimension's names and sizes with::

  # some field's dimensions
  dims = fd.get_full_dimensions(fieldname='SST')
  
  # or, with the field object :
  field = fd.read_field('SST')
  size = field.get_dimsize('row')


Getting the global attributes::

  attributes = fd.read_global_attributes()
  for attr in attributes:
      print attr, fd.read_global_attribute(attr)

Getting some specific (normalised) attributes::

  # start time
  print fd.get_start_time()
  
  # end time
  print fd.get_end_time()


Reading values
==============

Geolocation information
-----------------------
if you don't use a data model, the geolocation information have to be queried like any field::

  lat = fd.read_values('lat')
  lon = fd.read_values('lon')
  z = fd.read_values('z')

  times = fd.read_values('time')

  # or, using fields
  field = fd.read_field('time')
  times = field.get_values()

  # convert to datetime object
  from netCDF4 import num2date
  field = fd.read_field('time')
  times = field.get_values()
  print num2date(times[:], field.units)

.. note::
  `lat`, `lon`, `z` and `time` are standardized geolocation field names and will
  work with any mapper (whatever internal naming was used in the native file format)

Fields
------
Getting the values of any field::

  data = fd.read_values('SST')

Getting a subset using slices on spatial dimensions::

  data = fd.read_values('SST', slices={'row':slice(10,20), 'cell':slice(30, 40)})

Using a data model
==================

The content of the file can be mapped into a data model which is convenient for
operations using these datamodel.

In the case of the WCT L2 file in nadir view, we will use the ``Swath`` model
as listed in the above table.::

  from cerbere.datamodel.swath import Swath
  swath = Swath()

Load the content from a file into the model, thanks to the mapper already
seen, using the `load` function::

  from cerbere.mapper.safeslfile import SAFESLFile
  fname = 'S3A_SL_2_WCT____20130621T101013_20130621T101053_20141201T092032_0039_009_022______MAR_O_NR_001.SEN3'
  # this creates the mapper object (equivalent to opening a file)
  fd = SAFESLFile(url=fname)
  
  swath.load(fd)

Read the lat, lon, z and times::

  lats = swath.get_lat()
  lons = swath.get_lon()
  z = swath.get_z()
  times = swath.get_times()

In above example, times are returned as numbers than can be converted to
datetime objects geting first the units::

  units = swath.get_time_units()
  
  import netCDF4
  times2 = netCDF4.num2date(times, units)

Note the following function can be used::

  times2 = swath.get_datetimes()

.. warning::

  The conversion from time to datetime objects can be very long and must be
  avoided at all cost.

Slices can be used to get a subset of data::

  data = fd.get_values('SST', slices={'row':slice(10,20), 'cell':slice(30, 40)})

A complete subset of the file can be obtained as follow::

  subset = swath.extract_subset(slices={'row': slice(10, 20),
                                        'cell': slice(30, 40)})

The result is a new `Swath` object, but smaller. It contains the same list of
fields as the original parent object.

It can be saved into netCDF using a `NCFile` mapper::

  newncf = NCFile(url=subsetfname, mode=WRITE_NEW, ncformat='NETCDF4')
  subset.save(newncf)

The data can also be displayed using the `cerbereutils` package::

  from cerbereutils.plot.mapping import CerMap
  
  m = CerMap(swath, 'SST', area=[-8., 48., 24., 60.])
  m.save('swath.png')

