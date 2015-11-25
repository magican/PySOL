Reading Sentinel-1 OCN data
===========================

As usual, you have to instantiate a mapper of the correct class for this type of files, here : `SAFEOCNNCFile`::

  from cerbere.mapper.safeocnncfile import SAFEOCNNCFile
  
  fname = '/home/cercache/project/mpc-sentinel1/data/esa/sentinel-1a/L2/WV/S1A_WV_OCN__2S/2015/320/S1A_WV_OCN__2SSV_20151116T203517_20151116T205112_008634_00C443_EE48.SAFE/measurement/s1a-wv1-ocn-vv-20151116t204729-20151116t204732-008634-00c443-051.nc'
  # this creates the mapper object (equivalent to opening a file)
  fd = SAFEOCNNCFile(url=fname)

There is a subtility here : OCN files are actually a collation of 
3 different products (WIND, WAVE, DOPPLER) with different spatial 
properties (grid size and resolution). So when instantiating the mapper,
you need to specify with the `product` keyword which product you want
to read (by default, if the keyword is not specified like in above example : WIND). 
If you want to read all products, you need to open 3 different mappers 
(these products should really have been stored in different files!)::

  from cerbere.mapper.safeocnncfile import SAFEOCNNCFile, WIND, WAVE, DOPPLER
  
  # opening the WAVE product
  fd = SAFEOCNNCFile(url=fname, product=WAVE)


Discovery of the file content
-----------------------------

Similarly to the NetCDF API, Cerbere provides methods to inspect the content of a file as if it was self described (which is not the case of all formats).

Getting the list of fields::

  fieldnames = fd.get_fieldnames()
  print fieldnames


.. note::
  Note that only the fields from the opened product (ex: WIND) are returned. Spatial and temporal fields are filtered out.

Getting a specific field and printing its properties::

  field = fd.read_field('owiWindSpeed')
  print field

Getting the names of the file dimensions::

  dims = fd.get_dimensions()
  print dims

Getting the dimensions of a given field::

  dims = fd.get_dimensions(fieldname='owiWindSpeed')
  
  # or, using the field object :
  field = fd.read_field('owiWindSpeed')
  dims = field.get_dimnames()

You get the dimension's names and sizes with::

  # some field's dimensions
  dims = fd.get_full_dimensions(fieldname='owiWindSpeed')
  
  # or, with the field object :
  field = fd.read_field('owiWindSpeed')
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
--------------

Geolocation information
+++++++++++++++++++++++
if you don't use a data model, the geolocation information have to be queried like any field::

  lat = fd.read_values('lat')
  lon = fd.read_values('lon')

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
  `lat`, `lon` and `time` are standardized geolocation field names and will work with any mapper (whatever internal naming was used in the native file format)

Fields
++++++

Getting the values of any field::

  data = fd.read_values('owiWindSpeed')

Getting a subset using slices on spatial dimensions::

  data = fd.read_values('owiWindSpeed', slices={'row':slice(10,20), 'cell':slice(30, 40)})

