==================
General principles
==================

Concept
=======

The main concept of **cerbere** is that each data acquisition corresponds to a well-known observation or sampling pattern, referred to as *data model* or *feature*. The package provides a set of classes implementing each data model, independent from the way the data are stored on disk or into a database (therefore storage format and structure). Storage (or *mapper*) classes are used to map the content of a file to the corresponding feature (e.g. the content of a grid file to a *grid* feature).

We assume indeed there is no practical reason why data corresponding to the same sampling pattern (or structure) would be represented differently. Having a set of predefined template for each feature type allows to write for once all usual generic operations such as display, extraction of values, remapping or resampling, saving to the same format conventions,...

Currently managed features include :

 * Grid
 * Swath
 * Image
 * GridTimeSeries
 * Trajectory
 * PointTimeSeries

.. Important::
  **cerbere** clearly separates the content from the format:
    * **content typing** : data structure or memory representation of the data,
      which should be unique for each feature, described in more details in 
      :doc:`datamodel`. This corresponds to the :mod:`~cerbere.datamodel`
      package.
    * **storage format** : the way the above structure is stored/mapped on file,
      which can be very different and exotic, described in more details in
      :doc:`storage`. This corresponds to the :mod:`~cerbere.mapper` package.

.. Note::
  the :mod:`~cerbere.mapper` package can be used independently, as a unified
  API to read the content from any data file.
  
  The :mod:`~cerbere.datamodel` requires the usage of the 
  :mod:`~cerbere.mapper` package to read or write the data into/from a
  :mod:`~cerbere.datamodel` object.

Using the `mapper` package
==========================

The following section describes the basic operations you can perform with 
**cerbere** to handle Earth Observation data files with 
:mod:`~cerbere.mapper` package. This package can be seen as a unified API to
access any data file content, whatever the format. There is one 
:mod:`~cerbere.mapper` class per data format.

Reading data from a file
------------------------

A :mod:`~cerbere.mapper` class must be available for the considered file.
Generic mapper classes are available, for instance :mod:`~cerbere.mapper.ncfile.NCFile`
for CF compliant NetCDF files.

.. Note::
   If no mapper class exists for a particular format, a new corresponding
   :mod:`~cerbere.mapper` class must be written.
   
   The complete list of existing mappers, and their compatibility with known
   datasets is listed in :doc:`compatibility`.

To read data from a file, first instantiate a mapper object of the
corresponding class, specifying the path to this file in the `url` argument:

.. doctest::

   >>> import mapper.ncfile
   >>> ncf = mapper.ncfile.NCFile(url="./test/GW_L2P_ALT_ENVI_GDR_20101210_120905_20101210_125912_097_196.nc")

.. Warning::
   This does not open the file. The file must be explicitly opened with :func:`open` function:

.. doctest::

   >>> f.open()

A mapper provides a set of methods to inspect the content of a file. They allow
to retrieve information from a file in the same way whatever its format.

Get the list of fields in a file (all but the geolocation fields) with *get_fieldnames*:

.. doctest::

   >>> f.get_fieldnames()

Get the dimensions (like in netCDF) of a file (Note that the geolocation dimension 
names returned are standardized)::

    f.get_dimensions()

Get the dimensions (like in netCDF) of a particular field::

    f.get_dimensions('sea_surface_temperature')
    
Get the size of a dimension (standard names can be used for geolocation dimensions)::

    f.get_dimsize('row')    # standard dimension name
    f.get_dimsize('ni')     # equivalent native name

Get a field and display it::

    field = f.read_field('sea_surface_temperature')
    print field

.. note::
    *Fields* are similar to variables in netcdf. A field consists of :
        * an attached *variable* describing the geophysical quantity provided by the field (together with a few descriptive attributes such standard name, etc...)
        * *attributes* further documenting the provided observation values (units,...) similar to the variable attributes in netCDF
        * an *array of values* (observations)
        * an optional array of *quality flags* (one for each observation value)
        * an optional array of *quality history*  (one for each observation value) documenting the reason why a value was flagged

Get the list of field attributes, as a dictionary::

    attr = read_field_attributes('sea_surface_temperature')
    print attr

Get the list of global attributes::

    attr = f.read_global_attributes()
    print attr


Then instantiate a data model class object and load the above mapper object. The data model object is now initialized with the content of the above file :

.. doctest::

 >>> import datamodel.trajectory
 >>> traj = datamodel.trajectory.Trajectory()
 >>> traj.load(ncf)



Saving data
-----------

Extracting sub-features
-----------------------

