=======
Storage
=======

Data features are stored into files in a wide range of ways, formats or structures. 
We standardize the access to any data format through specific interface classes, 
called *mappers*.

A **mapper** class behaves similarly to the netCDF API for instance, providing introspection
functions for a file content, with a little bit more intelligence : it also provides
some functions describing the structure of the data layed out in the file for an easy
mapping to a **datamodel** class.

This means the relationship between the file variables, metadata, or dimensions, and some feature's 
properties must be fully understood and explicited, following the abstract interface defined 
in the package's parent class *AbstractMapper*. 

**cerbere** comes with a set of predefined **mapper** classes but any new format 
(or formatting specificities) must be handled through a new **mapper** class to be written 
by the user.

This section describes how to understand and write a new **mapper** class.

Heritage
========

Any new mapper class must inherit from *AbstractMapper* class. This abstract class 
provides the list of methods of the mapper interface that need to be implemented.

**mapper** classes already exist for some standard formats such netCDF, HDF4. 
These formats are very generic but do not standardize the expression of feature information
(geolocation information,...) unless some convention is used (such as CF for netCDF). 
It may therefore be necessary to override these classes in a new more specialized inherited 
class.


Standardization methods
=======================

A datamodel class (for instance *Grid*) is expecting some information on the data 
structure (for instance the name and size of each dimension, the time and spatial coordinates). 
These information must be returned by a *mapper* class in a standard way.

Some internal methods (never called directly by a user and only used by *datamodel* classes) 
help realising the matching between the expected information by a *datamodel* class to instantiate 
and the information specified in the proprietary file format (following the data provider convention which is not 
standard).

*get_geolocation_field()* returns the internal name (in the proprietary file format) of any spatial and temporal 
coordinate field (*time*, *lat* or *lon*)::

    f.get_geolocation_field('time')
    f.get_geolocation_field('lat')
    f.get_geolocation_field('lon')

.. note::
    in self described format (netCDF, HDF) this will generally consist in providing the equivalent field 
    in the proprietary file format. If the field is not existing (for instance time information is sometimes 
    splitted into several fields for day, milliseconds, etc) or if the file format is not self-described (pure binary 
    file), virtual names are returned. They should be equivalent to the standard names.
 
*get_matching_dimname()* returns the internal name (in the proprietary file format) of any standard geolocation dimension::

    f.get_matching_dimname('row')
    f.get_matching_dimname('cell')

.. note::
    Refer to datamodel section for the list of standard geolocation dimension for each feature.

