================
Writing a mapper
================

This section describes how to write a mapper that will allow access to a new
format in read mode (see at the end of the section complementary information
in case you also want to use this mapper to save data). 

Creating a new mapper module
============================
Writing a mapper consists in writing a set of function that helps cerbere to
understand and access a file content. It must implements the interface defined
by :class:AbstractMapper class, and therefore inherit this class.

Create a new module with a basic structure as follow:

.. code-block: python

	"""
	.. module::cerbere.mapper.<your mapper module name>
	
	Mapper classs for <the format and/or product type handled by this mapper>
	
	:license: Released under GPL v3 license, see :ref:`license`.
	
	.. sectionauthor:: <your name>
	.. codeauthor:: <your name>
	"""
	
	# import parent class
	from cerbere.mapper import abstractmapper
	
	
	class <Your mapper class name>(abstractmapper.AbstractMapper):
	    """Mapper class to read <the format and/or product type handled by this mapper> files"""
	
	    def __init__(self, url=None, mode=abstractmapper.READ_ONLY, **kwargs):
	        """Initialize a <the format and/or product type handled by this mapper> file mapper"""
	        super(<Your mapper class name>, self).__init__(url=url, mode=mode, **kwargs)
	        return



Functions to override:
  __init__
  open()
  close
  get_matching_dimname
  get_standard_dimname
  read_field
  read_values
  get_bbox
  get_dimensions
  get_dimsize
  get_geolocation_field
  get_spatial_resolution_in_deg
  get_start_time
  get_end_time
  read_field_attributes
  read_fillvalue
  read_global_attribute
  read_global_attributes



Saving data with a mapper
=========================

If the file can be also opened in write mode:
  create_dim
  create_field
  write_field
  write_global_attributes
