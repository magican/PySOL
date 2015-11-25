# -*- coding: utf-8 -*-
"""
cerbere.datamodel.field
=======================

Classes for the handling the data fields

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
import logging
import copy
import numpy
from functools import partial

from .variable import Variable
import cerbere.mapper.slices

__all__ = ['QCLevel', 'QCDetail', 'Field']


class QCLevel(object):
    """Special layer for the storage of the quality level associated
     with a field.
    """
    def __init__(self, values=None, dimensions=None, levels=None,
                 meanings=None):
        self.values = values
        self.dimensions = dimensions
        self.levels = levels
        self.meanings = meanings
        return


class QCDetail(object):
    """Special layer for the storage of the quality control details associated
     with a field. This is intended to describe which quality steps have been
     applied and raised a warning/error (bit value = 1)
    """
    def __init__(self, values=None, dimensions=None, mask=None, meanings=None):
        self.values = values
        self.dimensions = dimensions
        self.mask = mask
        self.meanings = meanings
        return


class Field(object):
    """A Field contains data and associated metadata.

    This is equivalent to the notion of variable in netCDF (not to be confused
    with the definition of :class:`Variable` used in Cerbere.

        Args:
            variable (:class:`Variable`) : the physical variable measured and
                stored in the field

            dimensions (:class:`collections.OrderedDict`) : an ordered
                dictionary providing the name and size of each dimension

            values (:class:`numpy.ma.MaskedArray`) : the measured data.
                Optional (the field can be created with default values)

            datatype (:class:`numpy.dtype`) : the type of the data to be
                stored. If `values` is provided, the datatype does not need
                to be provided as the field datatype will be the datatype of
                the provide values

            fields (:class:`Field`) : the subfields composing the main field.
                This is intended to group for instance components of the same
                variable (such as a vector's northward and eastward components
                for wind, currents,...). This allows to relate these components
                to the same physical variable (e.g wind) and to a single 
                qc_levels and qc_details information.

            fillvalue : the default value to associate with missing data in the
                field's values. The fillvalue must be of the same type as
                `datatype` and `values`

            attributes (dictionary) : a dictionary of the metadata associated
                with the field's values

            units (string) : the units in which the data are given (if
                applicable)

            allocate (bool) : if True, pre-allocate the data arrays to store the values 
                (in case no values are yet provided)
    """
    def __init__(self, variable, dimensions, values=None, datatype=None,
                 fields=None, fillvalue=None, attributes={},
                 units=None, valid_min=None, valid_max=None,
                 qc_levels=None, qc_details=None, allocate=False):
        """
        """
        object.__init__(self)
        self.variable = variable
        self.dimensions = dimensions
        self.attributes = {}
        if fields:
            # Add components in case of a composite field :
            # each of these components must be itself a field
            self.components = fields
        else:
            # simple scalar field (no components)
            self.handler = None
            self.attributes = {}
            self.components = None
            self.valid_min = valid_min
            self.valid_max = valid_max
            self.units = units
            self._values = values
            if values is not None and datatype is None:
                self.datatype = values.dtype
            else:
                self.datatype = datatype
            self.fillvalue = fillvalue
            self.attributes = attributes
            self.qc_levels = qc_levels
            self.qc_details = qc_details
            if values is None and allocate:
                # Initialisation of an empty record
                if self.datatype is None:
                    raise Exception('Missing datatype')
                self._values = numpy.ma.masked_all(
                                tuple(self.dimensions.values()),
                                self.datatype
                                )
                if fillvalue:
                    self._values.set_fill_value(fillvalue)
        return

    def __str__(self):
        result = 'field : %s\n' % self.name
        result = result + '\tdimensions :\n'
        for dim in self.get_dimnames():
            result = result + '   \t  # %s : %s\n'\
                % (dim, self.get_dimsize(dim))
        result = result + '\tattributes :\n'
        for att in self.attributes:
            result = result + '   \t  # %s : %s\n'\
                % (att, self.attributes[att])
        result = result + '\tconventions :\n'
        result = result + '   \t  # Units : %s\n' % self.units
        result = result + '   \t  # Fillvalue : %s\n' % self.fillvalue
        return result

    @classmethod
    def format_slices(cls, slices, fielddims):
        """format the slices to an explicit list matching the list and order of
        a field dimensions

        the slices as provided as a dictionary where keys are the dimensions
        over which subsetting needs to be performed. Implicit dimensions where
        no subsetting is requested don't need to be specified. The completion
        to a comprehensive ordered list of slices for all field dimensions is
        performed by this function.

        Args:
            slices (dict): keys are the names and values are the slice object
                of the dimensions in which to subset

            fielddims (list): full list of the field dimension names

        Returns:
            list: a comprehensive list of slice objects for all field
                dimensions.
        """
        return cerbere.mapper.slices.format_slices(slices, fielddims)

    @classmethod
    def format_indices(cls, indices, fielddims):
        """format the indices to an explicit list of `slice` matching the list
        and order of a field dimensions

        the indices are provided as a dictionary where keys are the dimensions
        over which subsetting needs to be performed. Implicit dimensions where
        no subsetting is requested don't need to be specified. The completion
        to a comprehensive ordered list of slices for all field dimensions is
        performed by this function.

        Args:
            indices (dict): keys are the names and values are the indices
                of the dimensions in which to subset

            fielddims (list): full list of the field dimension names

        Returns:
            list: a comprehensive list of slice objects for all field
                dimensions.
        """
        return cerbere.mapper.slices.format_indices(indices, fielddims)

    def get_metadata(self):
        """returns all the metadata of the field"""
        return (self.units, self.fillvalue,
                self.valid_min, self.valid_max,
                copy.copy(self.attributes))

    def set_metadata(self, metadata):
        """set the metadata of the field"""
        units, fillvalue, valid_min, valid_max, attributes = metadata
        self.valid_min = valid_min
        self.valid_max = valid_max
        self.units = units
        self.fillvalue = fillvalue
        self.attributes = copy.copy(attributes)
        return
    metadata = property(get_metadata, set_metadata)

    def get_name(self):
        return self.variable.shortname

    def set_name(self, name):
        self.variable.shortname = name
    name = property(get_name, set_name)

    def is_composite(self):
        """
        return True if the field is a composite field, i.e. it is a composition
        of sub fields (vector components, real and imaginary part of a complex,
        ...
        """
        return self.components is not None

    def get_components(self):
        """Return the list of components of the field

        Components (or sub fields) are intended for non scalar fields (ex:
        vector like current or wind, real and imaginary part of a complex,...)
        """
        res = []
        if self.is_composite():
            for rec in self.components:
                res.extend(rec.get_components())
            return res
        else:
            return [self]

    def get_dimsize(self, dimname):
        return self.dimensions[dimname]

    def get_dimnames(self):
        return self.dimensions.keys()

    def attach_storage(self, handler):
        """Attach a mapper to the field.

        When a field is stored in a file, its data are not loaded until they
        are explicitly requested (call to ``get_values()``). The handler
        pointer is set by this function to locate where the data have to be
        read (or saved later for a new or updated file).

        Args:
            handler (AbstractMapper): the file mapper attached to the field.
        """
        self.handler = handler
        return

    def __fix_subset_offset(self, subset):
        """
        correct the slices of the subset where they are out of the field size
        """
        fixed_subset = []
        dim_names = self.dimensions.keys()
        for dim_i, dimslice in enumerate(subset):
            if not isinstance(dimslice, slice):
                raise Exception("Not a slice")
            size = self.dimensions[dim_names[dim_i]]
            if dimslice.start is not None:
                i0 = max(0, dimslice.start)
            else:
                i0 = None
            if dimslice.stop is not None:
                i1 = min(size, dimslice.stop)
            else:
                i1 = None
            fixed_subset.append(slice(i0, i1))
        return fixed_subset, subset

    def __pad_values(self, values, fixed_subset, original_subset):
        """
        pad the `values` array with fill values where subset indices are out
        of the field size.
        """
        padded_sizes = []
        padded_slices = []
        need_padding = False
        # Dimension adjust allows for ignoring if dimensions that are not
        # sliced (like a GHRSST time dimension)
        dimension_adjust = 0
        for dim_i, dimslice in enumerate(fixed_subset):
            if isinstance(dimslice, int):
                if dimslice == 0:
                    dimension_adjust += 1
                    continue
            padded_dimslice = original_subset[dim_i - dimension_adjust]
            if padded_dimslice.stop is not None and\
                    padded_dimslice.start is not None:
                padded_size = padded_dimslice.stop - padded_dimslice.start
                padded_sizes.append(padded_size)
                i0 = dimslice.start
                i1 = dimslice.stop
                if i0 != padded_dimslice.start or i1 != padded_dimslice.stop:
                    need_padding = True
                rel_i0 = i0 - padded_dimslice.start
                rel_i1 = padded_size - (padded_dimslice.stop - i1)
                padded_slices.append(slice(rel_i0, rel_i1))
            else:
                padded_sizes.append(values.shape[dim_i - dimension_adjust])
                padded_slices.append(padded_dimslice)
        if need_padding:
            new_values = numpy.ma.masked_all(
                                shape=tuple(padded_sizes),
                                dtype=values.dtype)
            if hasattr(values, 'fill_value'):
                if values.fill_value is not None:
                    new_values.set_fill_value(values.fill_value)
            new_values[padded_slices] = values[:]
            return new_values
        else:
            return values

    def get_values(self, slices=None, indices=None, cache=False, padding=False,
                   strong_cache=False, **kwargs):
        '''
        Args
            cache (bool): if cache is True, the data read from file are kept in
                memory. The full field data array is kept in cache, slices or
                indices are ignored by caching (though the result returned by
                this function will be a subset if slices or indices are
                provided. Use the `view` concept in mappers instead if you want
                frequent accesses to a subset of a file.
                Default is False.

            strong_cache (bool): *DEPRECATED* - same as `cache`

            padding (bool) : if True, pad the result with fill values where
                slices are out of the field size.
        '''
        if cache or strong_cache:
            # cache reads in the entire field, regardless of any
            # slicing, this is useful in some instances for lat / lon / time
            # fields, or when many many subsets are extracted from a netCDF
            # mapper with a Grid model.
            if self._values is None:
                self._values = self.handler.mapper.read_values(
                    self.handler.fieldname, None
                )
        if indices is not None and slices is not None:
            logging.error("Indices and slices provided : %s %s",
                          indices,
                          slices)
            raise Exception("slices and indices can not be provided at the"
                            " same time and are exclusive")
        if indices is not None:
            # transform indices
            subset = Field.format_indices(indices,
                                          self.dimensions)
        elif slices is not None:
            if not isinstance(slices, dict):
                raise Exception("slices must be a dictionary of python slice"
                                " objects.")
            subset = Field.format_slices(slices,
                                         self.dimensions)
        else:
            subset = None
        if subset is not None and padding:
            subset, original_subset = self.__fix_subset_offset(subset)
        if self._values is None:
            if self.handler and self.handler.is_saved():
                values = self.handler.mapper.read_values(
                    self.handler.fieldname,
                    subset,
                    **kwargs
                    )
                if padding:
                    values = self.__pad_values(values, subset, original_subset)
            else:
                logging.warning("Field %s has not yet any values",
                                self.get_name())
                return None
        else:
            if subset is None:
                values = self._values
            else:
                values = self._values[tuple(subset)]
                if padding:
                    values = self.__pad_values(values, subset, original_subset)
        if not isinstance(values, numpy.ndarray):
            # always return an array (even if one single element)
            values = numpy.array([values],)
        return values

    def set_values(self, values, slices=None, indices=None):
        if indices and slices:
            logging.error("Indices and slices provided : %s %s",
                          indices,
                          slices)
            raise Exception("slices and indices can not be \
                 provided at the same time and are exclusive")
        if indices:
            # transform indices
            subset = Field.format_indices(indices,
                                          self.dimensions)
        elif slices:
            subset = Field.format_slices(slices,
                                         self.dimensions)
        else:
            subset = None
        if subset is None:
            self._values = values
        else:
            self._values[tuple(subset)] = values
        if self.handler:
            self.handler.set_unsaved()
        return

    def is_saved(self):
        """
        Return True is the content of the field is saved on file and
        was not updated since
        """
        if self.handler is None:
            return False
        return self.handler.is_saved()

    @classmethod
    def compute(cls, operator, field1, field2=None, variable=None):
        """Perform an operation and returns the result as a field

        The operator may be for instance a numpy MaskedArray operator
        such as numpy.ma.anom, numpy.ma.corr,...

        To be used with caution.

        Args:
            operator (function) : the function to be called (ex: numpy.ma.anom)
            field1 (Field) : the field argument to the operator
            field2 (Field) : an optional 2nd field argument to the operator
            variable (Variable) : variable of the returned module field. If not
                provided, the returned field is created with a basic variable
                definition.

        Returns:
            Field: the result field
        """
        if variable is None:
            varname = 'result'
            variable = Variable(varname)
        if field2 is None:
            values = partial(operator(field1))
        else:
            values = partial(operator(field1, field2))
        field = Field(variable,
                      dimensions=copy.copy(field1.dimensions),
                      datatype=field1.datatype,
                      fillvalue=field1.fillvalue,
                      values=values,
                      units=field1.units)
        return field

    def clone(self):
        """Return a copy of the field without any data or reference to the
        source file"""
        field = Field(copy.copy(self.variable),
                      copy.copy(self.dimensions),
                      datatype=self.datatype,
                      fields=copy.copy(self.components),
                      fillvalue=self.fillvalue,
                      attributes=copy.copy(self.attributes),
                      units=self.units,
                      valid_min=self.valid_min,
                      valid_max=self.valid_max
                      )
        return field


def module(u, v, variable=None):
    """Return the module field from its two components

    The module is sqrt(u² + v²)

    Args:
        u (Field) : the eastward component
        v (Field) : the northward component
        variable (Variable) : variable of the returned module field. If not
            provided, the returned field is created with a basic variable
            definition.

    Returns:
        Field: the module field
    """
    values = numpy.ma.sqrt(numpy.ma.power(u.get_values(), 2)
                           + numpy.ma.power(v.get_values(), 2)
                           )
    if variable is None:
        if 'eastward_' in u.variable.shortname:
            varname = u.variable.shortname.replace('eastward_', '')
        else:
            varname = 'module'
        variable = Variable(varname)
    field = Field(variable,
                  dimensions=copy.copy(u.dimensions),
                  datatype=u.datatype,
                  fillvalue=u.fillvalue,
                  values=values,
                  units=u.units)
    return field
