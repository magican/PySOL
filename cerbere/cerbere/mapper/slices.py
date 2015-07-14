# -*- coding: utf-8 -*-
"""
.. module::cerbere.mapper.slices

Functions to handle slices when reading data from a :package:`mapper` or
:package:`datamodel` class.

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""


def get_shape_from_slice(slices):
    """Returns the shape of the result from a slicing operation.

    Args:
        slices (list<slice>): a list of slices

    Returns:
        tuple<int>: the shape of the expected result
    """
    shape = []
    for slc in slices:
        shape.append((slc.stop - slc.start) / slc.step)
    return tuple(shape)

def get_nice_slices(slices, dimsizes):
    """Return filled slices with valid positive limits.
    
    Fill (replace None values with numbers). Trim the slices beyond the array
    dimensions. Transform negative slices into positive values.
    """
    newslices = []
    if slices is None:
        for d in dimsizes:
            newslices.append(slice(0, d, 1))
        return newslices
    for sli, dimsize in zip(slices, dimsizes):
        start, stop = sli.start, sli.stop
        if start is not None and start < 0:
            start += dimsize
        elif start is None:
            start = 0
        if stop is not None and stop < 0:
            stop += dimsize
        elif stop is None:
            stop = dimsize
        if start is not None and start > dimsize:
            start = dimsize
        if stop is not None and stop > dimsize:
            stop = dimsize
        step = sli.step
        if step is None:
            step = 1
        newslices.append(slice(start, stop, step))   
    return newslices    
    

def fill_slices(slices, dimsizes):
    """Fill (replace None values with numbers) and check slices with respect
    to the dimensions of an array.    
    """
    # If you change the policy here, please make sure it is consistent
    # with slices management in read_values() and get_dimsize()
    if len(slices) != len(dimsizes):
        raise Exception('Unexpected slices list length : {}'
                        .format(slices))
    filled_slices = []
    for sli, dimsize in zip(slices, dimsizes):
        if sli.start is not None and abs(sli.start) > dimsize:
            raise Exception('Unexpected slice start : {}'.format(sli))
        if sli.stop is not None and abs(sli.stop) > dimsize:
            raise Exception('Unexpected slice stop : {}'.format(sli))
        start, stop, step = sli.indices(dimsize)
        if (step > 0 and start >= stop) or (step < 0 and start <= stop):
            raise Exception('Unexpected slice start Vs stop : {}'
                            .format(sli))
        # Special case : negative step until the first element
        if stop == -1:
            stop = None
        filled_slices.append(slice(start, stop, step))
    return filled_slices


def get_absolute_slices(view, slices, dimnames, dimsizes):
    """Returns absolute slices from the slices relative to a view."""
    # view slicing
    if view is not None:
        # make sure viewslice info is complete
        viewslices = format_slices(view, dimnames)
        viewslices = fill_slices(viewslices, dimsizes)
    # standard slicing
    if slices is not None:
        if view is None:
            finalslices = fill_slices(slices, dimsizes)
        else:
            # slices are relative to the opened view
            viewdimsizes = dimsizes
            finalslices = fill_slices(slices, viewdimsizes)
            for index, (sli, vwsli) in enumerate(zip(finalslices,
                                                     viewslices)):
                vw_first = vwsli.start
                first = sli.start
                if sli.stop is None:
                    last = 0
                else:
                    last = sli.stop - sli.step / abs(sli.step)
                first = vw_first + first * vwsli.step
                last = vw_first + last * vwsli.step
                start = first
                step = vwsli.step * sli.step
                stop = last + step / abs(step)
                if stop == -1:
                    stop = None
                finalslices[index] = slice(start, stop, step)
    elif view is not None:
        finalslices = viewslices
    else:
        finalslices = [slice(None) for i in range(len(dimnames))]
        finalslices = fill_slices(finalslices, dimsizes)
    return finalslices


def adjust_dimsize(view, dimname, dimsize):
    """correct the size of a dimension with respect to a view"""
    if view is not None and dimname in view:
        viewslice = fill_slices([view[dimname]],
                                [dimsize])[0]
        start, stop, step = viewslice.start, viewslice.stop, viewslice.step
        if stop is None:
            stop = -1
        dimsize = 1 + (abs(stop - start) - 1) / abs(step)
    return dimsize


def format_slices(slices, fielddims):
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
    if slices is None:
        return None
    formatted_slices = []
    dimensions = fielddims
    for _ in xrange(len(dimensions)):
        formatted_slices.append(slice(None, None, None))
    for k in slices.keys():
        try:
            idim = list(dimensions).index(k)
            formatted_slices[idim] = slices[k]
        except ValueError:
            pass
    return formatted_slices


def format_indices(indices, fielddims):
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
    if indices is None:
        return None
    formatted_indices = []
    dimensions = fielddims
    for i in xrange(len(dimensions)):
        formatted_indices.append(slice(None, None, None))
    for k in indices.keys():
        try:
            idim = list(dimensions).index(k)
            formatted_indices[idim] = indices[k]
        except ValueError:
            pass
    return tuple(formatted_indices)
