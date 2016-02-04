#!/usr/bin/env python
# coding=utf-8
"""
"""


import os
from collections import OrderedDict
import xml.etree.ElementTree as ET
import numpy as np
import pyproj
import glob
from netCDF4 import date2num

from cerbere.mapper.abstractmapper import AbstractMapper
from cerbere import READ_ONLY
from cerbere.mapper.slices import Slices
from cerbere.mapper.safemsil1cgranulefile import SAFEMSIL1CGranuleFile


def safemsil1c_stitching_groups(safe_path):
    """
    """
    groups = OrderedDict()
    gran_paths = glob.glob(os.path.join(safe_path, 'GRANULE', '*'))
    gran_paths.sort()
    for gran_path in gran_paths:
        gran_id = os.path.basename(gran_path)
        granmtd_id = '{}.xml'.format(gran_id[:-7].replace('MSI', 'MTD'))
        granmtd_path = os.path.join(gran_path, granmtd_id)
        tree = ET.parse(granmtd_path)
        root = tree.getroot()
        xmlns = {'n1': 'https://psd-12.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd',
                 'xsi': 'http://www.w3.org/2001/XMLSchema-instance'}
        proj_code = root.find('n1:Geometric_Info/Tile_Geocoding/HORIZONTAL_CS_CODE', xmlns).text
        if proj_code in groups:
            groups[proj_code].append(os.path.basename(gran_path))
        else:
            groups[proj_code] = [os.path.basename(gran_path)]
    return groups


class SAFEMSIL1CStitchedFile(AbstractMapper):
    """Abstract class for SAFE MSI L1C stitched files.

    Args:
        url (str): SAFE directory.

        granules (list): list of granules ID.

        native_resolution (str): 10m, 20m or 60m. The mapper exhibits only
            the bands relative to this native resolution.

        overview_index (int, optional): if set, the mapper will use the
            corresponding overview contained in JPEG 2000 files and will
            exhibits dimensions and geolocation for this overview. Index starts
            at 0 and ends at the number of overviews minus 1. To be used for
            quicklooks as it will be more noisy than averaging the full
            resolution. With None (default value), the full resolution is used.

        tight (bool, optional): if set, the mapper will exhibit the smallest
            bounding box surrounding the valid data. The default value is False.
    """
    def __init__(self, url, granules, native_resolution,
                 overview_index=None, tight=False, mode=READ_ONLY, **kwargs):
        """
        """
        if 'PRD_MSIL1C' not in os.path.basename(url):
            raise Exception('MSI L1C SAFE expected.')
        if native_resolution not in ['10m', '20m', '60m']:
            raise Exception('Wrong MSI L1C native resolution.')
        if mode != READ_ONLY:
            raise Exception('This mapper can only be used in read_only mode.')
        super(SAFEMSIL1CStitchedFile, self).__init__(url=url, mode=mode,
                                                     **kwargs)
        self._granules = granules
        self._native_resolution = native_resolution
        self._overview_index = overview_index
        self._tight = tight
        self._mappers = []
        self._offsets = []
        self._sizes = []
        self._attrs = OrderedDict()
        self._dims = OrderedDict()
        self._fields = OrderedDict()

    def get_geolocation_field(self, fieldname):
        """Return the equivalent field name in the file format for a standard
        geolocation field (lat, lon, time, z).

        Used for internal purpose and should not be called directly.

        Args:
            fieldname (str): name of the standard geolocation field (lat, lon
                or time)

        Return:
            str: name of the corresponding field in the native file format.
                Returns None if no matching is found
        """
        return fieldname

    def get_matching_dimname(self, dimname):
        """Return the equivalent name in the native format for a standard
        dimension.

        This is a translation of the standard names to native ones. It is used
        for internal purpose only and should not be called directly.

        The standard dimension names are:

        * x, y, time for :class:`~cerbere.datamodel.grid.Grid`
        * row, cell, time for :class:`~cerbere.datamodel.swath.Swath` or
          :class:`~cerbere.datamodel.image.Image`

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        Args:
            dimname (str): standard dimension name.

        Returns:
            str: return the native name for the dimension. Return `dimname` if
                the input dimension has no standard name.

        See Also:
            see :func:`get_standard_dimname` for the reverse operation
        """
        return dimname

    def get_standard_dimname(self, dimname):
        """
        Returns the equivalent standard dimension name for a
        dimension in the native format.

        This is a translation of the native names to standard ones. It is used
        for internal purpose and should not be called directly.

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        Args:
            dimname (string): native dimension name

        Return:
            str: the (translated) standard name for the dimension. Return
            `dimname` if the input dimension has no standard name.

        See Also:
            see :func:`get_matching_dimname` for the reverse operation
        """
        return dimname

    def open(self, datamodel=None, datamodel_geolocation_dims=None):
        """Open the file (or any other type of storage)

        Args:
            datamodel (str): type of feature read or written. Internal argument
                only used by the classes from :mod:`~cerbere.datamodel`
                package. Can be 'Grid', 'Swath', etc...

            datamodel_geolocation_dims (list, optional): list of the name of the
                geolocation dimensions defining the data model to be read in
                the file. Optional argument, only used by the datamodel
                classes, in case the mapper class can store different types of
                data models.

        Returns:
            a handler on the opened file
        """
        if self.is_opened():
            return self._handler
        # Get and open granule mappers
        for granule in self._granules:
            url = os.path.join(self.get_url(), 'GRANULE', granule)
            mapper = SAFEMSIL1CGranuleFile(url, self._native_resolution,
                                           overview_index=self._overview_index,
                                           tight=self._tight)
            mapper.open()
            self._mappers.append(mapper)
        # Define arbitrarily the first granule handler as the mapper handler
        handler = self._mappers[0]._handler
        self._handler = handler
        if datamodel is not None:
            self._feature_type = datamodel
        if datamodel_geolocation_dims is not None:
            self.datamodel_geolocation_dims = datamodel_geolocation_dims
        # Stitching grid
        code = [m.read_global_attribute('horizontal_cs_code') \
                for m in self._mappers]
        if min(code) != max(code):
            raise Exception('Granules do not share the same projection')
        gxdim = [m.read_global_attribute('xdim') for m in self._mappers]
        gydim = [m.read_global_attribute('ydim') for m in self._mappers]
        if min(gxdim) != max(gxdim) or min(gydim) != max(gydim):
            raise Exception('Granules do not share the same spacing.')
        xdim = gxdim[0]
        ydim = gydim[0]
        gnrows = [m.read_global_attribute('nrows') for m in self._mappers]
        gncols = [m.read_global_attribute('ncols') for m in self._mappers]
        gulx = [m.read_global_attribute('ulx') for m in self._mappers]
        guly = [m.read_global_attribute('uly') for m in self._mappers]
        glrx = [x0 + xdim * nx for x0, nx in zip(gulx, gnrows)]
        glry = [y0 + ydim * ny for y0, ny in zip(guly, gncols)]
        if xdim > 0:
            ulx = min(gulx)
            lrx = max(glrx)
        else:
            ulx = max(gulx)
            lrx = min(glrx)
        if ydim > 0:
            uly = min(guly)
            lry = max(glry)
        else:
            uly = max(guly)
            lry = min(glry)
        for x0, y0, nx, ny in zip(gulx, guly, gnrows, gncols):
            xoffset = int((x0 - ulx) / xdim)
            yoffset = int((y0 - uly) / ydim)
            self._offsets.append([yoffset, xoffset])
            self._sizes.append([ny, nx])
        nrows = int((lrx - ulx) / xdim)
        ncols = int((lry - uly) / ydim)
        # Construct attributes
        self._attrs = self._mappers[0]._attrs.copy()
        attr2list = ['tile_id', 'sensing_time', 'archiving_centre',
                     'archiving_time']
        for att in attr2list:
            self._attrs[att] = [m.read_global_attribute(att) \
                                for m in self._mappers]
        self._attrs['nrows'] = nrows
        self._attrs['ncols'] = ncols
        self._attrs['ulx'] = ulx
        self._attrs['uly'] = uly
        # Construct dimensions
        self._dims['time'] = 1
        self._dims['y'] = self._attrs['ncols']
        self._dims['x'] = self._attrs['nrows']
        # Construct fields
        for fn, f in self._mappers[0]._fields.iteritems():
            field = f.clone()
            for dimname in field.dimensions.keys():
                if dimname in self._dims.keys():
                    field.dimensions[dimname] = self._dims[dimname]
            self._fields[fn] = field
        for fname in self._fields:
            self._fields[fname].attach_storage(self.get_field_handler(fname))
        return handler

    def close(self):
        """Close handler on storage"""
        for m in self._mappers:
            m.close()
        self._mappers = []
        self._offsets = []
        self._sizes = []
        self._handler = None

    def read_values(self, fieldname, slices=None):
        """Read the data of a field.

        Args:
            fieldname (str): name of the field which to read the data from

            slices (list of slice, optional): list of slices for the field if
                subsetting is requested. A slice must then be provided for each
                field dimension. The slices are relative to the opened view (see
                :func:open) if a view was set when opening the file.

        Return:
            MaskedArray: array of data read. Array type is the same as the
                storage type.
        """
        field = self.read_field(fieldname)
        dimsizes = field.dimensions.values()
        slices = Slices(slices, dimsizes)
        if fieldname == 'time':
            start_time = date2num(self.get_start_time(), field.units)
            end_time = date2num(self.get_end_time(), field.units)
            time = (start_time + end_time) / 2
            values = np.ma.array([time], dtype=field.datatype)[slices]
        elif fieldname in ['lon', 'lat']:
            epsg = self.read_global_attribute('horizontal_cs_code')
            proj = pyproj.Proj(init=epsg)
            shp = slices.shape()
            y = np.tile(self.read_values('y', slices=[slices[0]])[:, np.newaxis],
                        (1, shp[1]))
            x = np.tile(self.read_values('x', slices=[slices[1]])[np.newaxis, :],
                        (shp[0], 1))
            lon, lat = proj(x, y, inverse=True)
            if fieldname == 'lon':
                values = lon.astype(field.datatype)
            elif fieldname == 'lat':
                values = lat.astype(field.datatype)
        elif fieldname in ['x', 'y']:
            if fieldname == 'x':
                spacing = self.read_global_attribute('xdim')
                start = self.read_global_attribute('ulx') + spacing / 2.
            elif fieldname == 'y':
                spacing = self.read_global_attribute('ydim')
                start = self.read_global_attribute('uly') + spacing / 2.
            inds = slices.indices_array(dtype=field.datatype)[0] * spacing + start
            values = np.ma.array(inds, dtype=field.datatype)
        else:
            shp = slices.shape()
            vals = np.zeros(shp, dtype=field.datatype)
            if field.datatype != np.dtype(np.bool):
                fillvalue = field.fillvalue
                count = np.zeros(shp, dtype='uint8')
                for m, off, siz in zip(self._mappers, self._offsets, self._sizes):
                    outsub, insub = slices.outsub_insub_slices(off, siz)
                    if outsub is None:
                        continue
                    subvals = m.read_values(fieldname, slices=insub)
                    if fillvalue is not None and fillvalue != 0 and \
                       subvals.mask is not np.ma.nomask:
                        subvals.data[subvals.mask] = 0
                    vals[outsub] += subvals
                    count[outsub] += np.uint8(1) - np.ma.getmask(subvals)
                over = np.where(count > 1)
                vals[over] /= count[over]
                if count.min() == 0:
                    mask = count == 0
                    if fillvalue is not None and fillvalue != 0:
                        vals[mask] = fillvalue
                    values = np.ma.array(vals, mask=mask)
                else:
                    values = np.ma.array(vals)
            else:
                for m, off, siz in zip(self._mappers, self._offsets, self._sizes):
                    outsub, insub = slices.outsub_insub_slices(off, siz)
                    if outsub is None:
                        continue
                    subvals = m.read_values(fieldname, slices=insub)
                    vals[outsub] += subvals
                values = vals
        return values

    def read_field(self, fieldname):
        """
        Return the :class:`cerbere.field.Field` object corresponding to
        the requested fieldname.

        The :class:`cerbere.field.Field` class contains all the metadata
        describing a field (equivalent to a variable in netCDF).

        Args:
            fieldname (str): name of the field

        Returns:
            :class:`cerbere.field.Field`: the corresponding field object
        """
        return self._fields[fieldname]

    def write_field(self, fieldname):
        """Writes the field data on disk.

        Args:
            fieldname (str): name of the field to write.
        """
        raise NotImplementedError

    def read_fillvalue(self, fieldname):
        """Read the fill value of a field.

        Args:
            fieldname (str): name of the field.

        Returns:
            number or char or str: fill value of the field. The type is the
                as the type of the data in the field.
        """
        return self.read_field(fieldname).fillvalue

    def create_field(self, field, dim_translation=None):
        """Creates a new field in the mapper.

        Creates the field structure but don't write yet its values array.

        Args:
            field (Field): the field to be created.

        See also:
            :func:`write_field` for writing the values array.
        """
        raise NotImplementedError

    def create_dim(self, dimname, size=None):
        """Add a new dimension.

        Args:
            dimname (str): name of the dimension.
            size (int): size of the dimension (unlimited if None)
        """
        raise NotImplementedError

    def get_fieldnames(self):
        """Returns the list of geophysical fields stored for the feature.

        The geolocation field names are excluded from this list.

        Returns:
            list<string>: list of field names
        """
        fieldnames = self._fields.keys()
        fieldnames.remove('time')
        fieldnames.remove('lon')
        fieldnames.remove('lat')
        return fieldnames

    def get_dimensions(self, fieldname=None):
        """Return the dimension's standard names of a file or a field in the
        file.

        Args:
            fieldname (str): the name of the field from which to get the
                dimensions. For a geolocation field, use the cerbere standard
                name (time, lat, lon), though native field name will work too.

        Returns:
            tuple<str>: the standard dimensions of the field or file.
        """
        if fieldname is None:
            dims = self._dims.keys()
        else:
            dims = self.read_field(fieldname).dimensions.keys()
        return tuple(dims)

    def get_dimsize(self, dimname):
        """Return the size of a dimension.

        Args:
            dimname (str): name of the dimension.

        Returns:
            int: size of the dimension.
        """
        return self._dims[dimname]

    def read_field_attributes(self, fieldname):
        """Return the specific attributes of a field.

        Args:
            fieldname (str): name of the field.

        Returns:
            dict<string, string or number or datetime>: a dictionary where keys
                are the attribute names.
        """
        return self.read_field(fieldname).attributes

    def read_global_attributes(self):
        """Returns the names of the global attributes.

        Returns:
            list<str>: the list of the attribute names.
        """
        return self._attrs.keys()

    def write_global_attributes(self, attrs):
        """Write the global attributes of the file.

        Args:
            attrs (dict<string, string or number or datetime>): a dictionary
                containing the attributes names and values to be written.
        """
        raise NotImplementedError

    def read_global_attribute(self, name):
        """Returns the value of a global attribute.

        Args:
            name (str): name of the global attribute.

        Returns:
            str, number or datetime: value of the corresponding attribute.
        """
        return self._attrs[name]

    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage.

        Returns:
            datetime: start time of the data in file.
        """
        return min([m.get_start_time() for m in self._mappers])

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage.

        Returns:
            datetime: end time of the data in file.
        """
        return max([m.get_end_time() for m in self._mappers])

    def get_bbox(self):
        """Returns the bounding box of the feature, as a tuple.

        Returns:
            tuple: bbox expressed as (lonmin, latmin, lonmax, latmax)
        """
        dimsizes = self.get_full_dimensions('lon').values()
        slices = [slice(None, None, dimsizes[0] - 1),
                  slice(None, None, dimsizes[1] - 1)]
        lon = self.read_values('lon', slices=slices)
        lat = self.read_values('lat', slices=slices)
        return (lon.min(), lat.min(), lon.max(), lat.max())

    def get_orbit_number(self):
        """In the case of a satellite orbit file, returns the orbit number.

        Returns:
            int: the orbit number
        """
        return self.read_global_attribute('sensing_orbit_number')


class SAFEMSIL1C10mStitchedFile(SAFEMSIL1CStitchedFile):
    """Mapper class for SAFE MSI L1C 10m stitched files.
    """
    def __init__(self, url, granules, **kwargs):
        super(SAFEMSIL1C10mStitchedFile, self).__init__(url, granules, '10m',
                                                        **kwargs)


class SAFEMSIL1C20mStitchedFile(SAFEMSIL1CStitchedFile):
    """Mapper class for SAFE MSI L1C 20m stitched files.
    """
    def __init__(self, url, granules, **kwargs):
        super(SAFEMSIL1C20mStitchedFile, self).__init__(url, granules, '20m',
                                                        **kwargs)


class SAFEMSIL1C60mStitchedFile(SAFEMSIL1CStitchedFile):
    """Mapper class for SAFE MSI L1C 60m stitched files.
    """
    def __init__(self, url, granules, **kwargs):
        super(SAFEMSIL1C60mStitchedFile, self).__init__(url, granules, '60m',
                                                        **kwargs)

