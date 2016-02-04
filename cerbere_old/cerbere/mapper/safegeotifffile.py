# encoding: utf-8
"""
cerbere.mapper.safegeotifffile
==============================

Mapper class for ESA SAFE GeoTiff files (for Sentinel-1)

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""


import os
import logging
from collections import OrderedDict
import xml.etree.ElementTree as ET
from datetime import datetime
import numpy as np
# NOTE : in scientific_toolbox_cloudphys_precise "virtual env",
# importing gdal before netCDF4 makes a 'netCDF4.Dataset(url)' crashes.
# Don't modify the order of the next three lines.
from netCDF4 import date2num, num2date
import gdal
from gdalconst import GA_ReadOnly
# \NOTE
from scipy import interpolate
import warnings

from .. import READ_ONLY
from .abstractmapper import AbstractMapper
from ..datamodel.field import Field
from ..datamodel.variable import Variable


# See `codes` variable in osgeo.gdal_array
GDAL2NP = {'Byte': np.dtype(np.uint8),
           'UInt16': np.dtype(np.uint16),
           'Int16': np.dtype(np.int16),
           'UInt32': np.dtype(np.uint32),
           'Int32': np.dtype(np.int32),
           'Float32': np.dtype(np.float32),
           'Float64': np.dtype(np.float64),
           'CInt16': np.dtype(np.complex64), # S1 SLC external
           'CInt32': np.dtype(np.complex64),
           'CFloat32': np.dtype(np.complex64), # S1 SLC internal
           'CFloat64': np.dtype(np.complex128)}
EARTHRADIUS = 6371009. # Mean earth radius [m]


class SAFEGeoTiffFile(AbstractMapper):
    """
    Generic storage class for Sentinel GeoTiff files.

    S-1 GeoTiff files are actually provided with additional metadata to be
    read in joined XML files (in `annotation` folder).
    """

    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        AbstractMapper.__init__(self, url=url, mode=mode, **kwargs)
        self._attributes = OrderedDict()
        self._lon = None
        self._lat = None
        self._lonlat_slices = None
        self.open()
        return

    def open(self, datamodel_geolocation_dims=None):
        """
        Open SAFE GeoTiff file and read metadata.
        """
        if self.is_opened():
            return self._handler
        if self.is_writable():
            raise NotImplementedError
        else:
            if not os.path.exists(self._url):
                raise Exception("File %s is not existing" % self._url)
        if (self._url is not None) and (self._mode is not None):
#             logging.debug("MODE : %s", self._mode)
            self._handler = gdal.Open(self._url, GA_ReadOnly)
            self._attributes['safe_name'] = self._get_safe_name()
            self._attributes.update(self._read_annotationdatasets())
            return self._handler
        else:
            return None

    def close(self):
        """
        Close SAFE GeoTiff file.
        """
        self._handler = None

    def get_geolocation_field(self, fieldname):
        """
        Return the equivalent field name in the file format for a standard
        geolocation field (lat, lon, time).

        Used for internal purpose and should not be called directly.

        :param fieldname: name of the standard geolocation field (lat, lon
            or time)
        :type fieldname: str

        :rtype: str or None
        :return: name of the corresponding field in the native file format.
            Returns None if no matching is found
        """
        return fieldname

    def get_matching_dimname(self, geodimname):
        """
        Return the equivalent name in the native format for a standard
        dimension.

        This is a translation of the standard names to native ones. It is used
        for internal purpose only and should not be called directly.

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        :arg dimname: standard dimension name
        :type dimname: str

        :rtype: str
        :return: return the native name for the dimension.
            return `dimname` if the input dimension has no standard name.

        *See*

        see :func:`get_standard_dimname` for the reverse operation
        """
        return geodimname

    def get_standard_dimname(self, geodimname):
        """
        Return the equivalent standard dimension name for a
        dimension in the native format.

        This is a translation of the native names to standard ones. It is used
        for internal purpose and should not be called directly.

        To be derived when creating an inherited data mapper class. This is
        mandatory for geolocation dimensions which must be standard.

        :arg dimname: native dimension name
        :type dimname: str

        :rtype: str
        :return: return the (translated) standard name for the dimension.
            return `dimname` if the input dimension has no standard name.

        *See*

        see :func:`get_matching_dimname` for the reverse operation
        """
        return geodimname

    def get_fieldnames(self):
        """
        Return the list of geophysical fields stored for the feature.

        :return: list of field names
        :rtype: list<str>
        """
        return ['time', 'lon', 'lat', 'incidence', 'elevation',
                'digital_number', 'gamma0_lut', 'beta0_lut', 'sigma0_lut',
                'gamma0', 'beta0', 'sigma0', 'complex', 'azimuth_time',
                'slant_range_time', 'doppler_centroid']

    def read_field_attributes(self, fieldname):
        """
        Return the specific storage attributes of a field
        (such as _FillValue, scale_factor, add_offset,...).

        :param fieldname: name of the field
        :type fieldname: str

        :return: a dictionary where keys are the attribute names.
        :rtype: Dict
        """
        attrs = {}
        return attrs

    def get_dimsize(self, dimname):
        """
        """
        if dimname == 'time':
            return 1
        elif dimname == 'row':
            return self.get_handler().RasterYSize
        elif dimname == 'cell':
            return self.get_handler().RasterXSize
        else:
            raise Exception('Unknown dimname : %s' % dimname)

    def get_dimensions(self, fieldname=None):
        """
        Return the dimension names of a file or a field in the file.

        :keyword fieldname: the field from which to get the dimension names.
            For a geolocation field, use the cerbere standard name
            (time, lat, lon), though native field name will work too.
        :type fieldname: str

        :return: the dimensions of the field or file.
        :rtype: tuple of strings
        """
        if fieldname is None:
            return ('time', 'row', 'cell')
        elif fieldname == 'time':
            return ('time',)
        elif fieldname in self.get_fieldnames():
            return ('row', 'cell')
        else:
            raise Exception('Unknown fieldname : %s' % fieldname)

    def get_full_dimensions(self, fieldname):
        """
        Return the dimension names and sizes of a field.

        :param fieldname: name of the field
        :type fieldname: str

        :rtype: OrderedDict
        :return: Ordered dictionary where keys are the dimension names
            and values are their respective size
        """
        dims = OrderedDict()
        if fieldname == 'time':
            dims['time'] = self.get_dimsize('time')
        elif fieldname in self.get_fieldnames():
            dims['row'] = self.get_dimsize('row')
            dims['cell'] = self.get_dimsize('cell')
        else:
            raise Exception('Unknown fieldname : %s' % fieldname)
        return dims

    def read_field(self, fieldname):
        """
        Return the :class:`cerbere.field.Field` object corresponding to
        the requested fieldname.

        The :class:`cerbere.field.Field` class contains all the metadata
        describing a field (equivalent to a variable in netCDF).

        :param fieldname: name of the field
        :type fieldname: str

        :return: the corresponding field object
        :rtype: :class:`cerbere.field.Field`
        """
        if fieldname not in self.get_fieldnames():
            raise Exception('Unknown fieldname : %s' % fieldname)
        # Variable defaults
        shortname = fieldname
        description = None
        authority = None
        standardname = None
        # Field defaults
        dims = self.get_full_dimensions(fieldname)
        datatype = None
        units = None
        # Create virtual fields
        if fieldname == 'time':
            description = 'time of measurement'
            datatype = float
            units = 'milliseconds since 1990-01-01T00:00:00'
        elif fieldname == 'lat':
            description = 'latitude'
            standardname = 'latitude'
            datatype = np.dtype(np.float32)
            units = 'degrees_north'
        elif fieldname == 'lon':
            description = 'longitude'
            standardname = 'longitude'
            datatype = np.dtype(np.float32)
            units = 'degrees_east'
        elif fieldname == 'incidence':
            description = 'incidence angle'
            datatype = np.dtype(np.float32)
            units = 'degrees'
        elif fieldname == 'elevation':
            description = 'elevation angle'
            datatype = np.dtype(np.float32)
            units = 'degrees'
        elif fieldname == 'digital_number':
            description = 'digital number'
            band = self.get_handler().GetRasterBand(1)
            typestr = gdal.GetDataTypeName(band.DataType)
            datatype = GDAL2NP[typestr]
        elif fieldname == 'gamma0_lut':
            description = 'interpolated gamma0 lookup table'
            datatype = np.dtype(np.float32)
        elif fieldname == 'beta0_lut':
            description = 'interpolated beta0 lookup table'
            datatype = np.dtype(np.float32)
        elif fieldname == 'sigma0_lut':
            description = 'interpolated sigma0 lookup table'
            datatype = np.dtype(np.float32)
        elif fieldname == 'gamma0':
            description = 'gamma0'
            datatype = np.dtype(np.float32)
        elif fieldname == 'beta0':
            description = 'beta0'
            datatype = np.dtype(np.float32)
        elif fieldname == 'sigma0':
            description = 'sigma0'
            datatype = np.dtype(np.float32)
        elif fieldname == 'complex':
            description = 'complex intensity'
        elif fieldname == 'azimuth_time':
            description = 'zero doppler azimuth time'
            units = 'milliseconds since 1990-01-01T00:00:00'
            datatype = np.dtype(np.float64)
        elif fieldname == 'slant_range_time':
            description = 'two way slant range time'
            units = 'seconds'
            datatype = np.dtype(np.float64)
        elif fieldname == 'doppler_centroid':
            description = 'doppler centroid'
            units = 'Hz'
            datatype = np.dtype(np.float32)
        else:
            raise NotImplementedError
        variable = Variable(shortname=shortname,
                            description=description,
                            authority=authority,
                            standardname=standardname)
        field = Field(variable, dims, datatype=datatype, units=units)
        field.attach_storage(self.get_field_handler(fieldname))
        return field

    def read_values(self, fieldname, slices=None, blocksize=None):
        """
        Read the values of a field.

        `slices` is optional. When provided, give for each dimension the
        corresponding python slice object to subset this dimension. Only the
        dimensions to be subsetted need to be provided, by default the full
        dimension length is read.

        .. code-block:: python

            # extracting a subset of a field with slices
            data = fd.read_values('sigma0', slices=[slice(10,20), slice(30,40)])

        :param fieldname: name of the field
        :type fieldname: str
        :param slices: dimensions slices, when reading a subset only for some
            dimensions
        :type slices: List<slice>
        """
        # Check fieldname
        if fieldname not in self.get_fieldnames():
            raise Exception('Unknown fieldname : %s' % fieldname)
        # Treat time special case
        if fieldname == 'time':
            time_units = self.read_field('time').units
            values = (date2num(self.get_start_time(), time_units)+ \
                      date2num(self.get_end_time(), time_units))/2.
            return values
        # Format slices
        dims = self.get_full_dimensions(fieldname)
        if dims.keys() != ['row', 'cell']:
            raise Exception('Unexpected dimensions')
        if slices is None:
            slices = [slice(0, dims['row'], 1), slice(0, dims['cell'], 1)]
        else:
            if isinstance(slices, list) == False:
                raise Exception('slices has to be a list of two slice objects')
            if len(slices) != 2:
                raise Exception('slices has to be a list of two slice objects')
            isslice = all([isinstance(slices[isl], slice) for isl in range(2)])
            if isslice == False:
                raise Exception('slices has to be a list of two slice objects')
            slices = list(slices) # copy
            for isl, dim in zip(range(2), dims.values()):
                (start, stop, step) = slices[isl].indices(dim)
                slices[isl] = slice(start, stop, step)
        # Read by blocks for memory issues
        if blocksize is not None:
            (az0, az1, azs) = slices[0].indices(dims['row'])
            azsize = np.ceil((az1-az0)/float(azs))
            (ra0, ra1, ras) = slices[1].indices(dims['cell'])
            rasize = np.ceil((ra1-ra0)/float(ras))
            if azsize*rasize > blocksize:
                blrasize = rasize
                blazsize = np.maximum(blocksize//blrasize, 1)
                nblocks = np.ceil(azsize/float(blazsize))
                lims = np.round(np.linspace(0, azsize, num=nblocks+1))
                blslices = list(slices) # copy
                for ibl in np.arange(nblocks):
                    (lim0, lim1) = (lims[ibl], lims[ibl+1]-1)
                    blaz0 = int(az0+lim0*azs)
                    blaz1 = int(az0+lim1*azs+1)
                    blslices[0] = slice(blaz0, blaz1, azs)
                    blvalues = self.read_values(fieldname, slices=blslices,
                                                blocksize=None)
                    if ibl == 0:
                        dtype = blvalues.dtype
                        values = np.zeros((azsize, rasize), dtype=dtype)
                    values[lim0:lim1+1, :] = blvalues
                return values
        # Read/interpolate values according to fieldname
        if fieldname in ['lon', 'lat']:
            line = np.mgrid[slices[0]].astype('int32')
            pixel = np.mgrid[slices[1]].astype('int32')
            nbursts = self.read_global_attribute('number_of_bursts')
            if nbursts != 0:
                geoloc = self._get_geolocation_grid()
                geo_line = np.copy(geoloc['line'][:, geoloc['npixels']/2])
                geo_line_time = np.copy(geoloc['azimuth_time'][:, geoloc['npixels']/2])
                azitimeint = (geo_line_time[-1]-geo_line_time[-2]) / \
                             (geo_line[-1]-geo_line[-2])
                lines_per_burst = self.read_global_attribute('lines_per_burst')
                if geo_line.size == nbursts+1: # General case
                    # Check if geo_line is as expected
                    geo_line_exp = np.arange(nbursts+1)*lines_per_burst
                    geo_line_exp[-1] -= 1
                    if (geo_line[0:-1] == geo_line_exp[0:-1]).all() == False:
                        raise Exception('Unexpected geolocation grid !')
                    if geo_line[-1] == geo_line_exp[-1]:
                        pass
                    elif (geo_line[-1] - 1) == geo_line_exp[-1]:
                        # Get correct azimuth time interval
                        # (no way to guess it from geolocation grid)
                        azitimeint = self.read_global_attribute('azimuth_time_interval') * 1000
                        # Last current line is start of burst n+1
                        # Add new line for end of this burst
                        # (for line_time to be inside geo_line_time range)
                        new_line = geo_line[-1] + lines_per_burst - 1
                        new_time = geo_line_time[-1] + (new_line-geo_line[-1])*azitimeint
                        geo_line = np.append(geo_line, new_line)
                        geo_line_time = np.append(geo_line_time, new_time)
                    else:
                        raise Exception('Unexpected geolocation grid !')
                    # Compute time for each line
                    burst = np.searchsorted(geo_line, line, side='right')-1
                    line_time = (line-geo_line[burst])*azitimeint + \
                                geo_line_time[burst]
                elif geo_line.size == nbursts: # Strange case
                    # Check if geo_line is as expected
                    geo_line_exp = np.arange(nbursts)*lines_per_burst
                    geo_line_exp[-1] -= 1
                    if (geo_line == geo_line_exp).all() == False:
                        raise Exception('Unexpected geolocation grid !')
                    # Last current line is end of burst n-1
                    # Modify it temporarily to the start of burst n
                    geo_line[-1] += 1
                    orig_geo_line_time = geo_line_time[-1]
                    dtime = geo_line_time[-3] - geo_line_time[-4]
                    geo_line_time[-1] = geo_line_time[-2] + dtime
                    # Add new line for end of burst n
                    # (for line_time to be inside geo_line_time range)
                    new_line = self.get_dimsize('row') - 1
                    new_time = geo_line_time[-1] + (new_line-geo_line[-1])*azitimeint
                    geo_line = np.append(geo_line, new_line)
                    geo_line_time = np.append(geo_line_time, new_time)
                    # Compute time for each line
                    burst = np.searchsorted(geo_line, line, side='right')-1
                    line_time = (line-geo_line[burst])*azitimeint + \
                                geo_line_time[burst]
                    # Rewrite end of burst n-1 because it is an actual GCP
                    geo_line[-2] -= 1
                    geo_line_time[-2] = orig_geo_line_time
                else:
                    raise Exception('Unexpected geolocation grid !')
                func = interpolate.interp1d(geo_line_time, geo_line,
                                            kind='linear')
                line = func(line_time)
            values = self._interpolate_lonlat(fieldname, slices, line, pixel)
        elif fieldname in ['incidence', 'elevation']:
            geoloc = self._get_geolocation_grid()
            xval, yval = geoloc['line'], geoloc['pixel']
            val = geoloc[fieldname+'_angle']
            # values = self._interpolate_values(xval, yval, val, slices,
            #                                  intmethod='griddata',
            #                                  method='cubic')
            # values = self._interpolate_values(xval, yval, val, slices,
            #                                  intmethod='interp2d',
            #                                  kind='cubic')
            values = self._interpolate_values(xval, yval, val, slices,
                                              intmethod='RectBivariateSpline',
                                              kx=3, ky=3)
        elif fieldname == 'digital_number':
            band = self.get_handler().GetRasterBand(1)
            values = band.ReadAsArray(int(slices[1].start),
                                      int(slices[0].start),
                                      int(slices[1].stop-slices[1].start),
                                      int(slices[0].stop-slices[0].start))
            values = values[::slices[0].step, ::slices[1].step]
            # attrs = self.read_field_attributes(fieldname)
            # if '_FillValue' in attrs:
            #     fill_value = attrs['_FillValue']
            # else:
            #     fill_value = None
            # if not fill_value is None:
            #     values = np.ma.array(values, fill_value=fill_value)
            # else:
            #     values = np.ma.array(values)
            # if 'scale_factor' in attrs:
            #     values = values * attrs['scale_factor']
            # if 'add_offset' in attrs:
            #     values = values + attrs['add_offset']
        elif fieldname in ['gamma0_lut', 'beta0_lut', 'sigma0_lut']:
            luts = self._get_lookup_tables()
            xval, yval = luts['line'], luts['pixel']
            val = luts[fieldname.replace('_lut', '')]
            # values = self._interpolate_values(xval, yval, val, slices,
            #                                  intmethod='griddata',
            #                                  method='linear')
            values = self._interpolate_values(xval, yval, val, slices,
                                              intmethod='RectBivariateSpline',
                                              kx=1, ky=1)
        elif fieldname == 'gamma0':
            dnum = self.read_values('digital_number', slices=slices)
            lut = self.read_values('gamma0_lut', slices=slices)
            values = (abs(dnum).astype('float32'))**2 / lut**2
        #elif fieldname in ['beta0', 'intensity']:
        elif fieldname == 'beta0':
            dnum = self.read_values('digital_number', slices=slices)
            lut = self.read_values('beta0_lut', slices=slices)
            values = (abs(dnum).astype('float32'))**2 / lut**2
        elif fieldname == 'sigma0':
            dnum = self.read_values('digital_number', slices=slices)
            lut = self.read_values('sigma0_lut', slices=slices)
            values = (abs(dnum).astype('float32'))**2 / lut**2
        # elif fieldname == 'amplitude':
        #     dnum = self.read_values('digital_number', slices=slices)
        #     lut = self.read_values('beta0_lut', slices=slices)
        #     values = abs(dnum).astype('float32') / lut
        elif fieldname == 'complex':
            dnum = self.read_values('digital_number', slices=slices)
            lut = self.read_values('sigma0_lut', slices=slices)
            values = dnum / lut
        elif fieldname == 'azimuth_time':
            geoloc = self._get_geolocation_grid()
            xval, yval = geoloc['line'], geoloc['pixel']
            val = geoloc['azimuth_time']
            values = self._interpolate_values(xval, yval, val, slices,
                                              intmethod='RectBivariateSpline',
                                              kx=1, ky=1)
        elif fieldname == 'slant_range_time':
            geoloc = self._get_geolocation_grid()
            xval, yval = geoloc['line'], geoloc['pixel']
            val = geoloc['slant_range_time']
            values = self._interpolate_values(xval, yval, val, slices,
                                              intmethod='RectBivariateSpline',
                                              kx=1, ky=3)
        elif fieldname == 'doppler_centroid':
            dce = self._get_doppler_centroid_estimates()
            xval = dce['azimuth_time'][:, 0]
            yval = dce['slant_range_time'][0, :]
            val = dce['frequency']
            if xval.size == 1:
                xval = np.append(xval, dce['azimuth_time_stop'][0])
                val = np.repeat(val, 2, axis=0)
            if yval.size == 1:
                msg = 'only one slant range for doppler centroid estimates ?'
                raise BaseException(msg)
            kx, ky = min(xval.size-1, 3), min(yval.size-1, 3)
            func = interpolate.RectBivariateSpline(xval, yval, val, kx=kx,
                                                   ky=ky)
            # xint = self.read_values('azimuth_time', slices=slices)[:, 0]
            # yint = self.read_values('slant_range_time', slices=slices)[0, :]
            mid = (self._attributes['number_of_lines']/2,
                   self._attributes['number_of_samples']/2)
            azislices = list(slices)
            azislices[1] = slice(mid[1], mid[1]+1, 1)
            xint = self.read_values('azimuth_time', slices=azislices)
            ranslices = list(slices)
            ranslices[0] = slice(mid[0], mid[0]+1, 1)
            yint = self.read_values('slant_range_time', slices=ranslices)
            values = func(xint, yint).astype('float32')
        else:
            raise NotImplementedError
        return values

    def read_global_attributes(self):
        """
        """
        return self._attributes

    def read_global_attribute(self, attr):
        """
        """
        return self._attributes[attr]

    def write_field(self, fieldname):
        """
        """
        raise NotImplementedError

    def read_fillvalue(self, fieldname):
        """
        """
        return None

    def create_field(self, field, dim_translation=None):
        """
        """
        raise NotImplementedError

    def create_dim(self, dimname, size=None):
        """
        """
        raise NotImplementedError

    def write_global_attributes(self, attrs):
        """
        Write the storage (file) global attributes.
        """
        raise NotImplementedError

    def get_start_time(self):
        """
        Return the minimum date of the file temporal coverage.
        """
        return datetime.strptime(self._attributes['start_time'],
                                 '%Y-%m-%dT%H:%M:%S.%f')

    def get_end_time(self):
        """
        Return the maximum date of the file temporal coverage.
        """
        return datetime.strptime(self._attributes['stop_time'],
                                 '%Y-%m-%dT%H:%M:%S.%f')

    def get_spatial_resolution_in_deg(self):
        """
        Return the average spatial resolution in degrees.
        """
        raise NotImplementedError

    def get_bbox(self):
        """
        Return the bounding box of the feature, as a tuple
        (lonmin, latmin, lonmax, latmax).
        """
        geoloc = self._get_geolocation_grid()
        return (geoloc['longitude'].min(), geoloc['latitude'].min(),
                geoloc['longitude'].max(), geoloc['latitude'].max())

    def _strtime2numtime(self, strtime, fmt='%Y-%m-%dT%H:%M:%S.%f'):
        """
        Convert string time to numeric time.
        """
        dtime = datetime.strptime(strtime, fmt)
        numtime = date2num(dtime, self.read_field('time').units)
        return numtime

    def _numtime2strtime(self, numtime, fmt='%Y-%m-%dT%H:%M:%S.%f'):
        """
        Convert numeric time to string time.
        """
        dtime = num2date(numtime, self.read_field('time').units)
        strtime = dtime.strftime(fmt)
        return strtime

    def _get_geolocation_grid(self):
        """
        """
        return self._attributes['geolocation_grid']

    def _get_doppler_centroid_estimates(self):
        """
        """
        return self._attributes['doppler_centroid_estimates']

    def _get_lookup_tables(self):
        """
        """
        return self._attributes['lookup_tables']

    # CHANGE NAME
    def _interpolate_lonlat(self, fieldname, slices, line, pixel):
        """
        """
        if slices != self._lonlat_slices:
            # Create GDAL transformer
            src_ds, dst_ds = self.get_handler(), None
            #options = ['']
            #options = ['MAX_GCP_ORDER=3']
            options = ['MAX_GCP_ORDER=-1']
            transformer = gdal.Transformer(src_ds, dst_ds, options)
            # Apply transformer
            #xy_grid = np.mgrid[slices[1], slices[0]].astype('int32')
            xy_grid = np.array((np.tile(pixel.reshape((-1, 1)), (1, line.size)),
                                np.tile(line.reshape((1, -1)), (pixel.size, 1))),
                               dtype='int32')
            xy_dims = xy_grid.shape[1:3]
            xy_grid = xy_grid.reshape((2, -1)).transpose()
            lonlat = transformer.TransformPoints(0, xy_grid)
            lonlat = np.array(lonlat[0], dtype='float32')
            self._lon = lonlat[:, 0].reshape(xy_dims).transpose()
            self._lat = lonlat[:, 1].reshape(xy_dims).transpose()
            self._lonlat_slices = slices
        if fieldname == 'lat':
            values = self._lat
        elif fieldname == 'lon':
            values = self._lon
        return values

    def _interpolate_values(self, xval, yval, val, slices,
                            intmethod='RectBivariateSpline', **kwargs):
        """
        """
        valtype = val.dtype
        if intmethod == 'griddata':
            xint, yint = np.mgrid[slices[0], slices[1]].astype('int32')
            values = interpolate.griddata((xval.flatten(), yval.flatten()),
                                          val.flatten(), (xint, yint), **kwargs)
        elif intmethod == 'interp2d':  # not working, pb with knots
            func = interpolate.interp2d(xval.flatten(), yval.flatten(),
                                        val.flatten(), **kwargs)
            xint = np.arange(slices[0].start, slices[0].stop,
                             slices[0].step, dtype='int32')
            yint = np.arange(slices[1].start, slices[1].stop,
                             slices[1].step, dtype='int32')
            values = func(yint, xint)
        elif intmethod == 'RectBivariateSpline':
            func = interpolate.RectBivariateSpline(xval[:, 0], yval[0, :],
                                                   val, **kwargs)
            xint = np.arange(slices[0].start, slices[0].stop,
                             slices[0].step, dtype='int32')
            yint = np.arange(slices[1].start, slices[1].stop,
                             slices[1].step, dtype='int32')
            values = func(xint, yint)
        values = values.astype(valtype)
        return values

    # Would be better to read it in a file but it is not in annotation files.
    # It is in manifest file we don't want to parse only for that.
    # So let's suppose SAFE format has to stay SAFE format ...
    def _get_safe_name(self):
        """
        Return SAFE name.
        """
        safename = os.path.basename(os.path.dirname(os.path.dirname(self._url)))
        return safename

    def _parse_annotationdataset(self, name):
        """
        Parse annotation data set.

        :param name: name of the annotation data set
            (product, calibration, noise)
        :type name: str
        """
        if name == 'product':
            adsdir = os.path.join(os.path.dirname(self._url),
                                  '../annotation/')
            adsfile = os.path.basename(self._url).replace('.tiff', '.xml')
        elif name == 'calibration':
            adsdir = os.path.join(os.path.dirname(self._url),
                                  '../annotation/calibration/')
            adsfile = 'calibration-' + \
                      os.path.basename(self._url).replace('.tiff', '.xml')
        elif name == 'noise':
            adsdir = os.path.join(os.path.dirname(self._url),
                                  '../annotation/calibration/')
            adsfile = 'noise-' + \
                      os.path.basename(self._url).replace('.tiff', '.xml')
        else:
            raise Exception('Wrong ads name : '+name)
        tree = ET.parse(os.path.join(adsdir, adsfile))
        ads = tree.getroot()
        return ads

    def _read_annotationdatasets(self):
        """
        Read annotation data sets.
        """
        dic = OrderedDict()
        # Parse annotation data sets
        pads = self._parse_annotationdataset('product')
        cads = self._parse_annotationdataset('calibration')
        # nads = self._parse_annotationdataset('noise')
        # Fill from annotation data sets
        ### {product,calibration,noise}/adsHeader
        ### status : ALL
        self._fill_adsheader(pads, dic)
        ### product/qualityInformation
        ### status : NONE
        ### product/generalAnnotation/productInformation
        ### status : ALL
        self._fill_productinformation(pads, dic)
        ### product/generalAnnotation/downlinkInformationList
        ### status : PARTIAL
        ### warning : possibly more than one prf for GRD products
        self._fill_downlinkinformation(pads, dic)
        ### product/generalAnnotation/orbitList
        ### status : ALL
        self._fill_orbit(pads, dic)
        hplat = np.sqrt(np.sum(dic['orbit_state_position']**2))
        vplat = np.sqrt(np.sum(dic['orbit_state_velocity']**2))
        vgrnd = vplat/hplat*EARTHRADIUS
        dic['platform_height'] = hplat
        dic['platform_velocity'] = vplat
        dic['ground_velocity'] = vgrnd.astype('float32')
        ### product/generalAnnotation/attitudeList
        ### status : NONE
        ### product/generalAnnotation/rawDataAnalysisList
        ### status : NONE
        ### product/generalAnnotation/replicaInformationList
        ### status : NONE
        ### product/generalAnnotation/noiseList
        ### status : NONE
        ### product/generalAnnotation/terrainHeightList
        ### status : NONE
        ### product/generalAnnotation/azimuthFmRateList
        ### status : ALL
        self._fill_azimuthfmrate(pads, dic)
        ### product/imageAnnotation/imageInformation
        ### status : PARTIAL (sliceList and imageStatistics miss)
        self._fill_imageinformation(pads, dic)
        dic['azimuth_ground_spacing'] = dic['azimuth_pixel_spacing']
        dic['range_ground_spacing'] = dic['range_pixel_spacing']
        if dic['product'] == 'SLC':
            incmid = dic['incidence_angle_mid_swath']*np.pi/180
            dic['range_ground_spacing'] /= np.sin(incmid)
        ### product/imageAnnotation/processingInformation
        ### status : PARTIAL
        ### todo : window and bandwidth infos in ./swathProcParamsList
        self._fill_processinginformation(pads, dic)
        ### product/dopplerCentroid
        ### status : ALL
        self._fill_dopplercentroid(pads, dic)
        ### product/antennaPattern
        ### status : NONE
        ### product/swathTiming
        ### status : ALL
        self._fill_swathtiming(pads, dic)
        ### product/geolocationGrid
        ### status : ALL
        self._fill_geolocationgrid(pads, dic)
        ### product/coordinateConversion
        ### status : NONE
        ### product/swathMerging
        ### status : NONE
        ### calibration/calibrationInformation
        ### calibration/calibrationVectorList
        ### status : ALL
        self._fill_calibration(cads, dic)
        ### noise/noiseVectorList
        ### status : NONE
        return dic

    @staticmethod
    def _fill_adsheader(ads, dic):
        """
        """
        adsh = ads.find('./adsHeader')
        dic['mission'] = adsh.find('./missionId').text
        dic['product'] = adsh.find('./productType').text
        dic['polarisation'] = adsh.find('./polarisation').text
        dic['mode'] = adsh.find('./mode').text
        dic['swath'] = adsh.find('./swath').text
        dic['start_time'] = adsh.find('./startTime').text
        dic['stop_time'] = adsh.find('./stopTime').text
        dic['absolute_orbit'] = int(adsh.find('./absoluteOrbitNumber').text)
        dic['mission_datatake_id'] = adsh.find('./missionDataTakeId').text
        dic['image_number'] = int(adsh.find('./imageNumber').text)

    @staticmethod
    def _fill_productinformation(pads, dic):
        """
        """
        ads = pads.find('./generalAnnotation/productInformation')
        dic['pass'] = ads.find('./pass').text
        dic['timeliness_category'] = ads.find('./timelinessCategory').text
        dic['platform_heading'] = float(ads.find('./platformHeading').text)
        dic['projection'] = ads.find('./projection').text
        dic['range_sampling_rate'] = float(ads.find('./rangeSamplingRate').text)
        dic['radar_frequency'] = float(ads.find('./radarFrequency').text)
        dic['azimuth_steering_rate'] = \
            float(ads.find('azimuthSteeringRate').text)

    @staticmethod
    def _fill_downlinkinformation(pads, dic):
        """
        """
        ads = pads.find('./generalAnnotation/downlinkInformationList')
        dlinf = ads.findall('./downlinkInformation')
        ndlinf = len(dlinf)
        prf = np.zeros(ndlinf, dtype='float32')
        for iprf in np.arange(ndlinf):
            prf[iprf] = np.float32(dlinf[iprf].find('./prf').text)
        dic['pulse_repetition_frequencies'] = prf
        if ndlinf == 1:
            dic['pulse_repetition_frequency'] = prf[0]

    #@staticmethod
    def _fill_orbit(self, pads, dic):
        """
        """
        vectors = pads.findall('./generalAnnotation/orbitList/orbit')
        nvect = len(vectors)
        osv = OrderedDict()
        osv['nlines'] = nvect
        osv['time'] = np.empty(nvect, dtype='float64')
        osv['frame'] = []
        osv['position'] = np.empty((nvect, 3), dtype='float32')
        osv['velocity'] = np.empty((nvect, 3), dtype='float32')
        for ivect, vector in enumerate(vectors):
            strtime = vector.find('./time').text
            osv['time'][ivect] = self._strtime2numtime(strtime)
            osv['frame'].append(vector.find('./frame').text)
            osv['position'][ivect, 0] = vector.find('./position/x').text
            osv['position'][ivect, 1] = vector.find('./position/y').text
            osv['position'][ivect, 2] = vector.find('./position/z').text
            osv['velocity'][ivect, 0] = vector.find('./velocity/x').text
            osv['velocity'][ivect, 1] = vector.find('./velocity/y').text
            osv['velocity'][ivect, 2] = vector.find('./velocity/z').text
        # osv['total_velocity'] = np.sqrt(osv['velocity'][:, 0]**2 + \
        #                                 osv['velocity'][:, 1]**2 + \
        #                                 osv['velocity'][:, 2]**2)
        dic['orbit_state_vectors'] = osv
        dic['orbit_state_position'] = osv['position'][nvect//2, :]
        dic['orbit_state_velocity'] = osv['velocity'][nvect//2, :]

    #@staticmethod
    def _fill_azimuthfmrate(self, pads, dic):
        """
        """
        coeffs = pads.findall('./generalAnnotation/azimuthFmRateList/azimuthFmRate')
        ncoeff = len(coeffs)
        afr = OrderedDict()
        afr['nlines'] = ncoeff
        afr['azimuth_time'] = np.empty(ncoeff, dtype='float64')
        afr['t0'] = np.empty(ncoeff, dtype='float32')
        afr['c0'] = np.empty(ncoeff, dtype='float32')
        afr['c1'] = np.empty(ncoeff, dtype='float32')
        afr['c2'] = np.empty(ncoeff, dtype='float32')
        for icoeff, coeff in enumerate(coeffs):
            strtime = coeff.find('./azimuthTime').text
            afr['azimuth_time'][icoeff] = self._strtime2numtime(strtime)
            afr['t0'][icoeff] = coeff.find('./t0').text
            poly1 = [coeff.find('./'+cname) for cname in ['c0', 'c1', 'c2']]
            poly2 = coeff.find('./azimuthFmRatePolynomial')
            if all([p is not None for p in poly1]): # old annotation
                polycoeff = [p.text for p in poly1]
            elif poly2 is not None: # new annotation (if not bug)
                polycoeff = poly2.text.split(' ')
            else:
                raise Exception('Could not find azimuth FM rate polynomial coefficients')
            afr['c0'][icoeff] = polycoeff[0]
            afr['c1'][icoeff] = polycoeff[1]
            afr['c2'][icoeff] = polycoeff[2]
        dic['azimuthfmrate_list'] = afr

    @staticmethod
    def _fill_imageinformation(pads, dic):
        """
        """
        ads = pads.find('./imageAnnotation/imageInformation')
        dic['first_line_time'] = ads.find('./productFirstLineUtcTime').text
        dic['last_line_time'] = ads.find('./productLastLineUtcTime').text
        dic['ascending_node_time'] = ads.find('./ascendingNodeTime').text
        dic['anchor_time'] = ads.find('./anchorTime').text
        dic['product_composition'] = ads.find('./productComposition').text
        dic['slice_number'] = int(ads.find('./sliceNumber').text)
        dic['slant_range_time'] = float(ads.find('./slantRangeTime').text)
        dic['pixel_value'] = ads.find('./pixelValue').text
        dic['output_pixels'] = ads.find('./outputPixels').text
        dic['range_pixel_spacing'] = float(ads.find('./rangePixelSpacing').text)
        dic['azimuth_pixel_spacing'] = \
            float(ads.find('./azimuthPixelSpacing').text)
        dic['azimuth_time_interval'] = \
            float(ads.find('./azimuthTimeInterval').text)
        dic['azimuth_frequency'] = float(ads.find('./azimuthFrequency').text)
        dic['number_of_samples'] = int(ads.find('./numberOfSamples').text)
        dic['number_of_lines'] = int(ads.find('./numberOfLines').text)
        dic['zero_dop_minus_acq_time'] = \
            float(ads.find('./zeroDopMinusAcqTime').text)
        dic['incidence_angle_mid_swath'] = \
            float(ads.find('./incidenceAngleMidSwath').text)

    @staticmethod
    def _fill_processinginformation(pads, dic):
        """
        """
        ads = pads.find('./imageAnnotation/processingInformation')
        dic['reference_range'] = float(ads.find('./referenceRange').text)

    #@staticmethod
    def _fill_dopplercentroid(self, pads, dic):
        """
        """
        dce = OrderedDict()
        estimates = pads.findall('./dopplerCentroid/dcEstimateList/dcEstimate')
        dce['nlines'] = len(estimates)
        dce['npixels'] = int(estimates[0].find('fineDceList').get('count'))
        dce['ngeocoeffs'] = \
            int(estimates[0].find('geometryDcPolynomial').get('count'))
        dce['ndatacoeffs'] = \
            int(estimates[0].find('dataDcPolynomial').get('count'))
        dims = (dce['nlines'], dce['npixels'])
        dce['azimuth_time'] = np.empty(dims, dtype='float64')
        dce['t0'] = np.empty(dce['nlines'], dtype='float64')
        dce['geo_polynom'] = np.empty((dce['nlines'], dce['ngeocoeffs']),
                                      dtype='float32')
        dce['data_polynom'] = np.empty((dce['nlines'], dce['ndatacoeffs']),
                                       dtype='float32')
        dce['data_rms'] = np.empty(dce['nlines'], dtype='float32')
        #dce['data_rms_threshold'] =
        dce['azimuth_time_start'] = np.empty(dce['nlines'], dtype='float64')
        dce['azimuth_time_stop'] = np.empty(dce['nlines'], dtype='float64')
        dce['slant_range_time'] = np.empty(dims, dtype='float64')
        dce['frequency'] = np.empty(dims, dtype='float32')
        for iline, estimate in enumerate(estimates):
            strtime = estimate.find('./azimuthTime').text
            dce['azimuth_time'][iline, :] = self._strtime2numtime(strtime)
            dce['t0'][iline] = estimate.find('./t0').text
            dce['geo_polynom'][iline, :] = \
                estimate.find('./geometryDcPolynomial').text.split()
            dce['data_polynom'][iline, :] = \
                estimate.find('./dataDcPolynomial').text.split()
            dce['data_rms'][iline] = estimate.find('./dataDcRmsError').text
            #dce['data_rms_threshold'] =
            strtime = estimate.find('./fineDceAzimuthStartTime').text
            dce['azimuth_time_start'][iline] = self._strtime2numtime(strtime)
            strtime = estimate.find('./fineDceAzimuthStopTime').text
            dce['azimuth_time_stop'][iline] = self._strtime2numtime(strtime)
            finedces = estimate.findall('./fineDceList/fineDce')
            for ipixel, finedce in enumerate(finedces):
                dce['slant_range_time'][iline, ipixel] = \
                    finedce.find('./slantRangeTime').text
                dce['frequency'][iline, ipixel] = \
                    finedce.find('./frequency').text
        dic['doppler_centroid_estimates'] = dce

    #@staticmethod
    def _fill_swathtiming(self, pads, dic):
        """
        """
        ads = pads.find('./swathTiming')
        dic['lines_per_burst'] = int(ads.find('./linesPerBurst').text)
        dic['samples_per_burst'] = int(ads.find('./samplesPerBurst').text)
        ads = ads.find('./burstList')
        dic['number_of_bursts'] = int(ads.get('count'))
        burstlist = OrderedDict()
        if dic['number_of_bursts'] != 0:
            bursts = ads.findall('./burst')
            nbursts = len(bursts)
            burstlist['nbursts'] = nbursts
            burstlist['azimuth_time'] = np.empty(nbursts, dtype='float64')
            burstlist['azimuth_anx_time'] = np.empty(nbursts, dtype='float64')
            burstlist['sensing_time'] = np.empty(nbursts, dtype='float64')
            burstlist['byte_offset'] = np.empty(nbursts, dtype='uint64')
            nlines = dic['lines_per_burst']
            shp = (nbursts, nlines)
            burstlist['first_valid_sample'] = np.empty(shp, dtype='int32')
            burstlist['last_valid_sample'] = np.empty(shp, dtype='int32')
            burstlist['valid_location'] = np.empty((nbursts, 4), dtype='int32')
            for ibur, burst in enumerate(bursts):
                strtime = burst.find('./azimuthTime').text
                burstlist['azimuth_time'][ibur] = self._strtime2numtime(strtime)
                burstlist['azimuth_anx_time'][ibur] = \
                    np.float64(burst.find('./azimuthAnxTime').text)
                strtime = burst.find('./sensingTime').text
                burstlist['sensing_time'][ibur] = self._strtime2numtime(strtime)
                burstlist['byte_offset'][ibur] = \
                    np.uint64(burst.find('./byteOffset').text)
                fvs = np.int32(burst.find('./firstValidSample').text.split())
                burstlist['first_valid_sample'][ibur, :] = fvs
                lvs = np.int32(burst.find('./lastValidSample').text.split())
                burstlist['last_valid_sample'][ibur, :] = lvs
                valind = np.where((fvs != -1) | (lvs != -1))[0]
                valloc = [ibur*nlines+valind.min(), fvs[valind].min(),
                          ibur*nlines+valind.max(), lvs[valind].max()]
                burstlist['valid_location'][ibur, :] = valloc
        dic['burst_list'] = burstlist

    #@staticmethod
    def _fill_geolocationgrid(self, pads, dic):
        """
        """
        ads = pads.find('./geolocationGrid/geolocationGridPointList')
        points = ads.findall('./geolocationGridPoint')
        npoints = len(points)
        geoloc = OrderedDict()
        geoloc['npoints'] = npoints
        geoloc['nlines'] = 0
        geoloc['npixels'] = 0
        geoloc['azimuth_time'] = np.empty(npoints, dtype='float64')
        geoloc['slant_range_time'] = np.empty(npoints, dtype='float64')
        geoloc['line'] = np.empty(npoints, dtype='int32')
        geoloc['pixel'] = np.empty(npoints, dtype='int32')
        geoloc['latitude'] = np.empty(npoints, dtype='float32')
        geoloc['longitude'] = np.empty(npoints, dtype='float32')
        geoloc['height'] = np.empty(npoints, dtype='float32')
        geoloc['incidence_angle'] = np.empty(npoints, dtype='float32')
        geoloc['elevation_angle'] = np.empty(npoints, dtype='float32')
        gridkeys = [k for k in geoloc if isinstance(geoloc[k], np.ndarray)]
        for ipt, point in enumerate(points):
            strtime = point.find('./azimuthTime').text
            geoloc['azimuth_time'][ipt] = self._strtime2numtime(strtime)
            geoloc['slant_range_time'][ipt] = \
                point.find('./slantRangeTime').text
            geoloc['line'][ipt] = point.find('./line').text
            geoloc['pixel'][ipt] = point.find('./pixel').text
            geoloc['latitude'][ipt] = point.find('./latitude').text
            geoloc['longitude'][ipt] = point.find('./longitude').text
            geoloc['height'][ipt] = point.find('./height').text
            geoloc['incidence_angle'][ipt] = \
                point.find('./incidenceAngle').text
            geoloc['elevation_angle'][ipt] = \
                point.find('./elevationAngle').text
        # Make grid
        lines = geoloc['line']
        geoloc['npixels'] = len(np.where(lines == lines[0])[0])
        pixels = geoloc['pixel']
        geoloc['nlines'] = len(np.where(pixels == pixels[0])[0])
        if geoloc['nlines']*geoloc['npixels'] != geoloc['npoints']:
            msg = 'nlines*npixels != npoints in geolocation grid'
            raise BaseException(msg)
        newshape = (geoloc['nlines'], geoloc['npixels'])
        for key in gridkeys:
            geoloc[key] = np.reshape(geoloc[key], newshape)
        # Check grid
        chk = geoloc['line'].min(axis=1) != geoloc['line'].max(axis=1)
        if chk.any():
            raise BaseException('geolocation lines do not form a grid')
        chk = geoloc['pixel'].min(axis=0) != geoloc['pixel'].max(axis=0)
        if chk.any():
            raise BaseException('geolocation pixels do not form a grid')
        linargsort = geoloc['line'][:, 0].argsort()
        if ((linargsort[1:] - linargsort[0:-1]) != 1).any():
            warnings.warn('Geolocation lines are not sorted -> do sort')
            for key in gridkeys:
                geoloc[key] = geoloc[key][linargsort, :]
        pixargsort = geoloc['pixel'][0, :].argsort()
        if ((pixargsort[1:] - pixargsort[0:-1]) != 1).any():
            warnings.warn('Geolocation pixels are not sorted -> do sort')
            for key in gridkeys:
                geoloc[key] = geoloc[key][:, pixargsort]
        dic['geolocation_grid'] = geoloc

    #@staticmethod
    def _fill_calibration(self, cads, dic):
        """
        """
        ads = cads.find('calibrationInformation')
        dic['absolute_calibration_constant'] = \
            float(ads.find('absoluteCalibrationConstant').text)
        luts = OrderedDict()
        vectors = cads.findall('./calibrationVectorList/calibrationVector')
        luts['nlines'] = len(vectors)
        luts['npixels'] = int(vectors[0].find('./pixel').get('count'))
        dims = (luts['nlines'], luts['npixels'])
        luts['azimuth_time'] = np.empty(dims, dtype='float64')
        luts['line'] = np.empty(dims, dtype='int32')
        luts['pixel'] = np.empty(dims, dtype='int32')
        luts['sigma0'] = np.empty(dims, dtype='float32')
        luts['beta0'] = np.empty(dims, dtype='float32')
        luts['gamma0'] = np.empty(dims, dtype='float32')
        luts['digital_number'] = np.empty(dims, dtype='float32')
        for iline, vector in enumerate(vectors):
            strtime = vector.find('./azimuthTime').text
            luts['azimuth_time'][iline, :] = self._strtime2numtime(strtime)
            luts['line'][iline, :] = vector.find('./line').text
            luts['pixel'][iline, :] = vector.find('./pixel').text.split()
            luts['sigma0'][iline, :] = \
                vector.find('./sigmaNought').text.split()
            luts['beta0'][iline, :] = \
                vector.find('./betaNought').text.split()
            luts['gamma0'][iline, :] = vector.find('./gamma').text.split()
            luts['digital_number'][iline, :] = \
                vector.find('./dn').text.split()
        dic['lookup_tables'] = luts
