#!/usr/bin/env python
# coding=utf-8
"""
"""


import os
from collections import OrderedDict
import xml.etree.ElementTree as ET
import numpy as np
import pyproj
from netCDF4 import date2num
from osgeo import gdal, ogr
from datetime import datetime
import tempfile
import shutil

from cerbere.mapper.abstractmapper import AbstractMapper
from cerbere import READ_ONLY
from cerbere.mapper.slices import Slices
from cerbere.datamodel.field import Field
from cerbere.datamodel.variable import Variable


class SAFEMSIL1CGranuleFile(AbstractMapper):
    """Abstract class for SAFE MSI L1C granule files.

    Args:
        url (str): the granule directory (i.e. a directory in the GRANULES
            directory of the SAFE product).

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

    NATRES2BANDS = {'10m': ['B02', 'B03', 'B04', 'B08'],
                    '20m': ['B05', 'B06', 'B07', 'B8A', 'B11', 'B12'],
                    '60m': ['B01', 'B09', 'B10']}
    BANDS = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07',
             'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']

    def __init__(self, url, native_resolution, overview_index=None,
                 tight=False, mode=READ_ONLY, **kwargs):
        """
        """
        if 'MSI_L1C_TL' not in os.path.basename(url):
            raise Exception('MSI L1C granule expected.')
        if native_resolution not in ['10m', '20m', '60m']:
            raise Exception('Wrong MSI L1C native resolution.')
        if mode != READ_ONLY:
            raise Exception('This mapper can only be used in read_only mode.')
        super(SAFEMSIL1CGranuleFile, self).__init__(url=url, mode=mode, **kwargs)
        self._native_resolution = native_resolution
        self._bands = self.NATRES2BANDS[native_resolution]
        self._overview_index = overview_index
        self._tight = tight
        self._tight_slices = None
        self._attrs = OrderedDict()
        self._dims = OrderedDict()
        self._fields = OrderedDict()
        self._band_handlers = {}
        self._mask_handlers = {}
        self._angle_grids = {}
        self._angle_handlers = {}
        self._initial_ul = None

    def _get_granule_path(self):
        """
        """
        return self.get_url()

    def _get_granule_metadata_path(self):
        """
        """
        granule_path = self._get_granule_path()
        granule_id = os.path.basename(granule_path)
        granmtd_id = '{}.xml'.format(granule_id[:-7].replace('MSI', 'MTD'))
        return os.path.join(granule_path, granmtd_id)

    def _get_safe_path(self):
        """
        """
        return os.path.dirname(os.path.dirname(self._get_granule_path()))

    def _get_safe_metadata_path(self):
        """
        """
        safe_path = self._get_safe_path()
        safe_id = os.path.basename(safe_path)
        safemtd_id = safe_id.replace('PRD_MSI', 'MTD_SAF').replace('.SAFE', '.xml')
        return os.path.join(safe_path, safemtd_id)

    def _get_band_path(self, band):
        """
        """
        granule_path = self._get_granule_path()
        granule_id = os.path.basename(granule_path)
        band_id = '{}_{}.jp2'.format(granule_id[:-7], band)
        return os.path.join(granule_path, 'IMG_DATA', band_id)

    def _get_mask_path(self, mask_type, band):
        """
        """
        granule_path = self._get_granule_path()
        granule_id = os.path.basename(granule_path)
        mask_prefix = granule_id[:-7].replace('MSI_L1C_TL',
                                              'MSK_{}'.format(mask_type.upper()))
        mask_id = '{}_{}_MSIL1C.gml'.format(mask_prefix, band)
        return os.path.join(granule_path, 'QI_DATA', mask_id)

    def _get_mask_handler(self, mask_type, band):
        """
        """
        mask_name = '{}_{}'.format(mask_type, band)
        if mask_name not in self._mask_handlers:
            mask_path = self._get_mask_path(mask_type, band)
            # Create a temporary directory and copy GML file.
            # (did not find how to prevent OGR from creating a .gfs file when
            # reading a GML file. The temporary directory is here to avoid to
            # mess original data up)
            dtemp = tempfile.mkdtemp()
            try:
                tmp_mask_path = os.path.join(dtemp, os.path.basename(mask_path))
                shutil.copyfile(mask_path, tmp_mask_path)
                # Read mask
                orig_ds = ogr.Open(tmp_mask_path)
                orig_layer = orig_ds.GetLayer(0)
                # Copy mask to memory
                drv = ogr.GetDriverByName('MEMORY')
                mem_ds = drv.CreateDataSource(mask_name)
                if orig_layer is not None:
                    mem_layer = mem_ds.CopyLayer(orig_layer, orig_layer.GetName())
                    # Modify mask
                    if mask_type == 'DETFOO':
                        fielddefn = ogr.FieldDefn('detector_index', ogr.OFTInteger)
                        mem_layer.CreateField(fielddefn)
                        for feature in mem_layer:
                            field_index = feature.GetFieldIndex('gml_id')
                            gml_id = feature.GetField(field_index)
                            detector_index = int(gml_id.split('-')[2])
                            feature.SetField('detector_index', detector_index)
                            mem_layer.SetFeature(feature)
                orig_ds = None
                self._mask_handlers[mask_name] = mem_ds
            except:
                raise
            finally:
                shutil.rmtree(dtemp)
        return self._mask_handlers[mask_name]

    def _get_band_handler(self, band):
        """
        """
        if band not in self._band_handlers:
            band_path = self._get_band_path(band)
            handler = gdal.Open(band_path)
            self._band_handlers[band] = handler
        return self._band_handlers[band]

    def _get_angle_handler(self, angle_type, band=None, detector=None,
                           nodatavalue=None):
        """
        """
        angle_name = angle_type
        if band is not None and detector is not None:
            angle_name = '{}_{:02d}_{}'.format(band, detector, angle_name)
        if angle_name not in self._angle_handlers:
            grid = self._angle_grids[angle_name]
            xdim = grid['col_step']
            if np.sign(xdim) != np.sign(self.read_global_attribute('xdim')):
                xdim = -xdim
            ydim = grid['row_step']
            if np.sign(ydim) != np.sign(self.read_global_attribute('ydim')):
                ydim = -ydim
            values = grid['values']
            drv = gdal.GetDriverByName('MEM')
            dset = drv.Create(angle_name, values.shape[1], values.shape[0], 1,
                              gdal.GDT_Float32)
            ulx = self._initial_ul[0]
            uly = self._initial_ul[1]
            geotf = [ulx - xdim / 2., xdim, 0,
                     uly - ydim / 2., 0, ydim]
            dset.SetGeoTransform(geotf)
            if nodatavalue is not None:
                dset.GetRasterBand(1).SetNoDataValue(nodatavalue)
                values[np.where(~np.isfinite(values))] = nodatavalue
            dset.GetRasterBand(1).WriteArray(values)
            self._angle_handlers[angle_name] = dset
        return self._angle_handlers[angle_name]

    def _get_xy_envelope(self):
        """
        """
        xmin, xmax, ymin, ymax = [], [], [], []
        xmlns = {'gml': 'http://www.opengis.net/gml/3.2'}
        for band in self._bands:
            detfoo_path = self._get_mask_path('DETFOO', band)
            tree = ET.parse(detfoo_path)
            root = tree.getroot()
            env = root.find('gml:boundedBy/gml:Envelope', xmlns)
            lowcorn = env.find('gml:lowerCorner', xmlns).text.split(' ')
            xmin.append(float(lowcorn[0]))
            ymin.append(float(lowcorn[1]))
            upcorn = env.find('gml:upperCorner', xmlns).text.split(' ')
            xmax.append(float(upcorn[0]))
            ymax.append(float(upcorn[1]))
        return [min(xmin), max(xmax), min(ymin), max(ymax)]

    def _parse_angle_entry(self, entry):
        """
        """
        col_step = float(entry.find('COL_STEP').text)
        row_step = float(entry.find('ROW_STEP').text)
        values = entry.findall('Values_List/VALUES')
        shape = (len(values), len(values[0].text.split(' ')))
        grid = np.zeros(shape, dtype='float32')
        for line, val in enumerate(values):
            grid[line, :] = np.array([float(e) for e in val.text.split(' ')])
        return {'col_step': col_step, 'row_step': row_step, 'values': grid}

    def _get_gdal_dataset_from_slices(self, slices, btype=gdal.GDT_Byte,
                                      nodatavalue=None, fill=None,
                                      use_step=False):
        """
        """
        xdim = self.read_global_attribute('xdim')
        ulx = self.read_global_attribute('ulx')
        ydim = self.read_global_attribute('ydim')
        uly = self.read_global_attribute('uly')
        drv = gdal.GetDriverByName('MEM')
        if use_step == False:
            offext = slices.indices_offsetextent()
            dset = drv.Create('', offext[1][1], offext[0][1], 1, btype)
            geotf = [ulx + xdim * offext[1][0], xdim, 0,
                     uly + ydim * offext[0][0], 0, ydim]
        else:
            shape = slices.shape()
            dset = drv.Create('', shape[1], shape[0], 1, btype)
            minmax = slices.indices_minmax()
            if slices[1].step > 0:
                ulcx = ulx + xdim * minmax[1][0] + xdim / 2.
            else:
                ulcx = ulx + xdim * minmax[1][1] + xdim / 2.
            if slices[0].step > 0:
                ulcy = uly + ydim * minmax[0][0] + ydim / 2.
            else:
                ulcy = uly + ydim * minmax[0][1] + ydim / 2.
            newxdim = xdim * slices[1].step
            newydim = ydim * slices[0].step
            newulx = ulcx - newxdim / 2.
            newuly = ulcy - newydim / 2.
            geotf = [newulx, newxdim, 0, newuly, 0, newydim]
        dset.SetGeoTransform(geotf)
        if nodatavalue is not None:
            dset.GetRasterBand(1).SetNoDataValue(nodatavalue)
        if fill is not None:
            dset.GetRasterBand(1).Fill(fill)
        return dset

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
        # Define arbitrarily the granule metadata XML file as the mapper handler
        granmtd_path = self._get_granule_metadata_path()
        handler = open(granmtd_path, self._mode)
        self._handler = handler
        if datamodel is not None:
            self._feature_type = datamodel
        if datamodel_geolocation_dims is not None:
            self.datamodel_geolocation_dims = datamodel_geolocation_dims
        # Parse SAFE metadata XML
        safemtd_path = self._get_safe_metadata_path()
        tree = ET.parse(safemtd_path)
        root = tree.getroot()
        xmlns = {'n1': 'https://psd-13.sentinel2.eo.esa.int/PSD/User_Product_Level-1C.xsd',
                 'xsi': 'http://www.w3.org/2001/XMLSchema-instance'}
        prod = root.find('n1:General_Info/Product_Info', xmlns)
        query = ['PRODUCT_START_TIME', 'PRODUCT_STOP_TIME', 'PROCESSING_LEVEL',
                 'PRODUCT_TYPE', 'PROCESSING_BASELINE', 'GENERATION_TIME']
        for q in query:
            elem = prod.find(q)
            self._attrs[elem.tag.lower()] = elem.text
        dttk = prod.find('Datatake')
        query = ['SPACECRAFT_NAME', 'DATATAKE_TYPE', 'DATATAKE_SENSING_START',
                 'SENSING_ORBIT_NUMBER', 'SENSING_ORBIT_DIRECTION']
        for q in query:
            elem = dttk.find(q)
            if q == 'SENSING_ORBIT_NUMBER':
                self._attrs[elem.tag.lower()] = int(elem.text)
            else:
                self._attrs[elem.tag.lower()] = elem.text
        image = root.find('n1:General_Info/Product_Image_Characteristics', xmlns)
        for elem in image.findall('Special_Values'):
            if elem.find('SPECIAL_VALUE_TEXT').text == 'NODATA':
                value = int(elem.find('SPECIAL_VALUE_INDEX').text)
                self._attrs['nodata_value'] = value
            elif elem.find('SPECIAL_VALUE_TEXT').text == 'SATURATED':
                value = int(elem.find('SPECIAL_VALUE_INDEX').text)
                self._attrs['saturated_value'] = value
        query = ['QUANTIFICATION_VALUE']
        for q in query:
            elem = image.find(q)
            if q == 'QUANTIFICATION_VALUE':
                self._attrs[elem.tag.lower()] = float(elem.text)
            else:
                self._attrs[elem.tag.lower()] = elem.text
        del tree
        # Parse granule metadata XML
        tree = ET.parse(handler)
        root = tree.getroot()
        xmlns = {'n1': 'https://psd-12.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd',
                 'xsi': 'http://www.w3.org/2001/XMLSchema-instance'}
        gnrl = root.find('n1:General_Info', xmlns)
        query = ['TILE_ID', 'SENSING_TIME', 'Archiving_Info/ARCHIVING_CENTRE',
                 'Archiving_Info/ARCHIVING_TIME']
        for q in query:
            elem = gnrl.find(q)
            self._attrs[elem.tag.lower()] = elem.text
        geom = root.find('n1:Geometric_Info', xmlns)
        query = ['Tile_Geocoding/HORIZONTAL_CS_NAME',
                 'Tile_Geocoding/HORIZONTAL_CS_CODE']
        for q in query:
            elem = geom.find(q)
            self._attrs[elem.tag.lower()] = elem.text
        for elem in geom.findall('Tile_Geocoding/Size'):
            if elem.get('resolution') != self._native_resolution[0:2]:
                continue
            self._attrs['nrows'] = int(elem.find('NROWS').text)
            self._attrs['ncols'] = int(elem.find('NCOLS').text)
        for elem in geom.findall('Tile_Geocoding/Geoposition'):
            if elem.get('resolution') != self._native_resolution[0:2]:
                continue
            self._attrs['ulx'] = float(elem.find('ULX').text)
            self._attrs['uly'] = float(elem.find('ULY').text)
            self._attrs['xdim'] = float(elem.find('XDIM').text)
            self._attrs['ydim'] = float(elem.find('YDIM').text)
        sun = geom.find('Tile_Angles/Sun_Angles_Grid')
        for anglename in ['Zenith', 'Azimuth']:
            grid_name = 'sun_{}'.format(anglename.lower())
            entry = sun.find(anglename)
            self._angle_grids[grid_name] = self._parse_angle_entry(entry)
        viewings = geom.findall('Tile_Angles/Viewing_Incidence_Angles_Grids')
        for viewing in viewings:
            bandid = viewing.attrib['bandId']
            band = self.BANDS[int(bandid)]
            if band not in self._bands:
                continue
            detector = int(viewing.attrib['detectorId'])
            for anglename in ['Zenith', 'Azimuth']:
                grid_name = '{}_{:02d}_viewing_{}'.format(band, detector,
                                                          anglename.lower())
                entry = viewing.find(anglename)
                pentry = self._parse_angle_entry(entry)
                if np.isfinite(pentry['values']).any():
                    self._angle_grids[grid_name] = pentry
        del tree
        self._initial_ul = [self._attrs['ulx'], self._attrs['uly']]
        # Modify coordinates if an overview is requested
        if self._overview_index is not None:
            dset = self._get_band_handler(self._bands[0])
            bnd = dset.GetRasterBand(1)
            if self._overview_index >= bnd.GetOverviewCount():
                raise Exception('Wrong overview index.')
            ovbnd = bnd.GetOverview(self._overview_index)
            nrows = ovbnd.XSize
            ncols = ovbnd.YSize
            xfac = int(np.floor(float(self._attrs['nrows']) / nrows))
            yfac = int(np.floor(float(self._attrs['ncols']) / ncols))
            self._attrs['nrows'] = nrows
            self._attrs['ncols'] = ncols
            self._attrs['xdim'] *= xfac
            self._attrs['ydim'] *= yfac
        # Refine coordinates and dimensions according to real footprint
        if self._tight == True:
            env = self._get_xy_envelope()
            xminind = (env[0] - self._attrs['ulx']) / self._attrs['xdim']
            xmaxind = (env[1] - self._attrs['ulx']) / self._attrs['xdim']
            if self._attrs['xdim'] > 0:
                xindmin = int(np.maximum(0, np.floor(xminind)))
                xindmax = int(np.minimum(self._attrs['nrows'], np.ceil(xmaxind)))
            else:
                xindmin = int(np.maximum(0, np.floor(xmaxind)))
                xindmax = int(np.minimum(self._attrs['nrows'], np.ceil(xminind)))
            yminind = (env[2] - self._attrs['uly']) / self._attrs['ydim']
            ymaxind = (env[3] - self._attrs['uly']) / self._attrs['ydim']
            if self._attrs['ydim'] > 0:
                yindmin = int(np.maximum(0, np.floor(yminind)))
                yindmax = int(np.minimum(self._attrs['ncols'], np.ceil(ymaxind)))
            else:
                yindmin = int(np.maximum(0, np.floor(ymaxind)))
                yindmax = int(np.minimum(self._attrs['ncols'], np.ceil(yminind)))
            # Create tight slices (equivalent to view slices in other mappers)
            slices = [slice(yindmin, yindmax, 1), slice(xindmin, xindmax, 1)]
            dimsizes = [self._attrs['ncols'], self._attrs['nrows']]
            self._tight_slices = Slices(slices, dimsizes)
            # Update attributes
            self._attrs['ulx'] += xindmin * self._attrs['xdim']
            self._attrs['nrows'] = xindmax - xindmin
            self._attrs['uly'] += yindmin * self._attrs['ydim']
            self._attrs['ncols'] = yindmax - yindmin
        # Construct dimensions and fields
        self._dims['time'] = 1
        self._dims['y'] = self._attrs['ncols']
        self._dims['x'] = self._attrs['nrows']
        self._fields['time'] = Field(
            Variable(shortname='time',
                     description='sensing time',
                     standardname='time'),
            OrderedDict([('time', self._dims['time'])]),
            datatype=np.dtype(np.int64),
            units='microseconds since 2000-01-01T00:00:00Z')
        self._fields['lon'] = Field(
            Variable(shortname='lon',
                     description='longitude',
                     standardname='longitude'),
            OrderedDict([('y', self._dims['y']), ('x', self._dims['x'])]),
            datatype=np.dtype(np.float32),
            units='degrees_east')
        self._fields['lat'] = Field(
            Variable(shortname='lat',
                     description='latitude',
                     standardname='latitude'),
            OrderedDict([('y', self._dims['y']), ('x', self._dims['x'])]),
            datatype=np.dtype(np.float32),
            units='degrees_north')
        self._fields['x'] = Field(
            Variable(shortname='x',
                     description='x coordinate of projection',
                     standardname='projection_x_coordinate'),
            OrderedDict([('x', self._dims['x'])]),
            datatype=np.dtype(np.float32),
            units='m')
        self._fields['y'] = Field(
            Variable(shortname='y',
                     description='y coordinate of projection',
                     standardname='projection_y_coordinate'),
            OrderedDict([('y', self._dims['y'])]),
            datatype=np.dtype(np.float32),
            units='m')
        for bnd in self._bands:
            name = '{}_digital_number'.format(bnd)
            self._fields[name] = Field(
                Variable(shortname=name),
                OrderedDict([('y', self._dims['y']), ('x', self._dims['x'])]),
                datatype=np.dtype(np.uint16),
                fillvalue=np.uint16(0))
            if os.path.exists(self._get_mask_path('DETFOO', bnd)):
                name = '{}_detector_index'.format(bnd)
                self._fields[name] = Field(
                    Variable(shortname=name),
                    OrderedDict([('y', self._dims['y']), ('x', self._dims['x'])]),
                    datatype=np.dtype(np.uint8),
                    fillvalue=np.uint8(0))
            is_zenith = [ang.startswith(bnd) and ang.endswith('zenith') \
                         for ang in self._angle_grids]
            if any(is_zenith):
                name = '{}_viewing_zenith_angle'.format(bnd)
                self._fields[name] = Field(
                    Variable(shortname=name),
                    OrderedDict([('y', self._dims['y']), ('x', self._dims['x'])]),
                    datatype=np.dtype(np.float32),
                    fillvalue=np.float32(1000),
                    units='degrees')
            is_azimuth = [ang.startswith(bnd) and ang.endswith('azimuth') \
                          for ang in self._angle_grids]
            if any(is_azimuth):
                name = '{}_viewing_azimuth_angle'.format(bnd)
                self._fields[name] = Field(
                    Variable(shortname=name),
                    OrderedDict([('y', self._dims['y']), ('x', self._dims['x'])]),
                    datatype=np.dtype(np.float32),
                    fillvalue=np.float32(1000),
                    units='degrees')
        if os.path.exists(self._get_mask_path('CLOUDS', 'B00')):
            self._fields['cloud_mask'] = Field(
                Variable(shortname='cloud_mask'),
                OrderedDict([('y', self._dims['y']), ('x', self._dims['x'])]),
                datatype=np.dtype(np.bool))
        if 'sun_zenith' in self._angle_grids:
            self._fields['sun_zenith_angle'] = Field(
                Variable(shortname='sun_zenith_angle'),
                OrderedDict([('y', self._dims['y']), ('x', self._dims['x'])]),
                datatype=np.dtype(np.float32),
                fillvalue=np.float32(1000),
                units='degrees')
        if 'sun_azimuth' in self._angle_grids:
            self._fields['sun_azimuth_angle'] = Field(
                Variable(shortname='sun_azimuth_angle'),
                OrderedDict([('y', self._dims['y']), ('x', self._dims['x'])]),
                datatype=np.dtype(np.float32),
                fillvalue=np.float32(1000),
                units='degrees')
        for fname in self._fields:
            self._fields[fname].attach_storage(self.get_field_handler(fname))
        return handler

    def close(self):
        """Close handler on storage"""
        if self.is_opened():
            self._handler.close()
            self._handler = None
        for band in self._band_handlers.keys():
            self._band_handlers[band] = None
            del self._band_handlers[band]
        for mask in self._mask_handlers.keys():
            self._mask_handlers[mask] = None
            del self._mask_handlers[mask]
        for angle in self._angle_handlers.keys():
            self._angle_handlers[angle] = None
            del self._angle_handlers[angle]

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
            # xy = np.zeros((shp[0] * shp[1], 2), dtype='float32')
            # xy[:, 1] = np.tile(self.read_values('y', slices=[slices[0]])[:, np.newaxis],
            #                    (1, shp[1])).flatten()
            # xy[:, 0] = np.tile(self.read_values('x', slices=[slices[1]])[np.newaxis, :],
            #                    (shp[0], 1)).flatten()
            # lonlat = proj(xy, inverse=True)
            # if fieldname == 'lon':
            #     values = lonlat[:, 0].astype(field.datatype).reshape(shp)
            # elif fieldname == 'lat':
            #     values = lonlat[:, 1].astype(field.datatype).reshape(shp)
        elif fieldname in ['x', 'y']:
            if fieldname == 'x':
                spacing = self.read_global_attribute('xdim')
                start = self.read_global_attribute('ulx') + spacing / 2.
            elif fieldname == 'y':
                spacing = self.read_global_attribute('ydim')
                start = self.read_global_attribute('uly') + spacing / 2.
            inds = slices.indices_array(dtype=field.datatype)[0] * spacing + start
            values = np.ma.array(inds, dtype=field.datatype)
        elif fieldname[4:] == 'digital_number':
            band = fieldname[0:3]
            handler = self._get_band_handler(band)
            bnd = handler.GetRasterBand(1)
            if self._overview_index is not None:
                bnd = bnd.GetOverview(self._overview_index)
            if self._tight == True:
                slices = slices.absolute_slices(self._tight_slices)
            offext = slices.indices_offsetextent()
            offset = [oe[0] for oe in offext]
            size = [oe[1] for oe in offext]
            vals = bnd.ReadAsArray(offset[1], offset[0], size[1], size[0])
            if slices[0].step != 1 or slices[1].step != 1:
                vals = vals[::slices[0].step, ::slices[1].step]
            # nodata = self.read_global_attribute('nodata_value')
            # satur = self.read_global_attribute('saturated_value')
            # values = np.ma.masked_where((vals == nodata) | (vals == satur),
            #                             vals, copy=False)
            values = np.ma.masked_equal(vals, field.fillvalue, copy=False)
        elif fieldname[4:] == 'detector_index':
            band = fieldname[0:3]
            source_ds = self._get_mask_handler('DETFOO', band)
            source_layer = source_ds.GetLayer(0)
            fillvalue = field.fillvalue
            if fillvalue != 0:
                fill = int(fillvalue)
            else:
                fill = None
            target_ds = self._get_gdal_dataset_from_slices(slices, fill=fill)
            gdal.RasterizeLayer(target_ds, [1], source_layer,
                                options=['ATTRIBUTE=detector_index'])
            vals = target_ds.GetRasterBand(1).ReadAsArray()
            target_ds = None
            if slices[0].step != 1 or slices[1].step != 1:
                vals = vals[::slices[0].step, ::slices[1].step]
            values = np.ma.masked_equal(vals, fillvalue, copy=False)
        elif fieldname[4:] in ['viewing_zenith_angle', 'viewing_azimuth_angle']:
            band = fieldname[0:3]
            angle_type = fieldname[4:].rstrip('_angle')
            det_ind = self.read_values('{}_detector_index'.format(band),
                                       slices=slices)
            hist, _ = np.histogram(det_ind.compressed(),
                                   np.arange(1, 14))
            vals = np.zeros(det_ind.shape, dtype=field.datatype)
            fillvalue = field.fillvalue
            vals.fill(fillvalue)
            for idet in range(1, 13):
                if hist[idet - 1] == 0:
                    continue
                ind = np.where(det_ind == idet)
                source_ds = self._get_angle_handler(angle_type, band=band,
                                                    detector=idet,
                                                    nodatavalue=float(fillvalue))
                indmin = [i.min() for i in ind]
                indmax = [i.max() + 1 for i in ind]
                rel_slices = Slices([slice(indmin[0], indmax[0]),
                                     slice(indmin[1], indmax[1])],
                                    det_ind.shape)
                abs_slices = rel_slices.absolute_slices(slices)
                target_ds = self._get_gdal_dataset_from_slices(abs_slices,
                                                               btype=gdal.GDT_Float32,
                                                               fill=float(fillvalue),
                                                               use_step=True)
                gdal.ReprojectImage(source_ds, target_ds, None, None,
                                    gdal.GRA_Bilinear)
                detvals = target_ds.GetRasterBand(1).ReadAsArray()
                target_ds = None
                vals[ind] = detvals[ind[0] - indmin[0], ind[1] - indmin[1]]
                # target_ds = self._get_gdal_dataset_from_slices(slices,
                #                                                btype=gdal.GDT_Float32,
                #                                                fill=float(fillvalue),
                #                                                use_step=True)
                # gdal.ReprojectImage(source_ds, target_ds, None, None,
                #                     gdal.GRA_Bilinear)
                # detvals = target_ds.GetRasterBand(1).ReadAsArray()
                # target_ds = None
                # vals[ind] = detvals[ind]
            values = np.ma.masked_equal(vals, fillvalue, copy=False)
            # TEST
            # band = fieldname[0:3]
            # angle_type = fieldname[4:].rstrip('_angle')
            # fillvalue = field.fillvalue
            # source_ds = self._get_angle_handler(angle_type, band=band,
            #                                     detector=12,
            #                                     nodatavalue=float(fillvalue))
            # target_ds = self._get_gdal_dataset_from_slices(slices,
            #                                                btype=gdal.GDT_Float32,
            #                                                fill=float(fillvalue),
            #                                                use_step=True)
            # gdal.ReprojectImage(source_ds, target_ds, None, None,
            #                     gdal.GRA_Bilinear)
            # vals = target_ds.GetRasterBand(1).ReadAsArray()
            # target_ds = None
            # values = np.ma.masked_equal(vals, fillvalue, copy=False)
            # import matplotlib.pyplot as plt
            # det_ind = self.read_values(band+'_detector_index')
            # plt.figure() ; plt.imshow(det_ind, interpolation='nearest') ; plt.colorbar()
            # source = np.ma.masked_equal(source_ds.GetRasterBand(1).ReadAsArray(), fillvalue)
            # plt.figure() ; plt.imshow(source, interpolation='nearest') ; plt.colorbar()
            # plt.figure() ; plt.imshow(values, interpolation='nearest') ; plt.colorbar()
            # plt.show()
            #import pdb ; pdb.set_trace()
            # TEST
        elif fieldname == 'cloud_mask':
            source_ds = self._get_mask_handler('CLOUDS', 'B00')
            source_layer = source_ds.GetLayer(0)
            if source_layer is None:
                return np.zeros(slices.shape(), dtype='bool')
            target_ds = self._get_gdal_dataset_from_slices(slices)
            gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[1],
                                options=['ALL_TOUCHED=TRUE'])
            values = target_ds.GetRasterBand(1).ReadAsArray().astype('bool')
            target_ds = None
            if slices[0].step != 1 or slices[1].step != 1:
                values = values[::slices[0].step, ::slices[1].step]
        elif fieldname in ['sun_zenith_angle', 'sun_azimuth_angle']:
            angle_type = fieldname.rstrip('_angle')
            fillvalue = field.fillvalue
            source_ds = self._get_angle_handler(angle_type,
                                                nodatavalue=float(fillvalue))
            target_ds = self._get_gdal_dataset_from_slices(slices,
                                                           btype=gdal.GDT_Float32,
                                                           fill=float(fillvalue),
                                                           use_step=True)
            gdal.ReprojectImage(source_ds, target_ds, None, None,
                                gdal.GRA_Bilinear)
            vals = target_ds.GetRasterBand(1).ReadAsArray()
            target_ds = None
            values = np.ma.masked_equal(vals, fillvalue, copy=False)
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
        strtime = self.read_global_attribute('sensing_time')
        return datetime.strptime(strtime, '%Y-%m-%dT%H:%M:%S.%fZ')

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage.

        Returns:
            datetime: end time of the data in file.
        """
        strtime = self.read_global_attribute('sensing_time')
        return datetime.strptime(strtime, '%Y-%m-%dT%H:%M:%S.%fZ')

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


class SAFEMSIL1C10mGranuleFile(SAFEMSIL1CGranuleFile):
    """Mapper class for SAFE MSI L1C 10m granule files.
    """
    def __init__(self, url, **kwargs):
        super(SAFEMSIL1C10mGranuleFile, self).__init__(url, '10m', **kwargs)


class SAFEMSIL1C20mGranuleFile(SAFEMSIL1CGranuleFile):
    """Mapper class for SAFE MSI L1C 20m granule files.
    """
    def __init__(self, url, **kwargs):
        super(SAFEMSIL1C20mGranuleFile, self).__init__(url, '20m', **kwargs)


class SAFEMSIL1C60mGranuleFile(SAFEMSIL1CGranuleFile):
    """Mapper class for SAFE MSI L1C 60m granule files.
    """
    def __init__(self, url, **kwargs):
        super(SAFEMSIL1C60mGranuleFile, self).__init__(url, '60m', **kwargs)
