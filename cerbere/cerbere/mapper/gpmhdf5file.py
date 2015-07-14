# encoding: utf-8
"""
cerbere.mapper.GPMHDF5File
==========================

Mapper class for JAXA GPM HDF5 files

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
from collections import OrderedDict
from datetime import datetime

import numpy

import cerbere.mapper.slices
from cerbere.mapper import hdf5file
from cerbere.datamodel import field
from cerbere.datamodel import variable

REFERENCE_TIME = datetime(2001, 1, 1, 0, 0, 0)

MS = "MS"
HS = "HS"

class GPMHDF5File(hdf5file.HDF5File):
    """Mapper class for JAXA GPM HDF5 files"""

    def __init__(self,
                 url=None,
                 mode=hdf5file.READ_ONLY,
                 subproduct=None,
                 **kwargs):
        """
        Args:
            subproduct (str): sub-product to open. Used only in the case of
                Ka L2 product which contains two underlying products with their
                own geolocation grid.
        """
        hdf5file.HDF5File.__init__(self, url=url, mode=mode, **kwargs)
        self.attributes = None
        self.__dim2name = None
        self.__name2dim = None
        self._subproduct = subproduct
        return

    def close(self):
        self.attributes = None
        self.__dim2name = None
        self.__name2dim = None

    def open(self,
             view=None,datamodel_geolocation_dims=None,datamodel=None):
        """Open the HDF5 file

        Args:
            view (dict, optional): a dictionary where keys are dimension names
                and values are slices. A view can be set on a file, meaning
                that only the subset defined by this view will be accessible.
                This view is expressed as any subset (see :func:`get_values`).
                For example::

                view = {'time':slice(0,0), 'lat':slice(200,300),
                'lon':slice(200,300)}

        Returns:
            an handler on the opened file
        """
        handler = hdf5file.HDF5File.open(self, view,datamodel_geolocation_dims=datamodel_geolocation_dims,datamodel=datamodel)
        if (self.read_global_attribute('AlgorithmID') == '2AKa' and
                self._subproduct is None):
            raise Exception("The subproduct has not been specified when"
                            "opening the file.")
        return handler

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
        if fieldname in ['lat', 'lon']:
            if self.read_global_attribute('AlgorithmID') == '2AKu':
                return {'lat': '/NS/Latitude',
                        'lon': '/NS/Longitude'}[fieldname]
            elif self.read_global_attribute('AlgorithmID') == '2AKa':
                return ('/%s/' % self._subproduct
                        + {'lat': 'Latitude',
                           'lon': 'Longitude'}[fieldname])
        else:
            return fieldname

    def get_matching_dimname(self, dimname):
        """Return the equivalent name in the native format for a standard
        dimension.

        This is a translation of the standard names to native ones. It is used
        for internal purpose only and should not be called directly.

        The standard dimension names are row, cell, time for
        :class:`~cerbere.datamodel.swath.Swath`

        Args:
            dimname (str): standard dimension name.

        Returns:
            str: return the native name for the dimension. Return `dimname` if
                the input dimension has no standard name.

        See Also:
            see :func:`get_standard_dimname` for the reverse operation
        """
        #call get_dimensions to be sure __name2dim is not None
        self.get_dimensions()
        if dimname in self.__name2dim:
            return self.__name2dim[dimname]
        else:
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
        if self.read_global_attribute('AlgorithmID') == '2AKu':
            match = {'nscan': 'row',
                     'nray': 'cell'}
        elif self.read_global_attribute('AlgorithmID') == '2AKa':
            if self._subproduct == "HS":
                match = {'nscan': 'row',
                         'nrayHS': 'cell'}
            else:
                match = {'nscan': 'row',
                         'nrayMS': 'cell'}
        else:
            raise Exception("Unknown product type: %s",
                            self.read_global_attribute('AlgorithmID'))
        if dimname in match:
            return match[dimname]
        else:
            return dimname

    def get_fieldnames(self):
        """Returns the list of geophysical fields stored for the feature.

        The geolocation field names are excluded from this list.

        Returns:
            list<string>: list of field names
        """
        fieldnames = hdf5file.HDF5File.get_fieldnames(self)
        if '/AlgorithmRuntimeInfo' in fieldnames:
            fieldnames.remove('/AlgorithmRuntimeInfo')
        if self.read_global_attribute('AlgorithmID') == '2AKa':
            filtered_fieldnames = []
            for fieldname in fieldnames:
                if self._subproduct == "HS" and '/MS/' not in fieldname:
                    filtered_fieldnames.append(fieldname)
                elif self._subproduct == "MS" and '/HS/' not in fieldname:
                    filtered_fieldnames.append(fieldname)
            return filtered_fieldnames
        else:
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
        # fill in dictionaries for translation
        if self.__dim2name is None and self.__name2dim is None:
            self.__dim2name = {}
            self.__name2dim = {}
            dims = hdf5file.HDF5File.get_dimensions(self)
            newdims = []
            for dim in dims:
                rank = int(dim[-1])
                dataset = dim[:-1]
                dimnames = (self.get_handler()[dataset]
                            .attrs.get('DimensionNames'))
                if dimnames is not None:
                    dimnames = dimnames.split(',')
                    newdim = self.get_standard_dimname(dimnames[rank])
                    self.__dim2name[dim] = newdim
                    self.__name2dim[newdim] = dim
                    if dimnames[rank] not in newdims:
                        newdims.append(newdim)
        if fieldname is None:
            return self.__name2dim.keys()
        else:
            # get field native name with full group path
            native_fieldname = self._get_native_fieldname(fieldname)
            if (native_fieldname not in self.fieldnames and
                    native_fieldname not in ['lat', 'lon', 'time']):
                raise Exception("Unknown fieldname in this product")
            # get dimension names from attribute
            attrdims = self.get_handler()[native_fieldname].attrs.get(
                'DimensionNames'
                )
            # replace with standard names
            attrdims = [self.get_standard_dimname(dim)
                        for dim in attrdims.split(',')]
            return tuple(attrdims)

    def read_global_attributes(self):
        """Returns the names of the global attributes.

        Returns:
            list<str>: the list of the attribute names.
        """
        if self.attributes is not None:
            return self.attributes.keys()
        else:
            self.attributes = OrderedDict([])
        grossattrs = self.get_handler().attrs.values()
        for attr in grossattrs:
            items = attr.strip(';\n').split(';\n')
            for item in items:
                key, val = item.split('=')
                self.attributes[key] = val
        return self.attributes.keys()

    def read_global_attribute(self, attr):
        """Returns the value of a global attribute.

        Args:
            name (str): name of the global attribute.

        Returns:
            str, number or datetime: value of the corresponding attribute.
        """
        if self.attributes is None:
            self.read_global_attributes()
        return self.attributes[attr]

    def read_field_attributes(self, fieldname):
        """Return the specific attributes of a field.

        Args:
            fieldname (str): name of the field.

        Returns:
            dict<string, string or number or datetime>: a dictionary where keys
                are the attribute names.
        """
        attributes = hdf5file.HDF5File.read_field_attributes(self, fieldname)
        newattrs = {}
        for attr in attributes:
            # remove standard attributes
            if attr not in ['Units',
                            'DimensionNames',
                            'CodeMissingValue']:
                items = attributes[attr].strip(';\n').split(';\n')
                if len(items) == 1:
                    newattrs[attr] = attributes[attr]
                else:
                    for item in items:
                        key, val = item.split('=')
                        newattrs[key] = val
        return newattrs

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
        if fieldname == 'time':
            # create a field for time
            vartime = variable.Variable(
                shortname='time',
                description='acquisition time',
                authority=self.get_naming_authority(),
                standardname='time'
                )
            fieldobj = field.Field(
                vartime,
                self.get_full_dimensions('lat'),
                datatype=numpy.dtype(numpy.int64),
                units='seconds since 2001-01-01 00:00:00'
                )
            fieldobj.attach_storage(self.get_field_handler(fieldname))
            return fieldobj
        fieldobj = hdf5file.HDF5File.read_field(self, fieldname)
        # get field native name
        native_fieldname = self._get_native_fieldname(fieldname)
        # get unit and fill value
        attrs = self.get_handler()[native_fieldname].attrs.keys()
        if 'Units' in attrs:
            fieldobj.units = (self.get_handler()[native_fieldname]
                              .attrs.get('Units'))
        if 'CodeMissingValue' in attrs:
            fillvalue = (self.get_handler()[native_fieldname]
                         .attrs.get('CodeMissingValue'))
            fillvalue = numpy.float(fillvalue)
            fieldobj.fillvalue = fillvalue
            if fieldobj._values is not None:
                self._values.set_fill_value(fillvalue)
        return fieldobj

    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage.

        Returns:
            datetime: start time of the data in file.
        """
        attr = self.read_global_attribute('StartGranuleDateTime')
        return datetime.strptime(attr[:-5], "%Y-%m-%dT%H:%M:%S")

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage.

        Returns:
            datetime: end time of the data in file.
        """
        attr = self.read_global_attribute('StopGranuleDateTime')
        return datetime.strptime(attr[:-5], "%Y-%m-%dT%H:%M:%S")

    def read_values(self, fieldname, slices=None):
        """Read the data of a field.

        Args:
            fieldname (str): name of the field which to read the data from

            slices (list of slice, optional): list of slices for the field if
                subsetting is requested. A slice must then be provided for each
                field dimension. The slices are relative to the opened view
                (see :func:open) if a view was set when opening the file.

        Return:
            MaskedArray: array of data read. Array type is the same as the
                storage type.
        """
        if fieldname == 'time':
            # adjust slice with respect to view
            # time is only defined along-track
            dims = self.get_full_dimensions('lat')
            newslices = cerbere.mapper.slices.get_absolute_slices(
                view=self.view,
                slices=slices,
                dimnames=dims.keys(),
                dimsizes=dims.values()
                )
            shape = cerbere.mapper.slices.get_shape_from_slice(newslices)
            newslices = [newslices[0]]
            if self.read_global_attribute('AlgorithmID') == '2AKu':
                substr = '/NS/'
            elif self.read_global_attribute('AlgorithmID') == '2AKa':
                substr = '/%s/' % self._subproduct
            year = hdf5file.HDF5File.read_values(self,
                                                 substr + 'ScanTime/Year',
                                                 newslices)
            day = hdf5file.HDF5File.read_values(self,
                                                substr + 'ScanTime/DayOfYear',
                                                newslices)
            sec = hdf5file.HDF5File.read_values(self,
                                                substr + 'ScanTime/SecondOfDay',
                                                newslices)
            year = year - REFERENCE_TIME.year + 2001
            minyear = year.min()
            maxyear = year.max()
            if minyear != maxyear:
                yeardiff = numpy.array(year.shape, dtype=numpy.int32)
                yeardiff[yeardiff == minyear] = (
                    datetime(minyear, 1, 1, 0, 0, 0)
                    - REFERENCE_TIME).total_seconds()
                yeardiff[yeardiff == maxyear] = (
                    datetime(maxyear, 1, 1, 0, 0, 0)
                    - REFERENCE_TIME).total_seconds()
            else:
                yeardiff = (datetime(minyear, 1, 1, 0, 0, 0)
                            - REFERENCE_TIME).total_seconds()
            total_seconds = yeardiff + (day-1) * 86400. + sec
            # reshape to (row, cell) array
            total_seconds = numpy.resize(total_seconds,
                                         (shape[1], shape[0])
                                         ).transpose()

            return total_seconds
        else:
            return hdf5file.HDF5File.read_values(self, fieldname, slices)

    def get_product_version(self):
        """return the product version"""
        return self.read_global_attribute('ProductVersion')

    def get_orbit_number(self):
        """In the case of a satellite orbit file, returns the orbit number.

        Returns:
            int: the orbit number
        """
        return self.read_global_attribute('GranuleNumber')
    
    def get_collection_id(self):
        if 'KAR' in self._url:
            res = 'NASA-JAXA-RAIN-GLO-30MIN-001-GPM-Ka'
        elif 'KUR' in self._url:
            res = 'NASA-JAXA-RAIN-GLO-30MIN-001-GPM-Ku'
        return res
    
    def get_cycle_number(self):
        res = self._url.split('_')[4]
        return res