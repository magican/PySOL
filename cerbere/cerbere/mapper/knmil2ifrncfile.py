# -*- coding: utf-8 -*-
"""
cerbere.mapper.ascatifrncfile
============================

Mapper class for IFREMER ASCAT / OSCAT netcdf files (converted from BUFR)

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
from datetime import datetime
from collections import OrderedDict

import numpy
from .. import READ_ONLY, WRITE_NEW
from .ncfile import NCFile
from ..datamodel.field import Field
from ..datamodel.variable import Variable

VIRTUALFIELD_DESCR = {'wind_speed_selection': 'speed',
                      'wind_dir_selection': 'direction'
                      }

VIRTUALFIELD_UNITS = {
                'wind_speed_selection': 'm s-1',
                'wind_dir_selection': 'degrees'
                }

VIRTUALFIELD_MATCHING = {
    'OSCATL2B-050': {
        'wind_speed_selection': 'Wind_speed_at_10_m__%d',
        'wind_dir_selection': 'Wind_direction_at_10_m__%d'
        },
    'ASCATL2B-050': {
        'wind_speed_selection': 'Wind_speed_at_10_m_%d',
        'wind_dir_selection': 'Wind_direction_at_10_m_%d'
        }
    }


class KNMIL2IFRNCFile(NCFile):
    """Mapper class for IFREMER KNMI netcdf files (converted from BUFR)

    The IFREMER ASCAT & OSCAT files are in netCDF format but not cf compliant
    (for instance for time). This class is inherited from the
    :class:`NCFile` mapper but extended to mimic the presence of a `time`
    variable. Last, it adds virtual variables to provide the best
    wind speed and direction solutions.
    """
    def __init__(self, url=None, mode=READ_ONLY,
                 ncformat='NETCDF4_CLASSIC', **kwargs):
        super(KNMIL2IFRNCFile, self).__init__(url=url,
                                              mode=mode,
                                              ncformat=ncformat,
                                              **kwargs)

    def get_geolocation_field(self, fieldname):
        matching = {'time': 'time',
                    'lon': 'Longitude',
                    'lat': 'Latitude'}
        if fieldname in matching:
            return matching[fieldname]
        else:
            return None

    def get_matching_dimname(self, geodimname):
        matching = {'time': 'time',
                    'row': 'measurement',
                    'cell': 'swath'}
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None

    def get_standard_dimname(self, geodimname):
        matching = {'time': 'time',
                    'measurement': 'row',
                    'swath': 'cell'}
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None

    def get_fieldnames(self):
        """Returns the list of geophysical fields stored for the feature"""
        fieldnames = super(KNMIL2IFRNCFile, self).get_fieldnames()
        # add virtual variables for best wind speed and direction
        fieldnames.extend(['wind_speed_selection', 'wind_dir_selection'])
        return fieldnames

    def read_field(self, fieldname):
        """
        Return the field, without its values.

        Actual values can be retrieved with read_values() method.
        """
        # special implementation case for time field which is not
        # available as a variable in IFREMER KNMI files
        if fieldname == 'time':
            # create a field for time
            variable = Variable(
                    shortname=fieldname,
                    description='time',
                    authority=self.get_naming_authority(),
                    standardname='time'
                    )
            field = Field(
                    variable,
                    OrderedDict([('row', self.get_dimsize('row')),
                                 ('cell', self.get_dimsize('cell'))
                                 ]),
                    datatype=numpy.dtype(numpy.float32),
                    units='seconds since 1990-01-01 00:00:00'
                    )
            field.attach_storage(self.get_field_handler(fieldname))
        elif fieldname in ['wind_speed_selection', 'wind_dir_selection']:
            # create a virtual field
            variable = Variable(
                    shortname=fieldname,
                    description='best wind %s solution'
                        % VIRTUALFIELD_DESCR[fieldname],
                    authority=self.get_naming_authority(),
                    standardname='wind %s' % VIRTUALFIELD_DESCR[fieldname]
                    )
            field = Field(
                    variable,
                    OrderedDict([('row', self.get_dimsize('row')),
                                 ('cell', self.get_dimsize('cell'))
                                 ]),
                    datatype=numpy.dtype(numpy.float32),
                    units=VIRTUALFIELD_UNITS[fieldname]
                    )
            field.attach_storage(self.get_field_handler(fieldname))
        else:
            field = NCFile.read_field(self, fieldname)
        return field

    def read_values(self, fieldname, slices=None):
        """
        """
        if fieldname == 'time':
            if self._handler is None:
                self._handler = self.get_handler()
            if slices:
                datevar = self.get_handler().variables['Date'][slices]
            else:
                datevar = self.get_handler().variables['Date'][:]
            if self.get_handler().variables['Date'].units\
                     != "julian day since 1990-01-01 0:0:0":
                raise Exception("Inconsistent start date")
            if slices:
                millisecvar = self.get_handler().variables['Millisec'][slices]
            else:
                millisecvar = self.get_handler().variables['Millisec'][:]

            millisecvar = millisecvar.astype('float64')
            datevar = datevar.astype('float64')
            times = datevar * 86400. + millisecvar / 1000.
            if slices is None:
                return times
            else:
                return times[slices]
        elif fieldname in ['wind_speed_selection', 'wind_dir_selection']:
            idxname = 'Index_of_selected_wind_vector'
            if slices:
                idx = self.get_handler().variables[idxname][slices]
            else:
                idx = self.get_handler().variables[idxname][:]
            varname = VIRTUALFIELD_MATCHING[self.get_collection_id()][fieldname]
            solutions = []
            for sol in range(4):
                solname = varname % sol
                if solname in self.get_handler().variables:
                    solutions.append(self.get_handler().variables[solname])
            values = numpy.ma.masked_all(
                        (self.get_dimsize('row'), self.get_dimsize('cell')),
                        dtype=numpy.dtype(numpy.float32)
                        )
            if slices is not None:
                values = values[slices]
            for sol in range(len(solutions)):
                ind = numpy.ma.where(idx == sol + 1)
                if slices:
                    solvalues = solutions[sol][slices]
                else:
                    solvalues = solutions[sol][:]
                values[ind] = solvalues[ind]
            return values
        else:
            return super(KNMIL2IFRNCFile, self).read_values(fieldname,
                                                            slices=slices)

    def get_start_time(self):
        """Return start of temporal coverage"""
        start = datetime.strptime(
            self.get_handler().Start_Date,
            '%Y%m%dT%H%M%S'
            )
        return start

    def get_end_time(self):
        """Return end of temporal coverage"""
        end = datetime.strptime(
            self.get_handler().Stop_Date,
            '%Y%m%dT%H%M%S'
            )
        return end

    def get_bbox(self):
        """
        return the bounding box of the feature, as a tuple
        (lonmin, latmin, lonmax, latmax)
        """
        return (float(self.get_handler().West_Longitude),
                float(self.get_handler().South_Latitude),
                float(self.get_handler().East_Longitude),
                float(self.get_handler().North_Latitude)
                )

    def get_collection_id(self):
        """return the identifier of the product collection"""
        return self.get_handler().Sensor_Id +\
                self.get_handler().Product_Id + "-050"
