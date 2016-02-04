# -*- coding: utf-8 -*-
"""
cerbere.datamodel.trajectory
============================

Model class for the trajectory feature

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import logging

from netCDF4 import num2date

from .abstractfeature import AbstractFeature

__all__ = ['Trajectory']

class Trajectory(AbstractFeature):
    """
    Model class for the trajectory feature
    """
    def __init__(self, identifier=None, title=None, 
                description=None, source=None, metadata=None,
                fields=None,
                longitudes=None, latitudes=None, depths=None, times=None):
        super(Trajectory, self).__init__(
                            identifier=identifier,
                            title=title,
                            description=description,
                            source=source,
                            metadata=metadata,
                            fields=fields,
                            longitudes=longitudes,
                            latitudes=latitudes,
                            depths=depths,
                            times=times
                            )
        return


    def get_geolocation_dimnames(self):
        """
        Returns the geolocation dimension names defining a swath
        """
        return ['time']

    def get_geolocation_field_dimnames(self, dimname):
        """
        Returns the dimensions of a geolocation field
        """
        return ('time',)


    def get_geolocation_dimsizes(self):
        """
        Returns the geolocation dimension sizes of the trajectory object
        """
        return {'time':self.get_geolocation_field('time').get_dimsize('time')}

    def extract_subset(self,
                       boundaries=None,
                       slices=None,
                       fields=None,
                       padding=False):
        """Extract a subset feature from the trajectory.

        The created subset is a new Trajectory object without any reference to
        the source.

        Args:
            boundaries (tuple): area of the subset to extract, defined as
                llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat.

            slices (dict): indices for /time/ dimension of the subset to
                extract from the source data. If None, the complete feature
                is extracted.

            fields (list): list of field names to extract. If None, all fields
                are extracted.

            padding (bool): Passed to extract_field method to ensure padding
                with _FillValues for points outside of the bounds of the this
                feature (used only in conjuncture with slices.
        """
        if boundaries and slices:
            raise Exception("Boundaries and slices can not be both provided.")
        if boundaries:
            # get corresponding slices
            raise NotImplementedError
        elif slices:
            if len(slices) != 1:
                raise Exception("only slices for the time dimension must be"
                                "provided")
        subtraj = Trajectory(
            latitudes=self.extract_field('lat',
                                         slices=slices,
                                         padding=padding),
            longitudes=self.extract_field('lon',
                                          slices=slices,
                                          padding=padding),
            times=self.extract_field('time',
                                     slices=slices,
                                     padding=padding),
            metadata=self.metadata)
        if fields is None:
            fields = self.get_fieldnames()
        for field in fields:
            subtraj.add_field(self.extract_field(field,
                                                 slices=slices,
                                                 padding=padding))
        return subtraj

    def save(self, output, attrs=None):
        """
        Save the trajectory to a storage (file,...).

        Args:
            output (:class:`~cerbere.mapper.abstractmapper.AbstractMapper`): storage
                object where to save the feature data.

            attrs (dict): the global metadata (attributes) of the feature, as
                a dictionary where keys are the attributes names.
                See STANDARD_ATTRIBUTES_VALUES in abstractmapper class to see
                a list of standard attributes
        """
        if output.is_writable():

            # creating dimensions
            output.create_dim( 'time', None )

            # creating metadata
            if attrs:
                globalattr = attrs
            else:
                globalattr = {}
            if not self._bbox:
                self._bbox = (
                    self.get_lon().min(), self.get_lat().min(),
                    self.get_lon().max(), self.get_lat().max(),
                    )
            globalattr['geospatial_lat_min'] = self._bbox[1]
            globalattr['geospatial_lon_min'] = self._bbox[0]
            globalattr['geospatial_lat_max'] = self._bbox[3]
            globalattr['geospatial_lon_max'] = self._bbox[2]
            if 'time_coverage_start' not in globalattr:
                times = self.get_times()
                mintime = num2date(times.min(),
                                   self.get_time_units())
                globalattr['time_coverage_start'] = (
                    mintime.strftime('%Y%m%dT%H%M%S'))
            if 'time_coverage_end' not in globalattr:
                maxtime = num2date(times.max(),
                                   self.get_time_units())
                globalattr['time_coverage_end'] = (
                    maxtime.strftime('%Y%m%dT%H%M%S'))
            globalattr['cdm_data_type'] = 'trajectory'
            output.write_global_attributes(globalattr)

            # creating fields
            for geof in self._geolocation_fields:
                if not self._geolocation_fields['time']:
                    raise Exception('No time information defined')
                if not self._geolocation_fields[geof]:
                    logging.debug('Missing geolocation variable : %s', geof)
                else:
                    output.create_field(self._geolocation_fields[geof])
            for dataf in self._fields.keys():
                output.create_field(self._fields[dataf])
        else:
            raise Exception("Mapper object is not writable")

        # saving fields
        for geof in self._geolocation_fields:
            if self._geolocation_fields[geof]:
                output.write_field(self._geolocation_fields[geof])
        for dataf in self._fields.keys():
            output.write_field(self._fields[dataf])
        output.sync()
