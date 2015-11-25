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

    def save(self, output=None, attrs=None):
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
        if output is None:
            mapper = self.get_mapper()
        else:
            mapper = output
        if mapper.is_writable():

            # creating dimensions
            mapper.create_dim('time', None)

            # create additional dimensions
            dims = ['time']
            for v in self._fields.keys():
                if v not in self._geolocation_fields:
                    for d in self._fields[v].dimensions:
                        if d not in dims:
                            mapper.create_dim(d, self._fields[v].get_dimsize(d))
                            dims.append(d)

            # creating metadata
            if attrs:
                globalattr = attrs
            else:
                globalattr = {}
            if 'title' not in globalattr and self.title is not None:
                globalattr['title'] = self.title
            if 'summary' not in globalattr and self.description is not None:
                globalattr['summary'] = self.description
            if 'id' not in globalattr and self.identifier is not None:
                globalattr['id'] = self.identifier
            lonmin, latmin, lonmax, latmax = self.get_bbox()
            globalattr['geospatial_lat_min'] = latmin
            globalattr['geospatial_lon_min'] = lonmin
            globalattr['geospatial_lat_max'] = latmax
            globalattr['geospatial_lon_max'] = lonmax
            if 'time_coverage_start' not in globalattr:
                tmptime = self.get_start_time()
                if tmptime is not None:
                    globalattr['time_coverage_start'] = tmptime
            if 'time_coverage_end' not in globalattr:
                tmptime = self.get_end_time()
                if tmptime is not None:
                    globalattr['time_coverage_end'] = tmptime
            globalattr['cdm_data_type'] = 'trajectory'
            mapper.write_global_attributes(globalattr)

            # creating fields
            if self._geolocation_fields['time'] is None:
                raise Exception('No time information defined')
            for geof in self._geolocation_fields:
                if self._geolocation_fields[geof] is None:
                    logging.warning('Missing geolocation variable : %s', geof)
                else:
                    mapper.create_field(self._geolocation_fields[geof])
            for dataf in self._fields.keys():
                if dataf not in self._geolocation_fields:
                    mapper.create_field(self._fields[dataf])
        else:
            raise Exception('Mapper object is not writable')

        # saving fields
        for geof in self._geolocation_fields:
            field = self._geolocation_fields[geof]
            if field is not None and not field.is_saved():
                mapper.write_field(self._geolocation_fields[geof])
        for dataf in self._fields.keys():
            if not self._fields[dataf].is_saved():
                mapper.write_field(self._fields[dataf])
        mapper.sync()
        return
