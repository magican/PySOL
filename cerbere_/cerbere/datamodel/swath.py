# -*- coding: utf-8 -*-
"""
cerbere.datamodel.swath
============================

Model class for the swath feature

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import logging
import collections
from netCDF4 import num2date

import numpy

from .abstractfeature import AbstractFeature
from ..geo.point import Point
#from ..utils.interpolation import Interpolation

__all__ = ['Swath']


class Swath(AbstractFeature):
    """
    Model class for the swath feature, ie a two-dimensional grid along
    the satellite track
    """
    def __init__(self, identifier=None, title=None,
                description=None, source=None, metadata=None,
                fields=None,
                longitudes=None, latitudes=None, depths=None, times=None,
                spatial_resolution=None):
        super(Swath, self).__init__(
                identifier=identifier,
                title=title,
                description=description,
                source=source,
                metadata=metadata,
                fields=fields,
                longitudes=longitudes,
                latitudes=latitudes,
                depths=depths,
                times=times,
                spatial_resolution=spatial_resolution
                )
        return

    def get_geolocation_dimnames(self):
        """
        Returns the geolocation dimension names defining a swath
        """
        return ['row', 'cell']

    def get_geolocation_field_dimnames(self, fieldname):
        """
        Returns the dimensions of a geolocation field
        """
        if fieldname == 'depth':
            return ('depth',)
        else:
            return ('row', 'cell',)

    def get_geolocation_dimsizes(self):
        """
        Returns the geolocation dimension sizes of the swath object
        """
        row, cell = self.get_lat().shape
        return collections.OrderedDict([('row', row), ('cell', cell)])

#     def _get_field_values(self, field, slices=None, indices=None, cache=False,
#                           interpolated=False, **kwargs):
#         # remove time dimension
# #         if 'time' in field.dimensions:
# #             if indices:
# #                 indices['time'] = 0
# #             else:
# #                 if slices is None:
# #                     slices = {}
# #                 slices['time'] = 0
#         values = field.get_values(slices, indices, cache, **kwargs)
#         if interpolated:
#             values = Interpolation.interpolate2d(values)
#         return values

    def extract_subset(
            self, boundaries=None, slices=None, fields=None, padding=False):
        """Extract a subset feature from the swath.

        The created subset is a new Swath object without any reference to
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
            raise Exception("Not yet implemented")
        elif slices:
            if len(slices) != 2:
                raise Exception(
                    "slices for the two spatial dimensions of "
                    "the swath must be provided"
                )
        sub_swath = Swath(
            latitudes=self.extract_field('lat', slices=slices, padding=padding),
            longitudes=self.extract_field(
                'lon', slices=slices, padding=padding
            ),
            times=self.extract_field('time', slices=slices, padding=padding),
            depths=self.extract_field('z', slices=slices, padding=padding),
            metadata=self.metadata
        )
        if fields is None:
            fields = self.get_fieldnames()
        elif not type(fields) is list:
            raise Exception("fields must be a list")
        for field in fields:
            sub_swath.add_field(
                self.extract_field(field, slices=slices, padding=padding)
            )
        return sub_swath

    def get_spatial_resolution(self):
        """Return the spatial resolution of the feature"""
        if self.spatial_resolution is None:
            i = self.get_geolocation_dimsizes()['cell'] / 2
            j = self.get_geolocation_dimsizes()['row'] / 2
            res = Point.get_distance(
                    Point(self.get_lon(
                                slices={'row': slice(j, j + 1, 1),
                                        'cell': slice(i, i + 1, 1)}
                                ),
                          self.get_lat(
                                slices={'row': slice(j, j + 1, 1),
                                        'cell': slice(i, i + 1, 1)}
                                )),
                    Point(self.get_lon(
                                {'row': slice(j + 1, j + 2, 1),
                                 'cell': slice(i + 1, i + 2, 1)}
                                ),
                          self.get_lat(
                                {'row': slice(j + 1, j + 2, 1),
                                 'cell': slice(i + 1, i + 2, 1)}
                                ))
                    ).ravel()[0] / numpy.sqrt(2)
            return res

    def extract_section(self, lat1, lon1, lat2, lon2, fieldnames=None):
        """
        """
        x, y = np.linspace(x0, x1, num), np.linspace(y0, y1, num)

        # Extract the values along the line, using cubic interpolation
        zi = scipy.ndimage.map_coordinates(z, np.vstack((x,y)))
        
        return result

    def is_xy_in_grid(self, x, y, varName=None):
        """Renvoie False si la valeur est hors zone.
        Dans le cas où un nom de variable (varName) est précisé,
        renvoie également False si la valeur à cette position est masquée.
        @param x:indice de longitude
        @param y:indice de latitude
        @param varName: string, optionnel, nom d'une des variables, pour obtenir une précision supplémentaire :
        l'état (masqué ou non) de la valeur (x,y)"""
        
        if (x>=0) and (y>=0) and (x<len(self.lonValues)) and (y<len(self.latValues)) :
            if varName is None:
                return True
            else:
                # Vérification plus poussée
                val = self.recordVOs[varName].getValues(x=x,y=y)
                if isinstance(val, numpy.ma.MaskedArray):
                    return not val.mask
                # Pour une version récente de numpy, le type numpy.ma.core.MaskedConstant serait également à vérifier
                else:
                    return True
        else:
            return False
        
    def is_latlon_in_grid(self, lon, lat):
        """
        Indique si le point de coordonnées (lon, lat), donné dans le repère de la grille, 
        appartient à l'espace qu'elle délimite.
        Tient compte de la rotondité de la Terre.
        
        La grille est vue comme un maillage, donc pour une grille aux longitudes sur [0 359],
        le point 359,5 est au dehors. Le point 380 est quant à lui inclus.
        
        Voir aussi la fonction intersection()
        @param lon:longitude
        @param lat:latitude
        @return: True si le point est à l'intérieur de la grille (vue comme un maillage)
        """
        # Sous forme de grille de cases, le point 359,5 serait dans la case centrée sur 359
        # si la résolution est supérieure à 1° (la dernière case s'étend au moins de 358,5 à 359,5
        
        
        # 360° en longitudes, donc le point à 400° retombe sur 40°
        lon = lon % 360
        selflon = self.lonValues.values[[0,-1]]  # Premier et dernier élément
        selflat = self.latValues.values[[0,-1]]
        
        # lonmin est rapporté sur [0, 360[.
        lonmin = selflon[0] % 360
        # lonmax est rapporté sur [lonmin, lonmin + 360[.
        lonmax = selflon[1] % 360
        if lonmax<lonmin:
            lonmax += 360
        # tout comme lon.
        if lon<lonmin:
            lon += 360
        # Les latitudes n'ont pas d'équivalent : les points au delà de 90°, en partant de l'équateur
        # n'ont pas d'existance sur le globe.
            
        # Correction de la précision des valeurs : (pour que le point -24.99 appartienne à [-24.9899997711, 24.9899997711]
        if (float("%e" % lonmin)<=lon) and (lon<=float("%e" % lonmax)) and (float("%e" % selflat[0])<=lat) and (lat<=float("%e" % selflat[-1])):
            return True
        else:
            return False
            
        
    def xy2latlon(self, x, y):
        """
        Renvoie la latitude correspondant à l'indice y
        et la longitude correspondant à l'indice x.
        
        Pas d'amélioration de précision ici. Une valeur comme -74.589996 signifie -74.59
        si la résolution n'est pas de l'ordre de 10^-5 degré.
        @return: (lat, lon)
        """
        lon = self.lonValues.values[x]
        lat = self.latValues.values[y]
        return (lat, lon)

    def latlon2xy(self, lat, lon):
        """
        Se base sur getNearestIndex pour renvoyer les indices
        (x pour lon, y pour lat) en un seul appel de fonction.

        Si la latitude ou longitude est hors grille, le point bordure est renvoyé
        @return: (x,y)"""
        x = self.getNearestIndex(lat,'lat')
        y = self.getNearestIndex(lon,'lon')
        return (x, y)

    def save(self, output, attrs=None):
        """
        Save the swath to a storage (file,...)

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
            output.create_dim('row', self.get_geolocation_dimsizes()['row'])
            output.create_dim('cell', self.get_geolocation_dimsizes()['cell'])
            dims = ['row', 'cell']
            # create additional dimensions
            for v in self._fields.keys():
                if v not in self._geolocation_fields:
                    for d in self._fields[v].dimensions:
                        if d not in dims:
                            output.create_dim(d, self._fields[v].get_dimsize(d))
                            dims.append(d)
            # creating metadata
            if attrs:
                globalattr = attrs
            else:
                globalattr = {}
            if not self._bbox:
                self._bbox = (
                            self.get_lon().min(),
                            self.get_lat().min(),
                            self.get_lon().max(),
                            self.get_lat().max()
                            )
            globalattr['geospatial_lat_min'] = self._bbox[1]
            globalattr['geospatial_lon_min'] = self._bbox[0]
            globalattr['geospatial_lat_max'] = self._bbox[3]
            globalattr['geospatial_lon_max'] = self._bbox[2]

            # Calculate the time bounds of the output file
            times = self.get_times()
            time0 = num2date(times.min(),
                             self.get_geolocation_field('time').units)
            time1 = num2date(times.max(),
                             self.get_geolocation_field('time').units)
            if 'time_coverage_start' not in globalattr:
                globalattr['time_coverage_start'] = \
                    min(time0, time1).strftime('%Y%m%dT%H%M%S')
            if 'time_coverage_end' not in globalattr:
                globalattr['time_coverage_end'] = \
                    max(time0, time1).strftime('%Y%m%dT%H%M%S')

            # Add the global attributes
            globalattr['cdm_data_type'] = 'swath'
            output.write_global_attributes(globalattr)

            # creating records
            for geof in self._geolocation_fields:
                if not self._geolocation_fields['time']:
                    raise Exception('No time information defined')
                if not self._geolocation_fields[geof]:
                    logging.warning('Missing geolocation variable : %s', geof)
                else:
                    output.create_field(self._geolocation_fields[geof])
            for dataf in self._fields.keys():
                if dataf not in self._geolocation_fields:
                    output.create_field(self._fields[dataf])
        # saving records
        for geof in self._geolocation_fields:
            if self._geolocation_fields[geof]:
                output.write_field(self._geolocation_fields[geof])
        for dataf in self._fields.keys():
            output.write_field(self._fields[dataf])
        output.sync()
