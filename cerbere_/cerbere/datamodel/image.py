# -*- coding: utf-8 -*-
"""
cerbere.datamodel.image
============================

Model class for the image feature

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

import logging
import collections
import copy

import numpy

from .abstractfeature import AbstractFeature
from ..geo.point import Point

__all__ = ['Image']


class Image(AbstractFeature):
    """
    Model class for the image feature, ie a two-dimensional grid along the satellite track.

    The main difference with a image is that there is only one time 
    value associated with the feature (instead of a time value per
    pixel in a image).
    """
    def __init__(self, identifier=None, title=None, 
                description=None, source=None, metadata=None,
                fields=None,
                longitudes=None, latitudes=None, depths=None, times=None,
                spatial_resolution=None):
        super(Image, self).__init__(
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
        Returns the geolocation dimension names defining a image
        """
        return ['time', 'row', 'cell']

    def get_geolocation_field_dimnames(self, fieldname):
        """
        Returns the dimensions of a geolocation field
        """
        if fieldname == 'depth':
            return ('depth',)
        elif fieldname == 'time':
            return ('time',)
        else:
            return ('row', 'cell',)

    def get_geolocation_dimsizes(self):
        """
        Returns the geolocation dimension sizes of the image object
        """
        row, cell = self.get_lat().shape
        return collections.OrderedDict([('row', row),
                                        ('cell', cell)
                                        ])

    def get_values(self, fieldname, slices=None, indices=None, cache=False):
        '''
         if cache is True, data read from storage are kept in memory
         '''
        field = self._fields[fieldname]
        # remove time dimension
        if 'time' in field.dimensions:
            if slices is None:
                slices = {}
            slices['time'] = 0
        res = field.get_values(slices, indices, cache)
        return res

    def extract_subset(self, boundaries=None, slices=None, fields=None):
        """
        extract a subset from the current feature. The created subset is a new
        object without any reference to the source.

        Parameters
        ----------
        boundaries : tuple
            area of the subset to extract, defined as llcrnrlon, llcrnrlat,
            urcrnrlon, urcrnrlat

        slices : list of slices
            indices for x and y dimensions (or lon, lat) of the subset to
            extract from the source data

        fields (list)
            list of field names to extract
        """
        if boundaries and slices:
            raise Exception("Boundaries and slices can not be both provided.")
        if boundaries:
            # get corresponding slices
            raise Exception("Not yet implemented")
        elif slices:
            if len(slices) != 2:
                raise Exception ("slices for the two spatial dimensions of the image must be provided")
        subimage = Image(latitudes=self.extract_field('lat', slices=slices),
                        longitudes=self.extract_field('lon', slices=slices),
                        times=self.extract_field('time', slices=slices),
                        metadata=copy.copy(self.metadata)
                        )
        if fields is None:
            fields = self.get_fieldnames()
        elif not type(fields) is list:
            raise Exception("fields must be a list")
        for fieldname in fields:
            logging.debug("extract %s", fieldname)
            subimage.add_field(self.extract_field(fieldname, slices=slices))
        return subimage

    def get_spatial_resolution(self):
        """Return the spatial resolution of the feature"""
        if self.resolution is None:
            i = self.get_geolocation_dimsizes()['cell'] / 2
            j = self.get_geolocation_dimsizes()['row'] / 2
            dist = Point.get_distance(
                Point(
                    self.get_lon(
                        {'row': slice(j, j + 1, 1)},
                        {'cell': slice(i, i + 1, 1)}
                    ),
                    self.get_lat(
                        {'row': slice(j, j + 1, 1)},
                        {'cell': slice(i, i + 1, 1)}
                    )
                ),
                Point(
                    self.get_lon(
                        {'row': slice(j + 1, j + 2, 1)},
                        {'cell': slice(i + 1, i + 2, 1)}
                    ),
                    self.get_lat(
                        {'row': slice(j + 1, j + 2, 1)},
                        {'cell': slice(i + 1, i + 2, 1)}
                    )
                )
            )
            res = dist / 2.
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


    def save(self, output, attrs=None ):
        """
        Save the image to a storage (file,...)

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
            output.create_dim('time', None)
            output.create_dim('row', self.get_geolocation_dimsizes()['row'])
            output.create_dim('cell', self.get_geolocation_dimsizes()['cell'])
            dims = ['row', 'cell']
            for v in self._fields.keys():
                if not v in self._geolocation_fields:
                    for d in self._fields[v].dimensions:
                        if not d in dims:
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
                            self.get_lat().max(),
                            )
            
            globalattr['geospatial_lon_min'] = self._bbox[0]
            globalattr['geospatial_lat_min'] = self._bbox[1]
            globalattr['geospatial_lon_max'] = self._bbox[2]
            globalattr['geospatial_lat_max'] = self._bbox[3]
            #globalattr['geospatial_lat_resolution'] = self.get_lat(slice('x':x/2)) - self.get_lat()[0]
            #globalattr['geospatial_lon_resolution'] = self.get_lon()[1] - self.get_lon()[0]
            time0 = self.get_datetimes()[0]
            time1 = self.get_datetimes()[0]
            if not 'time_coverage_start' in globalattr:
                globalattr['time_coverage_start'] = \
                        min(time0, time1).strftime('%Y%m%dT%H%M%S')
            if not 'time_coverage_end' in globalattr:
                globalattr['time_coverage_end'] = \
                        max(time0, time1).strftime('%Y%m%dT%H%M%S')
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
                if not dataf in self._geolocation_fields:
                    output.create_field(self._fields[dataf])
        # saving records
        for geof in self._geolocation_fields:
            if self._geolocation_fields[geof]:
                output.write_field(self._geolocation_fields[geof])
        for dataf in self._fields.keys():
            output.write_field(self._fields[dataf])
        return
