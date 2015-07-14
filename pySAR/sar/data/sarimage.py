#!/usr/bin/env python
# coding=utf-8
"""SAR image module

This module provides a SAR image class built on top of a cerbere (SAR) mapper.

Notes
-----
In order to have horizontal range and vertical azimuth for printing and
displaying, azimuth/range are the first/second dimensions in inputs and
outputs. This is opposite to IDL habits.

Examples
--------
Getting sea surface roughness from a Sentinel-1A file at 500m spacing.

>>> from sar.utils.factory import sarmapper
>>> from sar.data.sarimage import SARImage
>>> sarmp = sarmapper('s1a-*.tiff')
>>> sarim = SARImage(sarmp)
>>> spacing = sarim.meters2pixels(500)
>>> roughness = sarim.get_data('roughness', spacing=spacing)
"""


import numpy as np
from sar.utils.gmf import compute_roughness
#import pdb


DEFAULTSTYLE = {'lon':'point',
                'lat':'point',
                'incidence':'point',
                'elevation':'point',
                'digital_number':'mean',
                'gamma0_lut':'point',
                'beta0_lut':'point',
                'sigma0_lut':'point',
                'gamma0':'mean',
                'beta0':'mean',
                'sigma0':'mean',
                'complex':'mean',
                'azimuth_time':'point',
                'slant_range_time':'point',
                'doppler_centroid':'point',
                'roughness':'mean'}


class SARImage(object):
    """SAR image class.

    Attributes
    ----------
    _mapper : cerbere mapper instance
    """

    def __init__(self, mapper):
        """Init.

        Parameters
        ----------
        mapper : cerbere mapper instance
        """
        self._mapper = mapper
        return

    def get_infos(self):
        """Get all SAR informations as a ``dict``."""
        return self._mapper.read_global_attributes()

    def print_infos(self):
        """Print formatted SAR informations."""
        for key, value in self.get_infos().items():
            if isinstance(value, dict):
                value2print = 'not printed dict'
            else:
                value2print = value
            print "%30s : %30s %s" % (key, value2print, type(value))

    def get_info(self, name):
        """Get one SAR information given a key name."""
        return self.get_infos()[name]

    def get_datanames(self):
        """
        """
        datanames = self._mapper.get_fieldnames()
        if 'time' in datanames:
            datanames.remove('time')
        return datanames

    def get_data(self, name, extent=None, spacing=[1, 1], midazimuth=False,
                 midrange=False, style=None, blocksize=50000000):
        """Get numerical values of a SAR grid variable.

        Parameters
        ----------
        name : str
            Variable name
        extent : array_like, optional
            [azimuth_start, range_start, azimuth_stop, range_stop] pixel extent.
            Default to [0, 0, azimuth_size-1, range_size-1]. Note that the
            limits defined by extent are all included in output unlike python's
            slice notation.
        spacing : array_like, optional
            [azimuth, range] pixel spacing. Default to [1, 1]. If only one value
            is provided, spacing is assumed equal for azimuth and range.
        midazimuth : bool, optionnal
            Default to False.
        midrange : bool, optionnal
            Default to False.
        style : str, optionnal
            Setting style='point' means data are read/computed every spacing
            pixel. Setting style='mean' means data are read/computed at full
            resolution then averaged every spacing pixel.
            Note that when style='point', the values are approximately at the
            center of the bins used for averaging when style='mean'.
            By default, style is automatically guessed from the given variable
            name (see global DEFAULTSTYLE).
        blocksize : int, optionnal

        Returns
        -------
        values : ndarray
            [azimuth, range] values.
        """
        # Manage inputs and make slices for mapper
        if extent is None:
            extent = self.extent_max()
        (ext, spa) = self._format_extent_spacing(extent, spacing,
                                                 midazimuth=midazimuth,
                                                 midrange=midrange)
        if style is None:
            style = DEFAULTSTYLE[name]
        if style == 'point':
            start = ext[0:2] + spa//2
            stop = ext[2:4] - (spa-1) + spa//2 + 1
            step = spa
        elif style == 'mean':
            start = ext[0:2]
            stop = ext[2:4] + 1
            step = (1, 1)
        else:
            raise Exception('Unknown style :'+str(style))
        slices = [slice(int(start[0]), int(stop[0]), int(step[0])),
                  slice(int(start[1]), int(stop[1]), int(step[1]))]
        # Get data by blocks for memory issues
        # for style='point', blocksize is handled by mapper
        if blocksize is not None and style == 'mean' and spa.max() > 1:
            azsize = slices[0].stop - slices[0].start
            rasize = slices[1].stop - slices[1].start
            if azsize*rasize > blocksize:
                blrasize = rasize
                blazsize = np.maximum(blocksize//blrasize, spa[0])
                blazsize -= blazsize % spa[0]
                nblocks = np.ceil(azsize/float(blazsize))
                blext = np.array(ext, copy=True)
                for ibl in np.arange(nblocks):
                    blext[0] = ext[0]+ibl*blazsize
                    blext[2] = np.minimum(blext[0]+blazsize-1, ext[2])
                    blvalues = self.get_data(name, extent=blext, spacing=spa,
                                             midazimuth=False, midrange=False,
                                             style=style, blocksize=None)
                    if ibl == 0:
                        dtype = blvalues.dtype
                        values = np.empty((azsize, rasize)/spa, dtype=dtype)
                    az0 = ibl*blazsize/spa[0]
                    values[az0:az0+blvalues.shape[0], :] = blvalues
                return values
        # Read values
        if name == 'roughness':
            values = self._mapper.read_values('sigma0', slices=slices,
                                              blocksize=blocksize)
        else:
            values = self._mapper.read_values(name, slices=slices,
                                              blocksize=blocksize)
        # Rebin
        if style == 'mean' and spa.max() > 1:
            sha = (values.shape[0]/spa[0], spa[0],
                   values.shape[1]/spa[1], spa[1])
            values = values.reshape(sha).mean(-1).mean(1)
        # Roughness special case
        if name == 'roughness':
            inc = self.get_data('incidence', spacing=spacing, extent=extent,
                                midazimuth=midazimuth, midrange=midrange,
                                blocksize=blocksize)
            pol = self.get_info('polarisation')
            values = compute_roughness(values, inc, pol)
        return values

    def _format_extent_spacing(self, extent, spacing, midazimuth=False,
                               midrange=False):
        """Format (and check) extent and spacing."""
        # Check extent
        ext = np.round(extent).flatten()
        if ext.size != 4:
            raise Exception('extent must contain 4 elements')
        extmax = self.extent_max()
        if (ext[0:2] < extmax[0:2]).any() or (ext[2:4] > extmax[2:4]).any():
            exttmp = np.array(ext)
            ext[0:2] = np.maximum(ext[0:2], extmax[0:2])
            ext[2:4] = np.minimum(ext[2:4], extmax[2:4])
            print 'Warning : extent is outside SAR image, '+str(exttmp)+\
                ' becomes '+str(ext)
        if (ext[0:2] > ext[2:4]).any():
            raise Exception('extent[0:2] must be less or equal than '+\
                            'extent[2:4]')
        # Check spacing
        spa = np.round(spacing).flatten()
        if spa.size == 1:
            spa = np.repeat(spa[0], 2)
        elif spa.size == 2:
            pass
        else:
            raise Exception('spacing must contain 1 or 2 elements')
        if (spa < [1, 1]).any():
            spatmp = np.array(spa)
            spa = np.maximum(spa, [1, 1])
            print 'Warning : spacing too small, '+str(spatmp)+' becomes '+\
                str(spa)
        if (spa > ext[2:4]-ext[0:2]+1).any():
            spatmp = np.array(spa)
            spa = np.minimum(spa, ext[2:4]-ext[0:2]+1)
            print 'Warning : spacing too large, '+str(spatmp)+' becomes '+\
                str(spa)
        # Make extent to be spacing modulo
        ext[2:4] -= (ext[2:4]-ext[0:2]+1) % spa
        # 1D extent
        if midazimuth == True:
            dim = (ext[2]-ext[0]+1)/spa[0]
            ext[0:3:2] = ext[0] + (dim-1)//2*spa[0] + [0, spa[0]-1]
        if midrange == True:
            dim = (ext[3]-ext[1]+1)/spa[1]
            ext[1:4:2] = ext[1] + (dim-1)//2*spa[1] + [0, spa[1]-1]
        return (ext, spa)

    def extent_max(self):
        """Get extent for the whole SAR image."""
        return np.array((0, 0, self.get_info('number_of_lines')-1,
                         self.get_info('number_of_samples')-1))

    def extent_burst(self, burst, valid=True):
        """Get extent for a SAR image burst."""
        nbursts = self.get_info('number_of_bursts')
        if nbursts == 0:
            raise Exception('No bursts in SAR image')
        if burst < 0 or burst >= nbursts:
            raise Exception('Invalid burst index number')
        if valid == True:
            burst_list = self.get_info('burst_list')
            extent = np.copy(burst_list['valid_location'][burst, :])
        else:
            extent = self.extent_max()
            nlines = self.get_info('lines_per_burst')
            extent[0:3:2] = [nlines*burst, nlines*(burst+1)-1]
        return extent

    def extent_from_lonlat(self, lonlat, dist):
        """Get extent from a geographic point and a distance."""
        raise Exception('Not yet implemented !')

    def ground_spacing(self, range_index=None, extent=None):
        """Get SAR image ground spacing.

        Parameters
        ----------
        range_index : int, optionnal
            Range index (pixel) at which the range ground spacing is computed
            for SLC products case. By default, corresponds to middle swath.
        extent : array_like, optionnal
            If provided, range_index is computed from the middle range extent.

        Returns
        -------
        ground_spacing : ndarray
            [azimuth, range] ground spacing in meters.

        Notes
        -----
        For GRD products, range_index and extent are ignored.
        """
        ground_spacing = np.array((self.get_info('azimuth_pixel_spacing'),
                                   self.get_info('range_pixel_spacing')))
        if self.get_info('product') == 'SLC':
            ext = self.extent_max()
            if extent is not None:
                ext[1:4:2] = extent[1:4:2]
            else:
                if range_index is not None:
                    ext[1:4:2] = range_index
            inc = self.get_data('incidence', extent=ext, midazimuth=True,
                                midrange=True)
            ground_spacing[1] /= np.sin(inc*np.pi/180)
        return ground_spacing

    def pixels2meters(self, pixels, **kwargs):
        """Convert pixels to meters."""
        meters = pixels*self.ground_spacing(**kwargs)
        return meters

    def meters2pixels(self, meters, **kwargs):
        """Convert meters to pixels."""
        pixels = meters/self.ground_spacing(**kwargs)
        return pixels
