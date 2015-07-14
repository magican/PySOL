#!/usr/bin/env python
# coding=utf-8
"""
"""


import numpy as np
from numpy.fft import fftshift, fftfreq
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sar.render.palette import getColorMap


class SARXSpec(object):
    """
    """

    def __init__(self, spectra, infos):
        """
        """
        self.spectra = spectra
        self.infos = infos
        return

    def get_infos(self):
        """
        """
        return self.infos

    def print_infos(self):
        """
        """
        for key, value in self.get_infos().items():
            if isinstance(value, dict):
                value2print = 'not printed dict'
            else:
                value2print = value
            print "%30s : %30s %s" % (key, value2print, type(value))

    def get_info(self, name):
        """
        """
        return self.infos[name]

    def get_k(self):
        """
        """
        azisz = self.get_info('azimuth_size')
        azidk = self.get_info('azimuth_dk')
        ransz = self.get_info('range_size')
        randk = self.get_info('range_dk')
        k = (fftshift(fftfreq(azisz))*azisz*azidk,
             fftshift(fftfreq(ransz))*ransz*randk)
        return k

    def get_k1d(self, kmax=None):
        """
        """
        azisz = self.get_info('azimuth_size')
        azidk = self.get_info('azimuth_dk')
        kazi = fftshift(fftfreq(azisz))*azisz*azidk
        ransz = self.get_info('range_size')
        randk = self.get_info('range_dk')
        kran = fftshift(fftfreq(ransz))*ransz*randk
        if kmax is not None:
            kazi = kazi[abs(kazi) <= kmax]
            kran = kran[abs(kran) <= kmax]
        k1d = (kazi, kran)
        return k1d

    def get_k2d(self, kmax=None):
        """
        """
        (kazi, kran) = self.get_k1d(kmax=kmax)
        k2d = np.sqrt(np.tile(kazi.reshape(-1, 1), (1, kran.size))**2 + \
                      np.tile(kran.reshape(1, -1), (kazi.size, 1))**2)
        return k2d

    def get_data(self, kmax=None, tau=None, part=None):
        """
        """
        index = []
        # Filter k
        if kmax is None:
            index.append(slice(0, self.get_info('azimuth_size')))
            index.append(slice(0, self.get_info('range_size')))
        else:
            (kazi, kran) = self.get_k1d()
            ind = np.where(abs(kazi) <= kmax)
            index.append(slice(ind[0].min(), ind[0].max()+1))
            ind = np.where(abs(kran) <= kmax)
            index.append(slice(ind[0].min(), ind[0].max()+1))
        # Filter tau
        if tau is None:
            nlooks = self.get_info('nlooks')
            index.append(slice(0, nlooks))
        else:
            index.append(tau)
        # Get spectra
        spec = np.array(self.spectra[index], copy=True)
        # Filter part and return
        if part == 'real':
            return spec.real
        elif part == 'imag':
            return spec.imag
        else:
            return spec

    def display_data(self, kmax=2*np.pi/50, kmin=2*np.pi/400):
        """
        """
        # Get data
        specs = self.get_data(kmax=kmax)
        (kaz, kra) = self.get_k1d(kmax=kmax)
        k = self.get_k2d(kmax=kmax)
        klim = [2*np.pi/400, 2*np.pi/200, 2*np.pi/100]
        theta = 2*np.pi*np.linspace(0, 1, num=361)
        heading = self.get_info('platform_heading')
        north = (90+heading)*np.pi/180 # north dir in spectrum, rad
        ntau = specs.shape[2]
        # Make figure
        dpi = 100
        figsize = (9, 4 * (ntau - 1))
        fig = plt.figure(figsize=figsize, dpi=dpi)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.95,
                            wspace=0.15, hspace=0.15)
        for ics in range(1, ntau):
            for ipl in range(2):
                if ipl == 0:
                    var = specs[:, :, ics].real
                    vmax = abs(var[k >= kmin]).max()
                    vmin = 0
                    tit = r'%i-$\tau$ Cross Spectrum Re' % ics
                    cmap = getColorMap(rgbFile='wind.pal')
                elif ipl == 1:
                    var = specs[:, :, ics].imag
                    vmax = abs(var[k >= kmin]).max()
                    vmin = -vmax
                    tit = r'%i-$\tau$ Cross Spectrum Im' % ics
                    cmap = cm.PuOr
                ax = plt.subplot(ntau - 1, 2, 2 * ics + ipl - 1)
                axkmin = np.minimum(kaz[0], kra[0])
                axkmax = np.maximum(kaz[-1], kra[-1])
                ax.set_ylim([axkmin, axkmax])
                ax.set_xlim([axkmin, axkmax])
                plt.imshow(var, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax,
                           interpolation='nearest', aspect='equal',
                           extent=[kra[0], kra[-1], kaz[0], kaz[-1]])
                ax.set_title(tit)
                ax.set_axis_off()
                for kli in klim:
                    ax.plot(kli*np.cos(theta), kli*np.sin(theta), ':k')
                    klistr = '%im' % (np.round(2*np.pi/kli))
                    ax.text(0, -kli, klistr, ha='center', fontsize='small',
                            va='top')
                ax.plot([axkmin, axkmax], [0, 0], ':k')
                ax.text(axkmax, 0, 'Range', ha='right', va='bottom',
                        fontsize='small')
                ax.plot([0, 0], [axkmin, axkmax], ':k')
                ax.text(0, axkmax, 'Azimuth', ha='left', va='top',
                        fontsize='small')
                k400 = 2*np.pi/400
                ax.arrow(0, 0, k400*np.cos(north), k400*np.sin(north),
                         length_includes_head=True, width=0.0001,
                         facecolor='k')
                if np.sin(north) > 0:
                    va = 'bottom'
                else:
                    va = 'top'
                ax.text(k400*np.cos(north), k400*np.sin(north), 'N',
                        ha='center', va=va, fontsize='small')
                plt.colorbar(shrink=0.66, format='%.1e')
        return fig
