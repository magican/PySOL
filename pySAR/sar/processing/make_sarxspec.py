#!/usr/bin/env python
# coding=utf-8
"""
"""


from sar.utils.factory import sarimage
from sar.transform.sarimage2sarxspec import sarimage2sarxspec_loop
from netCDF4 import Dataset
import os
import matplotlib.pyplot as plt
import numpy as np
from sar.render.palette import getColorMap
from scipy.stats import scoreatpercentile
from datetime import datetime


# def make_roughness_png(sarim, sarxspec, pngname):
#     """
#     """
#     spacing_m = 40.
#     spacing = np.round(sarim.meters2pixels(spacing_m))
#     ssr = np.sqrt(sarim.get_data('roughness', spacing=spacing))
#     dim = ssr.shape
#     (vmin, vmax) = (np.inf, -np.inf)
#     for ibur in range(sarim.get_info('number_of_bursts')):
#         extb = sarim.extent_burst(ibur, valid=True)
#         ssrb = np.sqrt(sarim.get_data('roughness', extent=extb,
#                                       spacing=spacing))
#         vmin = np.minimum(vmin, scoreatpercentile(ssrb, 1))
#         vmax = np.maximum(vmax, scoreatpercentile(ssrb, 99))
#     # vmin = scoreatpercentile(ssr[dim[0]*.05:dim[0]*.95,
#     #                              dim[1]*.05:dim[1]*.95], 0.1)
#     # vmax = scoreatpercentile(ssr[dim[0]*.05:dim[0]*.95,
#     #                              dim[1]*.05:dim[1]*.95], 99.9)
#     dpi = 100.
#     figsize = (dim[1]/dpi, dim[0]/dpi)
#     fig = plt.figure(figsize=figsize, dpi=dpi)
#     fig.subplots_adjust(bottom=0, top=1, left=0, right=1, wspace=0,
#                         hspace=0)
#     ax = plt.subplot(1, 1, 1)
#     ax.set_xlim([0, dim[1]-1])
#     ax.set_ylim([0, dim[0]-1])
#     ax.imshow(ssr, origin='lower', cmap=plt.cm.Greys_r, vmin=vmin, vmax=vmax)
#     ax.set_axis_off()
#     for iazi in range(len(sarxspec)):
#         for iran in range(len(sarxspec[0])):
#             ext = np.array(sarxspec[iazi][iran].get_info('extent'),
#                            dtype='float32')
#             extx = ext[[1, 3, 3, 1, 1]]/spacing[1]
#             exty = ext[[0, 0, 2, 2, 0]]/spacing[0]
#             ax.plot(extx, exty, 'r')
#             ax.text((extx[0]+extx[1])/2, (exty[0]+exty[2])/2,
#                     str(iran)+'-'+str(iazi), ha='center', va='bottom',
#                     color='red')
#             nvar = "%.2f" % sarxspec[iazi][iran].get_info('normalized_variance')
#             ax.text((extx[0]+extx[1])/2, (exty[0]+exty[2])/2,
#                     'nv='+nvar, ha='center', va='top', color='red')
#     fig.savefig(pngname)
#     plt.close(fig)


def ax_txt_corner(ax, txt, pos, **kwargs):
    """
    """
    pos2set = {'tl': {'x':0, 'y':1, 'ha':'left', 'va':'top'},
               'tr': {'x':1, 'y':1, 'ha':'right', 'va':'top'},
               'bl': {'x':0, 'y':0, 'ha':'left', 'va':'bottom'},
               'br': {'x':1, 'y':0, 'ha':'right', 'va':'bottom'}}
    xlim = ax.get_xlim()
    x = pos2set[pos]['x']
    xpos = xlim[x] + 0.01*(xlim[1-x]-xlim[x])
    ylim = ax.get_ylim()
    y = pos2set[pos]['y']
    ypos = ylim[y] + 0.01*(ylim[1-y]-ylim[y])
    ha = pos2set[pos]['ha']
    va = pos2set[pos]['va']
    ax.text(xpos, ypos, txt, ha=ha, va=va, **kwargs)


def make_sarxspec_fig(sarxspec, tau=1, part='real',
                      kmax=2*np.pi/75, kmin=2*np.pi/400,
                      xspec_size=(256, 256),
                      uniq_vmax=True, vmax=None, north_oriented=False,
                      klim=[2*np.pi/400, 2*np.pi/200, 2*np.pi/100],
                      north_arrow=True, index_pos='tl', vmax_pos='tr',
                      nvar_pos='br', fontsize='x-small', cmap=None):
    """
    """
    # Constants
    naz, nra = (len(sarxspec), len(sarxspec[0]))
    heading00 = sarxspec[0][0].get_info('platform_heading')
    reverse = (north_oriented == True) and (np.cos(heading00*np.pi/180) < 0)
    if cmap is None:
        if part == 'real':
            cmap = getColorMap(rgbFile='wind.pal')
        elif part == 'imag':
            cmap = plt.cm.PuOr
    if vmax is not None:
        uniq_vmax = True # force uniq_vmax
    elif uniq_vmax == True:
        vmaxs, nvars = [], []
        for iaz in range(naz):
            for ira in range(nra):
                spec = sarxspec[iaz][ira].get_data(kmax=kmax, tau=tau, part=part)
                k = sarxspec[iaz][ira].get_k2d(kmax=kmax)
                vmaxs.append(abs(spec[k >= kmin]).max())
                nvars.append(sarxspec[iaz][ira].get_info('normalized_variance'))
        ind = np.where(np.array(nvars) < 2)
        if ind[0].size == 0:
            vmax = np.median(np.array(vmaxs))
        else:
            vmax = np.array(vmaxs)[ind].max()
    theta = 2*np.pi*np.linspace(0, 1, num=361)
    # Create figure
    dpi = 100.
    fig_size = np.array((nra, naz), dtype='float')*xspec_size
    fig_size_inch = fig_size/dpi
    fig = plt.figure(figsize=fig_size_inch, dpi=dpi, facecolor='w')
    # fig.subplots_adjust(left=0., right=1., wspace=0.,
    #                     bottom=0., top=1., hspace=0.)
    # Make figure
    for iaz in range(naz):
        for ira in range(nra):
            # Get data
            spec = sarxspec[iaz][ira].get_data(kmax=kmax, tau=tau, part=part)
            k = sarxspec[iaz][ira].get_k2d(kmax=kmax)
            (kaz, kra) = sarxspec[iaz][ira].get_k1d(kmax=kmax)
            heading = sarxspec[iaz][ira].get_info('platform_heading')
            north = (90+heading)*np.pi/180 # north dir in spectrum, rad
            if uniq_vmax == False:
                vmax = abs(spec[k >= kmin]).max()
            if part == 'real':
                vmin = 0.
            elif part == 'imag':
                vmin = -vmax
            if reverse == True:
                spec = spec[::-1, ::-1]
                kaz, kra = -kaz[::-1], -kra[::-1]
                north += np.pi
            # Create axes and plot
            # if reverse == True:
            #     ax = fig.add_subplot(naz, nra, (nra-1-ira)+1+iaz*nra)
            # else:
            #     ax = fig.add_subplot(naz, nra, ira+1+(naz-1-iaz)*nra)
            if reverse == True:
                left = (nra-1-ira)*xspec_size[0]/fig_size[0]
                bottom = (naz-1-iaz)*xspec_size[1]/fig_size[1]
            else:
                left = ira*xspec_size[0]/fig_size[0]
                bottom = iaz*xspec_size[1]/fig_size[1]
            rect = [left, bottom, xspec_size[0]/fig_size[0],
                    xspec_size[1]/fig_size[1]]
            ax = fig.add_axes(rect)
            axkmin = np.minimum(kaz[0], kra[0])
            axkmax = np.maximum(kaz[-1], kra[-1])
            ax.set_ylim([axkmin, axkmax])
            ax.set_xlim([axkmin, axkmax])
            imsh = ax.imshow(spec, cmap=cmap, vmin=vmin, vmax=vmax,
                             extent=[kra[0], kra[-1], kaz[0], kaz[-1]],
                             origin='lower', interpolation='nearest',
                             aspect='equal')
            ax.set_axis_off()
            ax.plot([axkmin, axkmax], [0, 0], ':k')
            # ax.text(axkmax, 0, 'Ra', ha='right', va='top',
            #         fontsize='x-small')
            ax.plot([0, 0], [axkmin, axkmax], ':k')
            # ax.text(0, axkmax, 'Az', ha='left', va='top',
            #         fontsize='x-small')
            if klim is not None:
                for kli in klim:
                    ax.plot(kli*np.cos(theta), kli*np.sin(theta), ':k')
                    klistr = '%im' % (np.round(2*np.pi/kli))
                    ax.text(0, -kli, klistr, ha='center', va='top',
                            fontsize=fontsize)
            if north_arrow == True:
                k400 = 2*np.pi/400
                ax.arrow(0, 0, k400*np.cos(north), k400*np.sin(north),
                         length_includes_head=True, width=0.0001, facecolor='k')
            if index_pos is not None:
                txt = str(ira)+'-'+str(iaz)
                ax_txt_corner(ax, txt, index_pos, fontsize=fontsize)
            if vmax_pos is not None:
                if uniq_vmax == False:
                    txt = 'max=%.2e' % vmax
                else:
                    txt = 'max=%.2e' % abs(spec[k >= kmin]).max()
                ax_txt_corner(ax, txt, vmax_pos, fontsize=fontsize)
            if nvar_pos is not None:
                nvar = sarxspec[iaz][ira].get_info('normalized_variance')
                txt = 'nv=%.2f' % nvar
                ax_txt_corner(ax, txt, nvar_pos, fontsize=fontsize)
    return fig


# def make_sarxspec_png(sarxspec, realpngname, imagpngname):
#     """
#     """
#     kmax = 2*np.pi/90
#     kmin = 2*np.pi/400
#     nazi = len(sarxspec)
#     nran = len(sarxspec[0])
#     dpi = 100.
#     # figsize = (280*nran/dpi, 225*nazi/dpi)
#     # figre = plt.figure(figsize=figsize, dpi=dpi)
#     # figre.subplots_adjust(top=1-25./225./nazi, hspace=25./(225-25./nazi),
#     #                       right=1-25./280/nran, wspace=25./(280-25./nran),
#     #                       bottom=0., left=0.)
#     # figim = plt.figure(figsize=figsize, dpi=dpi)
#     # figim.subplots_adjust(top=1-25./225./nazi, hspace=25./(225-25./nazi),
#     #                       right=1-25./280/nran, wspace=25./(280-25./nran),
#     #                       bottom=0., left=0.)
#     figsize = (255*nran/dpi, 225*nazi/dpi)
#     figre = plt.figure(figsize=figsize, dpi=dpi)
#     figre.subplots_adjust(top=1-25./225./nazi, hspace=25./(225-25./nazi),
#                           right=1., wspace=0.,
#                           bottom=0., left=0.)
#     figim = plt.figure(figsize=figsize, dpi=dpi)
#     figim.subplots_adjust(top=1-25./225./nazi, hspace=25./(225-25./nazi),
#                           right=1., wspace=0.,
#                           bottom=0., left=0.)
#     k100 = 2*np.pi/100
#     k200 = k100/2
#     k400 = k100/4
#     k600 = k100/6
#     theta = 2*np.pi*np.linspace(0, 1, num=361)
#     cmapre = getColorMap(rgbFile='wind.pal')
#     cmapim = plt.cm.PuOr
#     vmaxre = -np.inf
#     vmaxim = -np.inf
#     for iazi in range(nazi):
#         for iran in range(nran):
#             nvar = sarxspec[iazi][iran].get_info('normalized_variance')
#             if nvar < 2:
#                 spec = sarxspec[iazi][iran].get_data()[:, :, 1]
#                 (kazi, kran) = sarxspec[iazi][iran].get_k()
#                 k = np.sqrt(np.tile(kazi.reshape(-1, 1), (1, kran.size))**2 + \
#                             np.tile(kran.reshape(1, -1), (kazi.size, 1))**2)
#                 vmaxre = np.maximum(vmaxre, spec[abs(k) > kmin].real.max())
#                 vmaxim = np.maximum(vmaxim, abs(spec[abs(k) > kmin].imag).max())
#     for iazi in range(nazi):
#         for iran in range(nran):
#             # Handle wavenumbers
#             spec = sarxspec[iazi][iran].get_data()[:, :, 1]
#             (kazi, kran) = sarxspec[iazi][iran].get_k()
#             indkazi = (np.where(abs(kazi) <= kmax))[0]
#             spec = spec[indkazi, :]
#             kazi = kazi[indkazi]
#             indkran = (np.where(abs(kran) <= kmax))[0]
#             spec = spec[:, indkran]
#             kran = kran[indkran]
#             k = np.sqrt(np.tile(kazi.reshape(-1, 1), (1, kran.size))**2 + \
#                         np.tile(kran.reshape(1, -1), (kazi.size, 1))**2)
#             heading = sarxspec[iazi][iran].get_info('platform_heading')
#             north = (90+heading)*np.pi/180 # north dir in spectrum, rad
#             dkazi, dkran = (kazi[1]-kazi[0], kran[1]-kran[0])
#             #integ = "%.2e" % (np.sum(abs(spec[abs(k) > kmin]))*dkazi*dkran)
#             nvar = "%.2f" % sarxspec[iazi][iran].get_info('normalized_variance')
#             # Make figures
#             for ifig in range(2):
#                 if ifig == 0:
#                     fig = figre
#                     sp = spec.real
#                     # vmax = spec[abs(k) > kmin].real.max()
#                     # vmin = spec[abs(k) > kmin].real.min()
#                     vmax = vmaxre
#                     vmin = 0.
#                     cmap = cmapre
#                 elif ifig == 1:
#                     fig = figim
#                     sp = spec.imag
#                     # vmax = abs(spec[abs(k) > kmin].imag).max()
#                     # vmin = -vmax
#                     vmax = vmaxim
#                     vmin = -vmax
#                     cmap = cmapim
#                 ax = fig.add_subplot(nazi, nran, iran+1+(nazi-1-iazi)*nran)
#                 axkmin = np.minimum(kazi[0], kran[0])
#                 axkmax = np.maximum(kazi[-1], kran[-1])
#                 ax.set_ylim([axkmin, axkmax])
#                 ax.set_xlim([axkmin, axkmax])
#                 imsh = ax.imshow(sp, origin='lower', cmap=cmap,
#                                  vmin=vmin, vmax=vmax,
#                                  interpolation='nearest', aspect='equal',
#                                  extent=[kran[0], kran[-1], kazi[0], kazi[-1]])
#                 ax.set_axis_off()
#                 ax.plot(k600*np.cos(theta), k600*np.sin(theta), ':k')
#                 ax.plot(k400*np.cos(theta), k400*np.sin(theta), ':k')
#                 ax.text(0, -k400, '400m', ha='center', fontsize='x-small',
#                         va='top')
#                 ax.plot(k200*np.cos(theta), k200*np.sin(theta), ':k')
#                 ax.text(0, -k200, '200m', ha='center', fontsize='x-small',
#                         va='top')
#                 ax.plot(k100*np.cos(theta), k100*np.sin(theta), ':k')
#                 ax.text(0, -k100, '100m', ha='center', fontsize='x-small',
#                         va='bottom')
#                 ax.plot([axkmin, axkmax], [0, 0], ':k')
#                 ax.text(axkmax, 0, 'Ra', ha='right', va='top',
#                         fontsize='x-small')
#                 ax.plot([0, 0], [axkmin, axkmax], ':k')
#                 ax.text(0, axkmax, 'Az', ha='left', va='top',
#                         fontsize='x-small')
#                 ax.arrow(0, 0, k400*np.cos(north), k400*np.sin(north),
#                          length_includes_head=True, width=0.0001,
#                          facecolor='k')
#                 # ax.set_title(str(iran)+'-'+str(iazi)+' int='+integ,
#                 #              fontsize='small', color='blue')
#                 ax.set_title(str(iran)+'-'+str(iazi)+' nv='+nvar,
#                              fontsize='small', color='blue')
#                 #cb = fig.colorbar(imsh, shrink=0.5, format='%.1e')
#                 # for t in cb.ax.get_yticklabels():
#                 #     t.set_fontsize('x-small')
#                 cb = fig.colorbar(imsh, shrink=0.5, format='',
#                                   ticks=[vmin, vmax])
#     figre.savefig(realpngname)
#     plt.close(figre)
#     figim.savefig(imagpngname)
#     plt.close(figim)


def make_sarxspec_netcdf(sarxspec, ncname):
    """
    """
    # Get infos
    naz, nra = (len(sarxspec), len(sarxspec[0]))
    nkaz00 = sarxspec[0][0].get_info('azimuth_look_size')
    nkaz = nkaz00 + (nkaz00 % 2)
    kazmin = sarxspec[0][0].get_info('azimuth_size')/2 - nkaz/2
    kazmax = sarxspec[0][0].get_info('azimuth_size')/2 + nkaz/2
    nkra00 = sarxspec[0][0].get_info('range_look_size')
    nkra = nkra00 + (nkra00 % 2)
    kramin = sarxspec[0][0].get_info('range_size')/2 - nkra/2
    kramax = sarxspec[0][0].get_info('range_size')/2 + nkra/2
    # Create netcdf dataset and global attributes
    rootgrp = Dataset(ncname, mode='w', format='NETCDF3_CLASSIC')
    # Create dimensions
    oswRaSize = rootgrp.createDimension('oswRaSize', nra)
    oswAzSize = rootgrp.createDimension('oswAzSize', naz)
    oswRaKSize = rootgrp.createDimension('oswRaKSize', nkra)
    oswAzKSize = rootgrp.createDimension('oswAzKSize', nkaz)
    # Create variables and their attributes
    oswLon = rootgrp.createVariable('oswLon', 'f4', ('oswAzSize', 'oswRaSize'))
    oswLon.units = "degree (decimal)"
    oswLon.long_name = "longitude"
    oswLat = rootgrp.createVariable('oswLat', 'f4', ('oswAzSize', 'oswRaSize'))
    oswLat.units = "degree (decimal)"
    oswLat.long_name = "latitude"
    oswGroundRngSize = rootgrp.createVariable('oswGroundRngSize', 'f4', ('oswAzSize', 'oswRaSize'))
    oswGroundRngSize.units = "m"
    oswGroundRngSize.long_name = "ground range size of estimation area"
    oswAziSize = rootgrp.createVariable('oswAziSize', 'f4', ('oswAzSize', 'oswRaSize'))
    oswAziSize.units = "m"
    oswAziSize.long_name = "azimuth size of estimation area"
    oswIncidenceAngle = rootgrp.createVariable('oswIncidenceAngle', 'f4', ('oswAzSize', 'oswRaSize'))
    oswIncidenceAngle.units = "degree"
    oswIncidenceAngle.long_name = "incidence angle"
    oswHeading = rootgrp.createVariable('oswHeading', 'f4', ('oswAzSize', 'oswRaSize'))
    oswHeading.units = "degree"
    oswHeading.long_name = "platform heading"
    oswNv = rootgrp.createVariable('oswNv', 'f4', ('oswAzSize', 'oswRaSize'))
    oswNv.units = "adimentional"
    oswNv.long_name = "normalized variance"
    oswCrossSpectraRe = rootgrp.createVariable('oswCrossSpectraRe', 'f4', ('oswAzSize', 'oswRaSize', 'oswAzKSize', 'oswRaKSize'))
    oswCrossSpectraRe.units = "" # "(m/rad)^2"
    oswCrossSpectraRe.long_name = "real part of cross spectra cartesian grid"
    oswCrossSpectraIm = rootgrp.createVariable('oswCrossSpectraIm', 'f4', ('oswAzSize', 'oswRaSize', 'oswAzKSize', 'oswRaKSize'))
    oswCrossSpectraIm.units = "" # "(m/rad)^2"
    oswCrossSpectraIm.long_name = "imaginary part of cross spectra cartesian grid"
    oswRaK = rootgrp.createVariable('oswRaK', 'f4', ('oswRaSize', 'oswRaKSize'))
    oswRaK.units = "rad/m"
    oswRaK.long_name = "array of the range wavenumbers for cartesian spectra"
    oswAzK = rootgrp.createVariable('oswAzK', 'f4', ('oswAzKSize',))
    oswAzK.units = "rad/m"
    oswAzK.long_name = "array of the azimuth wavenumbers for cartesian spectra"
    # Write variables
    for iaz in range(naz):
        for ira in range(nra):
            if nkaz00 != sarxspec[iaz][ira].get_info('azimuth_look_size'):
                raise Exception('Azimuth size is not constant for all XSpecs')
            if nkra00 != sarxspec[iaz][ira].get_info('range_look_size'):
                raise Exception('Range size is not constant for all XSpecs')
            oswLon[iaz, ira] = sarxspec[iaz][ira].get_info('longitude')
            oswLat[iaz, ira] = sarxspec[iaz][ira].get_info('latitude')
            oswGroundRngSize[iaz, ira] = sarxspec[iaz][ira].get_info('range_dist')
            oswAziSize[iaz, ira] = sarxspec[iaz][ira].get_info('azimuth_dist')
            oswIncidenceAngle[iaz, ira] = sarxspec[iaz][ira].get_info('incidence')
            oswHeading[iaz, ira] = sarxspec[iaz][ira].get_info('platform_heading')
            oswNv[iaz, ira] = sarxspec[iaz][ira].get_info('normalized_variance')
            specs = sarxspec[iaz][ira].get_data(tau=1)
            oswCrossSpectraRe[iaz, ira, :, :] = specs[kazmin:kazmax, kramin:kramax].real
            oswCrossSpectraIm[iaz, ira, :, :] = specs[kazmin:kazmax, kramin:kramax].imag
    for ira in range(nra):
        (kaz, kra) = sarxspec[0][ira].get_k1d()
        oswRaK[ira, :] = kra[kramin:kramax]
    oswAzK[:] = kaz[kazmin:kazmax]
    rootgrp.close()


def make_sarxspec(sar_fname, output_dir, output_type, **kwargs):
    """
    """
    # Check output_type before any computations
    if isinstance(output_type, basestring):
        output_type = [output_type]
    for otype in output_type:
        if otype not in ['netcdf', 'plot']:
            raise Exception(str(otype)+' is not an output type')
    # Init SAR image
    sarim = sarimage(sar_fname)
    # Handle output filename pattern
    dtime = datetime.strptime(sarim.get_info('start_time'),
                             '%Y-%m-%dT%H:%M:%S.%f')
    year = str(dtime.year)
    doy = str(dtime.timetuple().tm_yday)
    safe = os.path.basename(os.path.dirname(os.path.dirname(sar_fname)))
    output_dname = os.path.join(output_dir, year, doy, safe)
    base = os.path.splitext(os.path.basename(sar_fname))[0]
    output_fname = os.path.join(output_dname, base+'-xspec')
    # Compute SAR cross-spectra
    sarxspec = sarimage2sarxspec_loop(sarim, **kwargs)
    # Outputs
    if os.path.exists(output_dname) == False:
        os.makedirs(output_dname)
    if 'netcdf' in output_type:
        make_sarxspec_netcdf(sarxspec, output_fname+'.nc')
    if 'plot' in output_type:
        mode = sarim.get_info('mode')
        if mode == 'WV':
            fig = sarxspec[0][0].display_data()
            fig.savefig(output_fname+'_realimag_1tau2tau.png')
            plt.close(fig)
        else:
            fig = make_sarxspec_fig(sarxspec, tau=1, part='real')
            fig.savefig(output_fname+'_real_1tau.png')
            plt.close(fig)
            fig = make_sarxspec_fig(sarxspec, tau=1, part='imag')
            fig.savefig(output_fname+'_imag_1tau.png')
            plt.close(fig)
    # if 'roughness-plot':
    #     make_roughness_png(sarim, sarxspec, output_name+'-roughness.png')


### S1 L2 OCN netcdf template ###
# [gilles@amanita l2_lop_standalone]$ ncdump -h s1a-s3-osw-vv-20140507t140856-20140507t140921-000492-005427-001.nc
# netcdf s1a-s3-osw-vv-20140507t140856-20140507t140921-000492-005427-001 {
# dimensions:
# 	oswRaSize = 4 ;
# 	oswAzSize = 10 ;
# 	oswWavenumberBinSize = 30 ;
# 	oswAngularBinSize = 36 ;
# 	oswPartitions = 2 ;
# variables:
# 	float oswLon(oswAzSize, oswRaSize) ;
# 		oswLon:units = "degree (decimal)" ;
# 		oswLon:long_name = "longitude" ;
# 	float oswLat(oswAzSize, oswRaSize) ;
# 		oswLat:units = "degree (decimal)" ;
# 		oswLat:long_name = "latitude" ;
# 	float oswGroundRngSize(oswAzSize, oswRaSize) ;
# 		oswGroundRngSize:units = "m" ;
# 		oswGroundRngSize:long_name = "ground range size of estimation area" ;
# 	float oswAziSize(oswAzSize, oswRaSize) ;
# 		oswAziSize:units = "m" ;
# 		oswAziSize:long_name = "azimuth size of estimation area" ;
# 	float oswIncidenceAngle(oswAzSize, oswRaSize) ;
# 		oswIncidenceAngle:units = "degree" ;
# 		oswIncidenceAngle:long_name = "incidence angle" ;
# 	float oswHeading(oswAzSize, oswRaSize) ;
# 		oswHeading:units = "degree" ;
# 		oswHeading:long_name = "platform heading" ;
# 	float oswInten(oswAzSize, oswRaSize) ;
# 		oswInten:units = "adimentional" ;
# 		oswInten:long_name = "image intensity" ;
# 	float oswNv(oswAzSize, oswRaSize) ;
# 		oswNv:units = "adimentional" ;
# 		oswNv:long_name = "normalized variance" ;
# 	float oswSkew(oswAzSize, oswRaSize) ;
# 		oswSkew:units = "adimentional" ;
# 		oswSkew:long_name = "image skewness" ;
# 	float oswKurt(oswAzSize, oswRaSize) ;
# 		oswKurt:units = "adimentional" ;
# 		oswKurt:long_name = "image kurtosis" ;
# 	float oswPolSpec(oswAzSize, oswRaSize, oswAngularBinSize, oswWavenumberBinSize) ;
# 		oswPolSpec:units = "m^4" ;
# 		oswPolSpec:long_name = "estimated swell spectrum on log-polar grid" ;
# 	float oswPartitions(oswAzSize, oswRaSize, oswAngularBinSize, oswWavenumberBinSize) ;
# 		oswPartitions:units = "adimentional" ;
# 		oswPartitions:long_name = "wave partition numbers" ;
# 	float oswQualityCrossSpectraRe(oswAzSize, oswRaSize, oswAngularBinSize, oswWavenumberBinSize) ;
# 		oswQualityCrossSpectraRe:units = "(m/rad)^2" ;
# 		oswQualityCrossSpectraRe:long_name = "real part of cross spectra log-polar grid" ;
# 	float oswQualityCrossSpectraIm(oswAzSize, oswRaSize, oswAngularBinSize, oswWavenumberBinSize) ;
# 		oswQualityCrossSpectraIm:units = "(m/rad)^2" ;
# 		oswQualityCrossSpectraIm:long_name = "imaginary part of cross spectra log-polar grid" ;
# 	float oswK(oswWavenumberBinSize) ;
# 		oswK:units = "rad/m" ;
# 		oswK:long_name = "array of the wavenumbers for polar spectra" ;
# 	float oswPhi(oswAngularBinSize) ;
# 		oswPhi:units = "deg" ;
# 		oswPhi:long_name = "array of angular values for polar spectra" ;
# 	float oswLookSeparationTime(oswAzSize, oswRaSize) ;
# 		oswLookSeparationTime:units = "s" ;
# 		oswLookSeparationTime:long_name = "inter look separation time" ;
# 	float oswNrcs(oswAzSize, oswRaSize) ;
# 		oswNrcs:units = "dB" ;
# 		oswNrcs:long_name = "normalized radar cross section" ;
# 	float oswWindSeaHs(oswAzSize, oswRaSize) ;
# 		oswWindSeaHs:units = "m" ;
# 		oswWindSeaHs:long_name = "wind sea significant waveheights" ;
# 	float oswWaveAge(oswAzSize, oswRaSize) ;
# 		oswWaveAge:units = "adimentional" ;
# 		oswWaveAge:long_name = "inverse wave age" ;
# 	float oswSnr(oswAzSize, oswRaSize) ;
# 		oswSnr:units = "dB" ;
# 		oswSnr:long_name = "signal to noise ratio of cross spectra" ;
# 	float oswAzCutoff(oswAzSize, oswRaSize) ;
# 		oswAzCutoff:units = "m" ;
# 		oswAzCutoff:long_name = "azimuth cut-off wavelength" ;
# 	float oswRaCutoff(oswAzSize, oswRaSize) ;
# 		oswRaCutoff:units = "m" ;
# 		oswRaCutoff:long_name = "range cut-off wavelength" ;
# 	float oswSpecRes(oswAzSize, oswRaSize, oswAngularBinSize) ;
# 		oswSpecRes:units = "m" ;
# 		oswSpecRes:long_name = "azimuth cut-off wavelength as function of wave direction" ;
# 	float oswHs(oswAzSize, oswRaSize, oswPartitions) ;
# 		oswHs:units = "m" ;
# 		oswHs:long_name = "significant wave height for each partition" ;
# 	float oswWl(oswAzSize, oswRaSize, oswPartitions) ;
# 		oswWl:units = "m" ;
# 		oswWl:long_name = "Dominant wave length for each partition" ;
# 	float oswDirmet(oswAzSize, oswRaSize, oswPartitions) ;
# 		oswDirmet:units = "degrees clockwise from north" ;
# 		oswDirmet:long_name = "dominant wave direction for each partition" ;
# 	float oswAmbiFac(oswAzSize, oswRaSize, oswPartitions) ;
# 		oswAmbiFac:units = "adimentional" ;
# 		oswAmbiFac:long_name = "ambiguity ratio" ;
# 	float oswIconf(oswAzSize, oswRaSize, oswPartitions) ;
# 		oswIconf:units = "adimentional" ;
# 		oswIconf:long_name = "confidence of inversion for partition" ;
# 	float oswWindSpeed(oswAzSize, oswRaSize) ;
# 		oswWindSpeed:units = "m/s" ;
# 		oswWindSpeed:long_name = "wind speed" ;
# 	float oswWindDirection(oswAzSize, oswRaSize) ;
# 		oswWindDirection:units = "degree clockwise from range" ;
# 		oswWindDirection:long_name = "wind direction" ;
# 	float oswEcmwfWindSpeed(oswAzSize, oswRaSize) ;
# 		oswEcmwfWindSpeed:units = "m/s" ;
# 		oswEcmwfWindSpeed:long_name = "wind speed (ecmwf)" ;
# 	float oswEcmwfWindDirection(oswAzSize, oswRaSize) ;
# 		oswEcmwfWindDirection:units = "degree clockwise from north" ;
# 		oswEcmwfWindDirection:long_name = "wind direction (ecmwf)(where the wind comes from)" ;
# 	float oswLandCoverage(oswAzSize, oswRaSize) ;
# 		oswLandCoverage:units = "%" ;
# 		oswLandCoverage:long_name = "percentage of land coverage within a cell" ;
# 	float oswLandFlag(oswAzSize, oswRaSize) ;
# 		oswLandFlag:units = "adimentional" ;
# 		oswLandFlag:long_name = "land flag : 1 if presence of more than 10% of land coverage, 0 instead" ;
# 	float oswDepth(oswAzSize, oswRaSize) ;
# 		oswDepth:units = "m" ;
# 		oswDepth:long_name = "water detph" ;

# // global attributes:
# 		:TITLE = "Sentinel-1 OSW Component" ;
# 		:MISSION_NAME = "S1A" ;
# 		:SOURCE_PRODUCT = "S3" ;
# 		:LEVEL_1_SOURCE_NAME = "S1A_S3_SLC__1SSV_20140507T140856_20140507T140921_000492_0005FD_FF58.SAFE" ;
# 		:OSW_PROCESSING_UTC_TIME = "11-JUN-2014 08:16:28.000000" ;
# 		:OSW_PROCESSING_CENTER = "S-1 IPF" ;
# 		:OSW_PROCESSING_SOFTWARE = "s-1 osw V1.0" ;
# 		:oswAlgorithmVersion = "s-1 osw V1.0" ;
# 		:SOURCE_ACQUISITION_UTC_TIME = "07-MAY-2014 14:08:56.420751" ;
# 		:polarisation = "VV" ;
# 		:STATEVECTOR_UTC = "07-MAY-2014 14:08:54.421000" ;
# 		:STATEVECTOR_POS = -2588251.18795935, -5203343.65308994, 4016739.57503524 ;
# 		:STATEVECTOR_VEL = -3417.88044915425, -3007.93417692265, -6081.70610634732 ;
# 		:STATEVECTOR_ACC = 0., 0., 0. ;
# 		:firstMeasurementTime = "2014-05-07T14:08:56Z" ;
# 		:lastMeasurementTime = "2014-05-07T14:09:21Z" ;
# 		:OSW_CELL_SIZE_RANGE = 349 ;
# }
