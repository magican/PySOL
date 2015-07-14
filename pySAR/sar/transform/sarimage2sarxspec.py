#!/usr/bin/env python
# coding=utf-8
"""
"""


from sar.data.sarxspec import SARXSpec
from sar.render.palette import getColorMap
import numpy as np
from numpy.fft import fft2, fftshift, ifft2
from scipy.constants import c as clight
import matplotlib.pyplot as plt
from collections import OrderedDict as OD
import pdb


# def lowpass(inp, width, zeropadfac=2.5, rebinfac=0.3):
#     """Gaussian low pass filtering.
#     """
#     # see assert
#     # Handle dimensions
#     ndim = inp.ndim
#     inpshape = np.array(inp.shape, dtype='int32')
#     width = np.array(width, dtype='float32')
#     zeropadwidth = zeropadfac*width
#     rebinwidth = np.maximum(np.floor(rebinfac*width), 1.)
#     rinpshape = np.ceil(inpshape/rebinwidth).astype('int32')
#     zrinpshape = 2**np.ceil(np.log(rinpshape+zeropadwidth/rebinwidth)/np.log(2)).astype('int32')
#     prerinpshape = (rinpshape*rebinwidth).astype('int32')
#     # Make gaussian
#     gauss = []
#     for idim in np.arange(ndim):
#         k = fftfreq(zrinpshape[idim])*2*np.pi
#         d = (rebinwidth[idim]-1)*0.5/rebinwidth[idim]
#         g = np.exp(-0.5*(width[idim]/rebinwidth[idim]*k)**2 - 1j*k*d)
#         gauss.append(g)
#     # Apply low filter in fourier domain
#     if ndim == 1:
#         print 'TODO'
#     elif ndim == 2:
#         rpadwid = ((0, prerinpshape[0]-inpshape[0]),
#                    (0, prerinpshape[1]-inpshape[1]))
#         rshape = (rinpshape[0], rebinwidth[0], rinpshape[1], rebinwidth[1])
#         rinp = np.pad(inp, rpadwid, 'constant').reshape(rshape).mean(1).mean(-1)
#         zrpadwid = ((0, zrinpshape[0]-rinpshape[0]),
#                     (0, zrinpshape[1]-rinpshape[1]))
#         zrinp = np.pad(rinp, zrpadwid, 'constant')
#         lowinp = ifft2(fft2(zrinp)*np.outer(gauss[0], gauss[1])).real
#     return None


def estimate_dopplershift(profile):
    """Estimate Doppler shift from azimuth profile.
    """
    size = profile.size
    exp = np.exp(0 + 1j*2*np.pi*np.arange(size)/size)
    #ang = np.arctan(np.sum(profile*exp))
    pesum = np.sum(profile*exp)
    ang = np.arctan2(pesum.imag, pesum.real)
    if ang < 0:
        ang += 2*np.pi
    return np.round(ang/2/np.pi*size - size/2).astype('int32')


def sarimage2sarxspec(sarim, extent,
                      azi_periodo_size=1024, ran_periodo_size=None,
                      azi_spec_size=None, ran_spec_size=None,
                      azi_look_width=0.25, nlooks=3, look_sep=0.27,
                      ran_look_width=0.78, weighted=True, debug=None):
    """Transform a SARImage object into a SARXSpec object.
    """
    ###########################################################################
    # Get SAR inputs
    ###########################################################################
    image = sarim.get_data('complex', extent=extent)
    imshape = np.array(image.shape, dtype='int32')
    grdspacing = sarim.ground_spacing(extent=extent)
    prf = sarim.get_info('pulse_repetition_frequency')
    ###########################################################################
    # Deramp signal (TOPSAR mode only)
    ###########################################################################
    if sarim.get_info('number_of_bursts') != 0:
        # Get burst index
        nlb = sarim.get_info('lines_per_burst')
        bur0 = np.floor(extent[0]*1./nlb)
        bur1 = np.floor(extent[2]*1./nlb)
        if bur0 != bur1:
            raise Exception('extent covers more than one burst.')
        bur = bur0
        # Get inputs
        azitimeinterval = sarim.get_info('azimuth_time_interval')
        azisteeringrate = sarim.get_info('azimuth_steering_rate')*np.pi/180.
        vplat = sarim.get_info('platform_velocity')
        vgrnd = sarim.get_info('ground_velocity')
        radarfreq = sarim.get_info('radar_frequency')
        lambdaradar = clight/radarfreq
        azitime = (np.arange(nlb)-(nlb-1)/2.)*azitimeinterval
        azitime = azitime[extent[0]-bur*nlb:extent[2]+1-bur*nlb]
        rantime = sarim.get_data('slant_range_time', extent=extent,
                                 midazimuth=True)[0, :]
        afr_list = sarim.get_info('azimuthfmrate_list')
        if afr_list['nlines'] != sarim.get_info('number_of_bursts'):
            print 'WARNING : azimuth FM rate list size != number of bursts.'
        afr_t0 = afr_list['t0'][bur]
        afr_c0 = afr_list['c0'][bur]
        afr_c1 = afr_list['c1'][bur]
        afr_c2 = afr_list['c2'][bur]
        # Compute and apply deramping function
        # IDL Fab
        nu = 2.*np.sqrt(vplat*vgrnd)/clight
        w0 = 2.*np.pi*radarfreq
        beta = (nu**2*w0)/rantime
        gamma = 4*np.pi*vplat*azisteeringrate/lambdaradar
        dr = beta*gamma/(beta+gamma)#*1.25
        deramping = np.exp(1j*(-0.5)*np.outer(azitime**2, dr))
        # /IDL Fab
        # S1 DOC [S1-AD-19]-S1-TN-MDA-52-7445_L1_DPM.pdf page 84
        ka = afr_c0 + afr_c1*(rantime-afr_t0) + \
             afr_c2*(rantime-afr_t0)**2 # azimuth (doppler) fm rate
        ks = -2.*vplat/lambdaradar*azisteeringrate # doppler centroid rate
        #kt = -ka*ks/(ks-ka) # target instantaneous bandwidth rate
        ### gamma=-2*np.pi*ks and beta=-2*np.pi*ka ?
        kt = -ka*ks/(ka+ks) # = beta*gamma/(beta+gamma)/2./np.pi
        deramping = np.exp(-1j*np.pi*np.outer(azitime**2, kt))
        # /S1 DOC
        image *= deramping
        if debug == 'deramping':
            fig = plt.figure()
            plt.imshow(abs(image/deramping), origin='lower', cmap=plt.cm.Greys_r)
            fig = plt.figure()
            plt.imshow(abs(image), origin='lower', cmap=plt.cm.Greys_r)
            plt.show()
            pdb.set_trace()
    ###########################################################################
    # Clip bright targets
    ###########################################################################
    nstd = 15.
    intimage = abs(image)**2
    #print (intimage.std()/intimage.mean())**2
    intimage /= intimage.mean()
    clipvalue = nstd*intimage.std()/intimage.mean()
    ind = np.where(intimage > clipvalue)
    if len(ind[0]) != 0:
        fac = np.sqrt(intimage[ind]/clipvalue)
        image[ind] /= fac
        intimage = abs(image)**2
    nvar = (intimage.std()/intimage.mean())**2
    del intimage
    #print nvar
    ###########################################################################
    # Set periodograms/looks/specs sizes and positions
    ###########################################################################
    lookwidth = np.array((azi_look_width, ran_look_width))
    resamppow = np.ceil(np.log(lookwidth)/np.log(2)+1).astype('int32')
    # Change from LOP
    # avoid to make lookshape less than half specshape
    resamppow -= 1
    # /Change from LOP
    # Set specshape from pershape
    if azi_periodo_size is not None or ran_periodo_size is not None:
        ratio = grdspacing[0]/grdspacing[1]
        if azi_periodo_size is None:
            azi_periodo_size = 2**np.round(np.log(ran_periodo_size/ratio)/np.log(2))
        elif ran_periodo_size is None:
            ran_periodo_size = 2**np.round(np.log(azi_periodo_size*ratio)/np.log(2))
        pershape = np.array((azi_periodo_size, ran_periodo_size)).astype('int32')
        specshape = np.floor(2.**(resamppow)*pershape+0.5).astype('int32')
    # Set pershape from specshape
    elif azi_spec_size is not None or ran_spec_size is not None:
        ratio = grdspacing[0]/grdspacing[1]*(2.**(resamppow[1]-resamppow[0]))
        if azi_spec_size is None:
            azi_spec_size = 2**np.round(np.log(azi_spec_size/ratio)/np.log(2))
        elif ran_spec_size is None:
            ran_spec_size = 2**np.round(np.log(azi_spec_size*ratio)/np.log(2))
        specshape = np.array((azi_spec_size, ran_spec_size), dtype='int32')
        pershape = np.floor(2.**(-resamppow)*specshape+0.5).astype('int32')
    else:
        raise Exception('No sizes at all given.')
    if (pershape > imshape).any() == True:
        raise Exception('Periodogramm is larger than image !')
    # Change from LOP
    #npers = (imshape+pershape//2)//(pershape//2)-1 # LOP style
    npers = np.round((imshape-pershape)/(pershape/2.)+1).astype('int32')
    # /Change from LOP
    perpos = (np.floor(np.linspace(0, imshape[0]-pershape[0], num=npers[0])+
                       0.5).astype('int32'),
              np.floor(np.linspace(0, imshape[1]-pershape[1], num=npers[1])+
                       0.5).astype('int32'))
    lookshape = np.floor(pershape*lookwidth+0.5).astype('int32')
    lookpos = (np.floor((np.linspace(1-nlooks, nlooks-1, num=nlooks)*look_sep+
                         1-azi_look_width)*.5*pershape[0]).astype('int32'),
               np.floor(0.5*(pershape[1]-lookshape[1])+0.5).astype('int32'))
    ###########################################################################
    # Get doppler shift
    ###########################################################################
    dopcentroid = sarim.get_data('doppler_centroid', extent=extent,
                                 midazimuth=True, midrange=True)[0, 0]
    dopshift = np.round(dopcentroid/prf*pershape[0]).astype('int32')
    ###########################################################################
    # Compute/get azimuth/range fourier profiles
    # NOT IMPLEMENTED
    ###########################################################################
    # azimuth_profile = ...
    # range_profile = ...
    # profile = np.outer(azimuth_profile, range_profile)
    ###########################################################################
    # Detrend image
    # NOT IMPLEMENTED / DISABLED
    ###########################################################################
    # filterwidth = 800.
    # azidetwidth = filterwidth/grdspacing[0]
    # randetwidth = filterwidth/grdspacing[1]*3/4
    # lowint = lowpass(abs(image)**2, (azidetwidth, randetwidth))
    ###########################################################################
    # Compute window
    ###########################################################################
    intensity = np.float64(abs(image)**2).mean()
    aziwindow = np.hanning(pershape[0]+2)[1:-1]
    ranwindow = np.hanning(pershape[1]+2)[1:-1]
    window = np.sqrt(np.outer(aziwindow, ranwindow)/intensity)
    ###########################################################################
    # Compute periodograms, extract looks and compute co/cross spectra
    ###########################################################################
    look = np.zeros(specshape, dtype='complex64')
    U = np.zeros(np.hstack((specshape, nlooks)), dtype='float32')
    V = np.zeros(nlooks, dtype='float32')
    looks = np.zeros(np.hstack((specshape, nlooks)), dtype='complex64')
    specs = np.zeros(np.hstack((specshape, nlooks, nlooks)), dtype='complex64')
    for appos in iter(perpos[0]):
        for rppos in iter(perpos[1]):
            sub = image[appos:appos+pershape[0], rppos:rppos+pershape[1]]
            per = fftshift(fft2(sub*window))
            #if True:
            if sarim.get_info('number_of_bursts') != 0:
                #dopshift = abs(per).mean(1).argmax()-pershape[0]/2
                dopshift = estimate_dopplershift(abs(per).mean(axis=1))
                # print (np.round(dopcentroid/prf*pershape[0]),
                #        abs(per).mean(1).argmax()-pershape[0]/2, dopshift)
            per = np.roll(per, -dopshift, axis=0)
            #per *= profile
            if debug == 'periodogram':
                fig = plt.figure()
                plt.imshow(abs(per), origin='lower', interpolation='nearest',
                           cmap=plt.cm.Greys_r)
                fig = plt.figure()
                plt.plot(abs(per).mean(axis=0))
                fig = plt.figure()
                plt.plot(abs(per).mean(axis=1))
                plt.show()
                pdb.set_trace()
            for n, alpos in enumerate(lookpos[0]):
                look[0:lookshape[0], 0:lookshape[1]] = \
                    per[alpos:alpos+lookshape[0],
                        lookpos[1]:lookpos[1]+lookshape[1]]
                tmp = abs(look)**2
                U[:, :, n] += tmp
                V[n] += (tmp**2).sum()
                looks[:, :, n] = fft2(abs(ifft2(look))**2)
            for n in np.arange(nlooks):
                for m in np.arange(n+1):
                    specs[:, :, n, m] += looks[:, :, n]*np.conj(looks[:, :, m])
    ###########################################################################
    # Remove speckle noise bias in cospectra
    # NOT IMPLEMENTED
    ###########################################################################
    hx = np.zeros(specshape[0], dtype='float32')
    hx[[0, 1, -1]] = [0.5, 0.25, 0.25]
    hy = np.zeros(specshape[1], dtype='float32')
    hy[[0, 1, -1]] = [0.5, 0.25, 0.25]
    # hxy = np.outer(hx**2, hy**2)/(hx[0]*hy[0])**2
    # for n in np.arange(nlooks):
    #     T = (fft2(abs(ifft2(U[:, :, n]))**2)/pershape.prod()).real
    #     T = T - (T[0, 0]-0.5*V[n]/pershape.prod())*hxy
    #     #remove_bias ...
    ###########################################################################
    # Remove DC-value
    ###########################################################################
    center = specs[0, 0, :, :].copy()
    specs[np.tile([0, 1, -1], 3), np.repeat([0, 1, -1], 3), :, :] = 0
    if debug == 'specs':
        fig = plt.figure()
        plt.imshow(fftshift(specs[:, :, 1, 0].real), origin='lower',
                   interpolation='nearest', cmap=getColorMap(rgbFile='wind.pal'))
        fig = plt.figure()
        plt.imshow(fftshift(specs[:, :, 1, 0].imag), origin='lower',
                   interpolation='nearest', cmap=plt.cm.PuOr)
        plt.show()
        pdb.set_trace()
    ###########################################################################
    # Weight
    ###########################################################################
    ccs = np.zeros(np.hstack((specshape, nlooks)), dtype='complex64')
    if weighted == True:
        # if weight is complex64 then intValue may contain inf value
        weight = np.zeros(specshape, dtype='complex128')
        for n in np.arange(nlooks):
            weight = weight + specs[:, :, n, n]
        weight = abs(np.sum(weight, axis=1))
        weight = np.tile(weight.reshape([specshape[0], 1]), [1, specshape[1]])
        intValue = np.zeros(nlooks, dtype='float32')
        for n in np.arange(nlooks):
            intValue[n] = (np.sum(weight*specs[:, :, n, n])/center[n, n]).real
        rho = np.sqrt(abs(intValue/intValue.max()))
        norm = np.zeros(nlooks, dtype='float32')
        for n in np.arange(nlooks):
            for m in np.arange(n+1):
                ccs[:, :, n-m] = ccs[:, :, n-m] + \
                                 specs[:, :, n, m]*rho[n]*rho[m]/center[n, m]
                norm[n-m] = norm[n-m] + (rho[n]*rho[m])**2
        hfac = np.outer(hx**2, hy**2).sum()/(hx[0]**2*hy[0]**2)
        #hfac = 1.
        for n in np.arange(nlooks):
            ccs[:, :, n] = ccs[:, :, n]/norm[n]/hfac
    else: # TMP Mean spectrum for each time sep TMP
        for n in np.arange(nlooks):
            for m in np.arange(n+1):
                ccs[:, :, n-m] = ccs[:, :, n-m] + specs[:, :, n, m]
        for n in np.arange(nlooks):
            ccs[:, :, n] = ccs[:, :, n]/(nlooks-n)
    ###########################################################################
    # Shift spectra
    ###########################################################################
    for n in np.arange(nlooks):
        ccs[:, :, n] = fftshift(ccs[:, :, n])
    if debug == 'ccs':
        fig = plt.figure()
        plt.imshow(ccs[:, :, 1].real, origin='lower', interpolation='nearest',
                   cmap=getColorMap(rgbFile='wind.pal'))
        fig = plt.figure()
        plt.imshow(ccs[:, :, 1].imag, origin='lower', interpolation='nearest',
                   cmap=plt.cm.PuOr)
        plt.show()
        pdb.set_trace()
    ###########################################################################
    # Create and return SARCCS object
    ###########################################################################
    dk = 2.*np.pi/grdspacing/pershape
    lon = sarim.get_data('lon', extent=extent, midrange=True, midazimuth=True)[0, 0]
    lat = sarim.get_data('lat', extent=extent, midrange=True, midazimuth=True)[0, 0]
    radist = (extent[3]-extent[1]+1)*grdspacing[1]
    azdist = (extent[2]-extent[0]+1)*grdspacing[0]
    inc = sarim.get_data('incidence', extent=extent, midrange=True, midazimuth=True)[0, 0]
    head = sarim.get_info('platform_heading')
    infos = OD([('extent', extent),
                ('azimuth_size', specshape[0]),
                ('range_size', specshape[1]),
                ('azimuth_periodo_size', pershape[0]),
                ('range_periodo_size', pershape[1]),
                ('azimuth_look_size', lookshape[0]),
                ('range_look_size', lookshape[1]),
                ('azimuth_look_width', azi_look_width),
                ('nlooks', nlooks),
                ('look_separation', look_sep),
                ('range_look_width', ran_look_width),
                ('azimuth_dk', dk[0]), ('range_dk', dk[1]),
                ('longitude', lon), ('latitude', lat),
                ('azimuth_dist', azdist), ('range_dist', radist),
                ('incidence', inc), ('platform_heading', head),
                ('normalized_variance', nvar)])
    sarxspec = SARXSpec(ccs, infos)
    return sarxspec


def sarimage2sarxspec_loop(sarim, azi_dist=10000., ran_dist=10000., **kwargs):
    """Transform a SARImage object into a list of SARXSpec objects.
    """
    product = sarim.get_info('product')
    if product != 'SLC':
        raise Exception('SLC product expected.')
    mode = sarim.get_info('mode')
    if mode == 'WV':
        extent = sarim.extent_max()
        sarxspec = [[sarimage2sarxspec(sarim, extent, **kwargs)]]
    elif mode in ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']:
        grdspacing = sarim.ground_spacing()
        azi_size = sarim.get_info('number_of_lines')
        nazi = np.round(azi_size*grdspacing[0]/azi_dist).astype('int32')
        azi_lim = np.round(np.linspace(0, azi_size, num=nazi+1))
        ran_size = sarim.get_info('number_of_samples')
        nran = np.round(ran_size*grdspacing[1]/ran_dist).astype('int32')
        ran_lim = np.round(np.linspace(0, ran_size, num=nran+1))
        sarxspec = []
        for iazi in range(nazi):
            sarxspec_range = []
            for iran in range(nran):
                extent = [azi_lim[iazi], ran_lim[iran], azi_lim[iazi+1]-1,
                          ran_lim[iran+1]-1]
                sarxspec_range.append(sarimage2sarxspec(sarim, extent, **kwargs))
            sarxspec.append(sarxspec_range)
    elif mode in ['IW', 'EW']:
        #(azi_dist, ran_dist) = (10000., 10000.) # meters
        nburst = sarim.get_info('number_of_bursts')
        grdspacing = sarim.ground_spacing()
        burst_list = sarim.get_info('burst_list')
        ran_start = burst_list['valid_location'][:, 1].min()
        ran_stop = burst_list['valid_location'][:, 3].max()
        ran_size = ran_stop-ran_start+1
        nran = np.round(ran_size*grdspacing[1]/ran_dist).astype('int32')
        ran_lim = np.round(np.linspace(ran_start, ran_stop+1, num=nran+1))
        sarxspec = []
        #for ibur in range(1):
        for ibur in range(nburst):
            burst_extent = sarim.extent_burst(ibur, valid=True)
            azi_start = burst_extent[0]
            azi_stop = burst_extent[2]
            azi_size = azi_stop-azi_start+1
            nazi = np.round(azi_size*grdspacing[0]/azi_dist).astype('int32')
            azi_lim = np.round(np.linspace(azi_start, azi_stop+1, num=nazi+1))
            for iazi in range(nazi):
                sarxspec_range = []
                for iran in range(nran):
                    extent = [azi_lim[iazi], ran_lim[iran],
                              azi_lim[iazi+1]-1, ran_lim[iran+1]-1]
                    #print ibur, iazi, iran, extent
                    sarxspec_range.append(sarimage2sarxspec(sarim, extent, **kwargs))
                sarxspec.append(sarxspec_range)
    else:
        raise Exception('mode not expected.')
    return sarxspec
