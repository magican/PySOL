#!/usr/bin/env python
# coding=utf-8
"""
"""


from sarmapper import SARMapper
from sarimage import SARImage
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.image import imread
import numpy as np
import os
from scipy.stats import scoreatpercentile
import sys
import glob
import fnmatch
import pdb
import pickle
from netCDF4 import date2num, num2date


def make_png(infile, outdir):
    """
    """
    # Get DN
    sarmp = SARMapper(infile)
    sarim = SARImage(sarmp)
    orig_spacing = sarim.get_original_spacing()
    # n = np.ceil(sarim.get_info('number_of_samples')/600.)
    # spacing = orig_spacing[1]*n
    # spacing = orig_spacing[1]*10
    im = abs(sarim.get_values('digital_number', step=[10, 10]))
    final_spacing = orig_spacing*[10, 10]
    sar_size = sarim.get_full_extent()[2:4]
    extent = [sar_size[0]/2, sar_size[1]/2, sar_size[0]/2+1, sar_size[1]/2+1]
    lon = sarim.get_values('lon', extent=extent)
    lat = sarim.get_values('lat', extent=extent)
    inc = sarim.get_values('incidence', extent=extent)
    heading = sarim.get_info('platform_heading')
    north = 90+heading # north dir in image
    im_num = sarim.get_info('image_number')
    dist = (np.array(im.shape)-1)*final_spacing/1000.
    if sarim.get_info('pass') == 'Descending':
        im = im[::-1, ::-1]
        north = north+180
    # ind = im.nonzero()
    # ind = np.where(im > 50)
    # im2 = im[ind]
    # vmin = im2.min()
    # vmax = im2.max()
    imsh = im.shape
    im2 = im[imsh[0]*0.1:imsh[0]*0.9, imsh[1]*0.1:imsh[1]*0.9]
    vmin = im2.min()
    vmax = im2.max()
    # Make figure
    dpi = 100
    imsize = np.array(im.shape[::-1])
    margin = np.array(((900-imsize[0])/2, 60))
    figsize = np.array((900, imsize[1]+2*margin[1]))
    imsizeN = imsize.astype('float')/figsize
    marginN = margin.astype('float')/figsize
    #print imsize, margin, figsize
    fig = plt.figure(figsize=figsize.astype('float')/dpi, dpi=dpi)
    # imax = fig.gca()
    # imax.set_position([marginN[0], marginN[1], imsizeN[0], imsizeN[1]])
    imax = fig.add_axes([marginN[0], marginN[1], imsizeN[0], imsizeN[1]])
    imax.set_xlim(0, dist[1])
    imax.set_ylim(0, dist[0])
    plt.imshow(im, origin='lower', interpolation='nearest', vmin=vmin,
                vmax=vmax, cmap=cm.get_cmap('Greys_r'),
                extent=[0, dist[1], 0, dist[0]], aspect='auto')
    imax.set_xlabel('range distance [km]')
    imax.set_ylabel('azimuth distance [km]')
    tit = '#%03i / lon=%.2f / lat=%.2f / inc=%.2f' % (im_num, lon, lat, inc)
    imax.set_title(tit)
    #imax.set_frame_on(False)
    #imax.set_axis_off()
    # Add colorbar
    cbax = fig.add_axes([1-0.75*marginN[0], .25, 20./figsize[0], .50])
    plt.colorbar(label='digital number', cax=cbax, format='%.1e')
    meanstr = r'$\mu$=%.2f' % im2.mean()
    cbax.text(0.5, -0.025, meanstr, ha='center', va='top')
    stdstr = r'$\sigma$=%.2f' % im2.std()
    cbax.text(0.5, -0.1, stdstr, ha='center', va='top')
    minstr = 'min=%.2f' % im2.min()
    cbax.text(0.5, 1.025, minstr, ha='center', va='bottom')
    maxstr = 'max=%.2f' % im2.max()
    cbax.text(0.5, 1.1, maxstr, ha='center', va='bottom')
    # Add north
    cpsizeN = (margin[0]-margin[1])/figsize.astype('float')
    cpax = fig.add_axes([0., 0.5-cpsizeN[1]/2, cpsizeN[0], cpsizeN[1]])
    plt.arrow(.5, .5, .5*np.cos(north*np.pi/180),
              .5*np.sin(north*np.pi/180), facecolor='black',
              width=0.01, length_includes_head=True,
              head_width=0.1, head_length=0.1)
    plt.annotate('North', (.5, .5), ha='center', va='top')
    cpax.set_axis_off()
    # Add Logos
    python_logo = '/local/home/gilles/Documents/logo/python-powered-w-70x28.png'
    #python_logo = '/home/cercache/users/gguitton/Documents/logo/python-powered-w-70x28.png'
    logo = imread(python_logo)
    plt.figimage(logo, 5, figsize[1]-logo.shape[0]-5)
    odl_logo = '/local/home/gilles/Documents/logo/oceandatalab-85x32.png'
    #odl_logo = '/home/cercache/users/gguitton/Documents/logo/oceandatalab-85x32.png'
    logo = imread(odl_logo)
    plt.figimage(logo, 5, 5)
    # Save as PNG
    infileid = os.path.basename(os.path.dirname(os.path.dirname(infile)))
    outdir2 = os.path.join(outdir, infileid)
    if os.path.exists(outdir2) == False:
        os.makedirs(outdir2)
    outfile = os.path.join(outdir2, os.path.basename(infile).replace('.tiff','-digital_number.png'))
    plt.savefig(outfile, dpi=dpi)#, bbox_inches='tight')
    plt.close()
    cmd = 'convert '+outfile+' -colors 256 '+outfile
    os.system(cmd)


def make_histo_db(infile, outdir):
    """
    """
    histodb = {}
    # Fill
    sarmp = SARMapper(infile)
    sarim = SARImage(sarmp)
    # Infos
    histodb['pol'] = sarim.get_info('polarisation')
    histodb['mode'] = sarim.get_info('mode')
    histodb['swath'] = sarim.get_info('swath')
    histodb['pass'] = sarim.get_info('pass')
    sar_size = np.array(sarim.get_full_extent()[2:4])
    extent = [sar_size[0]/2, sar_size[1]/2, sar_size[0]/2+1, sar_size[1]/2+1]
    histodb['lon'] = sarim.get_values('lon', extent=extent)[0, 0]
    histodb['lat'] = sarim.get_values('lat', extent=extent)[0, 0]
    histodb['inc'] = sarim.get_values('incidence', extent=extent)[0, 0]
    time = sarim.get_values('azimuth_time', extent=extent)[0, 0]
    time = num2date(time, sarmp.read_field('time').units).timetuple()
    histodb['time'] = time.tm_yday + time.tm_hour/24. + time.tm_min/1440. + \
                      time.tm_sec/86400.
    # Full
    dn = abs(sarim.get_values('digital_number'))
    inc = sarim.get_values('incidence')
    ele = sarim.get_values('elevation')
    histodb['dnmin'] = dn.min()
    histodb['dnmax'] = dn.max()
    histodb['incmin'] = inc.min()
    histodb['incmax'] = inc.max()
    histodb['elemin'] = ele.min()
    histodb['elemax'] = ele.max()
    # Inside (without 10% borders)
    in_size = sar_size*.8 - np.mod(sar_size*.8, 100)
    start = np.round((sar_size-in_size)/2)
    dn = dn[start[0]:start[0]+in_size[0], start[1]:start[1]+in_size[1]]
    inc = inc[start[0]:start[0]+in_size[0], start[1]:start[1]+in_size[1]]
    ele = ele[start[0]:start[0]+in_size[0], start[1]:start[1]+in_size[1]]
    histodb['dninmin'] = dn.min()
    histodb['dninmax'] = dn.max()
    histodb['incinmin'] = inc.min()
    histodb['incinmax'] = inc.max()
    histodb['eleinmin'] = ele.min()
    histodb['eleinmax'] = ele.max()
    h = np.histogram(dn, bins=np.arange(46351))
    histodb['dninhisto'] = h[0].astype('int32')
    histodb['dninbins'] = h[1].astype('int32')
    ind = np.where(dn <= 100)
    histodb['dninlow'] = dn[ind]
    histodb['incinlow'] = inc[ind]
    histodb['eleinlow'] = ele[ind]
    # 10-mean
    sha = (in_size[0]/10, 10, in_size[1]/10, 10)
    histodb['dnin10'] = np.sqrt(np.power(dn, 2).reshape(sha).mean(-1).mean(1))
    extent = [start[0]+5, start[1]+5, start[0]+in_size[0]-5+1, start[1]+in_size[1]-5+1]
    histodb['incin10'] = sarim.get_values('incidence', extent=extent, step=(10,10))
    histodb['elein10'] = sarim.get_values('elevation', extent=extent, step=(10,10))
    # 50-mean
    sha = (in_size[0]/50, 50, in_size[1]/50, 50)
    histodb['dnin50'] = np.sqrt(np.power(dn, 2).reshape(sha).mean(-1).mean(1))
    extent = [start[0]+25, start[1]+25, start[0]+in_size[0]-25+1, start[1]+in_size[1]-25+1]
    histodb['incin50'] = sarim.get_values('incidence', extent=extent, step=(50,50))
    histodb['elein50'] = sarim.get_values('elevation', extent=extent, step=(50,50))
    # 100-mean
    sha = (in_size[0]/100, 100, in_size[1]/100, 100)
    histodb['dnin100'] = np.sqrt(np.power(dn, 2).reshape(sha).mean(-1).mean(1))
    extent = [start[0]+50, start[1]+50, start[0]+in_size[0]-50+1, start[1]+in_size[1]-50+1]
    histodb['incin100'] = sarim.get_values('incidence', extent=extent, step=(100,100))
    histodb['elein100'] = sarim.get_values('elevation', extent=extent, step=(100,100))
    # Save
    infileid = os.path.basename(os.path.dirname(os.path.dirname(infile)))
    outdir2 = os.path.join(outdir, infileid)
    if os.path.exists(outdir2) == False:
        os.makedirs(outdir2)
    outfile = os.path.join(outdir2, os.path.basename(infile).replace('.tiff','-histodb.pkl'))
    pickle.dump(histodb, open(outfile, 'wb'), 2) 


def make_histo_pngs():
    """
    """
    dire = "/local/home/data/S1/DN_analysis/WV/DN_histodb"
    files = []
    for root, dirnames, filenames in os.walk(dire):
        for filename in fnmatch.filter(filenames, "*.pkl"):
            files.append(os.path.join(root, filename))
    outdir = '/local/home/data/S1/DN_analysis/WV/DN_histo/'
    # WV1 HH
    f = fnmatch.filter(files, '*/s1a-wv1-slc-hh-*')
    tit = 'WV1 HH %i imagettes' % len(f)
    print tit
    make_histo_png(f, os.path.join(outdir, 's1a-wv1-slc-hh'), title=tit)
    # WV1 VV
    f = fnmatch.filter(files, '*/s1a-wv1-slc-vv-*')
    tit = 'WV1 VV %i imagettes' % len(f)
    print tit
    make_histo_png(f, os.path.join(outdir, 's1a-wv1-slc-vv'), title=tit)
    # WV2 HH
    f = fnmatch.filter(files, '*/s1a-wv2-slc-hh-*')
    tit = 'WV2 HH %i imagettes' % len(f)
    print tit
    make_histo_png(f, os.path.join(outdir, 's1a-wv2-slc-hh'), title=tit)
    # WV2 VV
    f = fnmatch.filter(files, '*/s1a-wv2-slc-vv-*')
    tit = 'WV2 VV %i imagettes' % len(f)
    print tit
    make_histo_png(f, os.path.join(outdir, 's1a-wv2-slc-vv'), title=tit)
    # By imagette
    # for f in files:
    #     fid = os.path.basename(os.path.dirname(f))
    #     fbase = os.path.basename(f).replace('-histodb.pkl','')
    #     tit = fbase
    #     print tit
    #     pngbase = os.path.join(outdir, fid, fbase)
    #     make_histo_png([f], pngbase, title=tit)
    #     #pdb.set_trace()


def make_histo_png(infiles, pngbase, title=''):
    """
    """
    if os.path.exists(os.path.dirname(pngbase)) == False:
        os.makedirs(os.path.dirname(pngbase))
    # Read db
    dninbins = np.arange(46350)
    dninhisto = np.zeros(46350)
    dninlow = np.zeros(0)
    eleinlow = np.zeros(0)
    timeinlow = np.zeros(0)
    elemin = 100
    elemax = -100
    dnin10 = np.zeros(0)
    elein10 = np.zeros(0)
    timein10 = np.zeros(0)
    dnin50 = np.zeros(0)
    elein50 = np.zeros(0)
    timein50 = np.zeros(0)
    dnin100 = np.zeros(0)
    elein100 = np.zeros(0)
    timein100 = np.zeros(0)
    for infile in infiles:
        histodb = pickle.load(open(infile, 'rb'))
        dninhisto += histodb['dninhisto']
        dninlow = np.hstack((dninlow, histodb['dninlow']))
        eleinlow = np.hstack((eleinlow, histodb['eleinlow']))
        timeinlow = np.hstack((timeinlow, np.repeat(histodb['time'], 
                                                    histodb['dninlow'].size)))
        #print histodb['elemin'], histodb['elemax']
        elemin = min(elemin, histodb['elemin'])
        elemax = max(elemax, histodb['elemax'])
        dnin10 = np.hstack((dnin10, histodb['dnin10'].flatten()))
        elein10 = np.hstack((elein10, histodb['elein10'].flatten()))
        timein10 = np.hstack((timein10, np.repeat(histodb['time'], 
                                                    histodb['dnin10'].size)))
        dnin50 = np.hstack((dnin50, histodb['dnin50'].flatten()))
        elein50 = np.hstack((elein50, histodb['elein50'].flatten()))
        timein50 = np.hstack((timein50, np.repeat(histodb['time'], 
                                                    histodb['dnin50'].size)))
        dnin100 = np.hstack((dnin100, histodb['dnin100'].flatten()))
        elein100 = np.hstack((elein100, histodb['elein100'].flatten()))
        timein100 = np.hstack((timein100, np.repeat(histodb['time'], 
                                                    histodb['dnin100'].size)))
    # Histo DN
    plt.figure(figsize=(12, 8), dpi=100)
    plt.suptitle(title)
    plt.subplot(221)
    plt.plot(dninbins[0:3000], dninhisto[0:3000])
    plt.xlabel('DN')
    plt.ylabel('Count')
    h10 = np.histogram(dnin10, bins=np.arange(3001))
    plt.subplot(222)
    plt.plot(h10[1][0:-1], h10[0])
    plt.xlabel('10x10 DN')
    plt.ylabel('Count')
    h50 = np.histogram(dnin50, bins=np.arange(3001))
    plt.subplot(223)
    plt.plot(h50[1][0:-1], h50[0])
    plt.xlabel('50x50 DN')
    plt.ylabel('Count')
    h100 = np.histogram(dnin100, bins=np.arange(3001))
    plt.subplot(224)
    plt.plot(h100[1][0:-1], h100[0])
    plt.xlabel('100x100 DN')
    plt.ylabel('Count')
    outfile = pngbase+'-histo_dn.png'
    plt.savefig(outfile, dpi=100)#, bbox_inches='tight')
    plt.close()
    cmd = 'convert '+outfile+' -colors 256 '+outfile
    os.system(cmd)
    # Histo DN zoom
    plt.figure(figsize=(12, 8), dpi=100)
    plt.suptitle(title)
    plt.subplot(221)
    plt.plot(dninbins[0:500], dninhisto[0:500])
    plt.xlabel('DN')
    plt.ylabel('Count')
    plt.subplot(222)
    h10min = np.where(h10[0] != 0)[0][0]
    h10min = np.floor(h10min/100.)*100
    h50min = np.where(h50[0] != 0)[0][0]
    h50min = np.floor(h50min/100.)*100
    plt.plot(h10[1][h50min:h50min+500], h10[0][h50min:h50min+500])
    plt.xlabel('10x10 DN')
    plt.ylabel('Count')
    plt.subplot(223)
    plt.plot(h50[1][h50min:h50min+500], h50[0][h50min:h50min+500])
    plt.xlabel('50x50 DN')
    plt.ylabel('Count')
    plt.subplot(224)
    plt.plot(h100[1][h50min:h50min+500], h100[0][h50min:h50min+500])
    plt.xlabel('100x100 DN')
    plt.ylabel('Count')
    outfile = pngbase+'-histo_dn_zoom.png'
    plt.savefig(outfile, dpi=100)#, bbox_inches='tight')
    plt.close()
    cmd = 'convert '+outfile+' -colors 256 '+outfile
    os.system(cmd)
    # Histo DN/elevation
    plt.figure(figsize=(12, 8), dpi=100)
    plt.suptitle(title)
    plt.subplot(221)
    plt.hist2d(dninlow, eleinlow, bins=[200, 200],
               range=[[0,3000],[elemin,elemax]], cmin=1)
    plt.xlabel('DN')
    plt.ylabel('Elevation Angle')
    plt.colorbar(label='Count')
    plt.subplot(222)
    plt.hist2d(dnin10, elein10, bins=[200, 200],
               range=[[0,3000],[elemin,elemax]], cmin=1)
    plt.xlabel('10x10 DN')
    plt.ylabel('Elevation Angle')
    plt.colorbar(label='Count')
    plt.subplot(223)
    plt.hist2d(dnin50, elein50, bins=[200, 100],
               range=[[0,3000],[elemin,elemax]], cmin=1)
    plt.xlabel('50x50 DN')
    plt.ylabel('Elevation Angle')
    plt.colorbar(label='Count')
    plt.subplot(224)
    plt.hist2d(dnin100, elein100, bins=[200, 50],
               range=[[0,3000],[elemin,elemax]], cmin=1)
    plt.xlabel('100x100 DN')
    plt.ylabel('Elevation Angle')
    plt.colorbar(label='Count')
    outfile = pngbase+'-histo_elev.png'
    plt.savefig(outfile, dpi=100)#, bbox_inches='tight')
    plt.close()
    cmd = 'convert '+outfile+' -colors 256 '+outfile
    os.system(cmd)
    # Histo DN/elevation zoom
    plt.figure(figsize=(12, 8), dpi=100)
    plt.suptitle(title)
    plt.subplot(221)
    plt.hist2d(dninlow, eleinlow, bins=[100, 200],
               range=[[0,100],[elemin,elemax]], cmin=1)
    plt.xlabel('DN')
    plt.ylabel('Elevation Angle')
    plt.colorbar(label='Count')
    plt.subplot(222)
    plt.hist2d(dnin10, elein10, bins=[200, 200],
               range=[[h50min,h50min+500],[elemin,elemax]], cmin=1)
    plt.xlabel('10x10 DN')
    plt.ylabel('Elevation Angle')
    plt.colorbar(label='Count')
    plt.subplot(223)
    plt.hist2d(dnin50, elein50, bins=[200, 100],
               range=[[h50min,h50min+500],[elemin,elemax]], cmin=1)
    plt.xlabel('50x50 DN')
    plt.ylabel('Elevation Angle')
    plt.colorbar(label='Count')
    plt.subplot(224)
    plt.hist2d(dnin100, elein100, bins=[200, 50],
               range=[[h50min,h50min+500],[elemin,elemax]], cmin=1)
    plt.xlabel('100x100 DN')
    plt.ylabel('Elevation Angle')
    plt.colorbar(label='Count')
    outfile = pngbase+'-histo_elev_zoom.png'
    plt.savefig(outfile, dpi=100)#, bbox_inches='tight')
    plt.close()
    cmd = 'convert '+outfile+' -colors 256 '+outfile
    os.system(cmd)
    # Histo DN/time
    if len(infiles) > 1:
        timemin = np.floor(timein10.min())
        timemax = np.ceil(timein10.max())
        ntime = (timemax-timemin)*2
        plt.figure(figsize=(12, 8), dpi=100)
        plt.suptitle(title)
        plt.subplot(221)
        plt.hist2d(dninlow, timeinlow, bins=[200, ntime],
                   range=[[0,3000],[timemin,timemax]], cmin=1)
        plt.xlabel('DN')
        plt.ylabel('Day of Year')
        plt.colorbar(label='Count')
        plt.subplot(222)
        plt.hist2d(dnin10, timein10, bins=[200, ntime],
                   range=[[0,3000],[timemin,timemax]], cmin=1)
        plt.xlabel('10x10 DN')
        plt.ylabel('Day of Year')
        plt.colorbar(label='Count')
        plt.subplot(223)
        plt.hist2d(dnin50, timein50, bins=[200, ntime],
                   range=[[0,3000],[timemin,timemax]], cmin=1)
        plt.xlabel('50x50 DN')
        plt.ylabel('Day of Year')
        plt.colorbar(label='Count')
        plt.subplot(224)
        plt.hist2d(dnin100, timein100, bins=[200, ntime],
                   range=[[0,3000],[timemin,timemax]], cmin=1)
        plt.xlabel('100x100 DN')
        plt.ylabel('Day of Year')
        plt.colorbar(label='Count')
        outfile = pngbase+'-histo_time.png'
        plt.savefig(outfile, dpi=100)#, bbox_inches='tight')
        plt.close()
        cmd = 'convert '+outfile+' -colors 256 '+outfile
        os.system(cmd)
        # Histo DN/time zoom
        plt.figure(figsize=(12, 8), dpi=100)
        plt.suptitle(title)
        plt.subplot(221)
        plt.hist2d(dninlow, timeinlow, bins=[100, ntime],
                   range=[[0,100],[timemin,timemax]], cmin=1)
        plt.xlabel('DN')
        plt.ylabel('Day of Year')
        plt.colorbar(label='Count')
        plt.subplot(222)
        plt.hist2d(dnin10, timein10, bins=[200, ntime],
                   range=[[h50min,h50min+500],[timemin,timemax]], cmin=1)
        plt.xlabel('10x10 DN')
        plt.ylabel('Day of Year')
        plt.colorbar(label='Count')
        plt.subplot(223)
        plt.hist2d(dnin50, timein50, bins=[200, ntime],
                   range=[[h50min,h50min+500],[timemin,timemax]], cmin=1)
        plt.xlabel('50x50 DN')
        plt.ylabel('Day of Year')
        plt.colorbar(label='Count')
        plt.subplot(224)
        plt.hist2d(dnin100, timein100, bins=[200, ntime],
                   range=[[h50min,h50min+500],[timemin,timemax]], cmin=1)
        plt.xlabel('100x100 DN')
        plt.ylabel('Day of Year')
        plt.colorbar(label='Count')
        outfile = pngbase+'-histo_time_zoom.png'
        plt.savefig(outfile, dpi=100)#, bbox_inches='tight')
        plt.close()
        cmd = 'convert '+outfile+' -colors 256 '+outfile
        os.system(cmd)
    # For MPC presentation
    if len(infiles) > 1:
        plt.figure(figsize=(6, 4), dpi=100)
        plt.suptitle(title+' (zoom)')
        plt.hist2d(dnin50, elein50, bins=[200, 100],
                   range=[[h50min,h50min+500],[elemin,elemax]], cmin=1)
        plt.xlabel('50x50 DN')
        plt.ylabel('Elevation Angle')
        plt.colorbar(label='Count')
        outfile = pngbase+'-mpc1.png'
        plt.savefig(outfile, dpi=100, bbox_inches='tight')
        plt.close()
        cmd = 'convert '+outfile+' -colors 256 '+outfile
        os.system(cmd)
        h50max = np.where(h50[0] != 0)[0].max()
        h50max = np.ceil(h50max/100.)*100
        plt.figure(figsize=(6, 4), dpi=100)
        plt.suptitle(title)
        plt.hist2d(dnin50, elein50, bins=[200, 100],
                   range=[[h50min,h50max],[elemin,elemax]], cmin=1)
        plt.xlabel('50x50 DN')
        plt.ylabel('Elevation Angle')
        plt.colorbar(label='Count')
        outfile = pngbase+'-mpc2.png'
        plt.savefig(outfile, dpi=100, bbox_inches='tight')
        plt.close()
        cmd = 'convert '+outfile+' -colors 256 '+outfile
        os.system(cmd)
        plt.figure(figsize=(6, 4), dpi=100)
        plt.suptitle(title+' (zoom)')
        plt.hist2d(dnin50, elein50, bins=[200, 100],
                   range=[[h50min,h50min+200],[elemin,elemax]], cmin=1)
        plt.xlabel('50x50 DN')
        plt.ylabel('Elevation Angle')
        plt.colorbar(label='Count')
        outfile = pngbase+'-mpc3.png'
        plt.savefig(outfile, dpi=100, bbox_inches='tight')
        plt.close()
        cmd = 'convert '+outfile+' -colors 256 '+outfile
        os.system(cmd)


if __name__ == "__main__":

    #infiles = ["/home/cercache/project/mpc-sentinel1/data/sample_data/wv_sl2/L2/ASA_WV_SL2__1_SH_20051128T002716_20051128T002717_019580_000000_F661.SAFE/measurement/asa-is2-slc-hh-20051127t235946-20051127t235947-019580-000000-001.tiff"]
    #infiles = ["/home/cercache/project/mpc-sentinel1/ipf/v220/S1IPF_V00220/working/S1S_CLNS_WV_60_second_test/S1A_WV_SL2__1_SH_20120101T043418_20120101T043421_001772_000001_0C3C.SAFE/measurement/s1a-wv1-slc-hh-20120101t043305-20120101t043308-001772-000001-001.tiff"]
    #outdir = "/local/home/gilles/tmp/pngtest"

    inp = sys.argv[1]
    # if os.path.isdir(inp):
    #     infiles = glob.glob(os.path.join(inp, '*.tiff'))
    # else:
    #     infiles = [inp]
    exclude = ('S1A_WV_SLC__1SSH_20140406T133433_20140406T133816_000039_7FFF80_0396.SAFE', # land
               'S1A_WV_SLC__1SSH_20140416T112316_20140416T113320_000183_0000FE_B67A.SAFE', # land
               'S1A_WV_SLC__1SSH_20140421T084408_20140421T084918_000255_00012B_DB30.SAFE', # doublon
               'S1A_WV_SLC__1SSV_20140419T090520_20140419T091031_000226_00011F_3F91.SAFE', # doublon
               'S1A_WV_SLC__1SSV_20140419T090520_20140419T091031_000226_00011F_92D6.SAFE') # doublon
    infiles = []
    for root, dirnames, filenames in os.walk(inp):
        if os.path.basename(os.path.dirname(root)) not in exclude:
            for filename in fnmatch.filter(filenames, '*.tiff'):
                infiles.append(os.path.join(root, filename))
    outdir = sys.argv[2]

    cnttot = len(infiles)
    for cnt, infile in enumerate(infiles):
        print cnt+1, cnttot
        print infile
        #pdb.set_trace()
        make_png(infile, os.path.join(outdir, 'DN_plot'))
        make_histo_db(infile, os.path.join(outdir, 'DN_histodb'))
        #pdb.set_trace()
