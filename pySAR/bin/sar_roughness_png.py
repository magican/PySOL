#!/usr/bin/env python
# coding=utf-8
"""
"""


from sar.utils.factory import sarimage
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.image import imread
import numpy as np
import os
from scipy.stats import scoreatpercentile
import sys
import glob


import pkg_resources
LOGO_PATH = pkg_resources.resource_filename('sar', 'share/logos')


def make_png(infile, outdir, vmax=None):
    """
    """
    # Get Sea Surface Roughness
    sarim = sarimage(infile)
    spacing_ra = np.ceil(sarim.get_info('number_of_samples')/600.)
    spacing_ra_m = sarim.pixels2meters(spacing_ra)[1]
    spacing = np.round(sarim.meters2pixels(spacing_ra_m))
    spacing_m = sarim.pixels2meters(spacing)
    ssr = sarim.get_data('roughness', spacing=spacing)
    lon = sarim.get_data('lon', midrange=True, midazimuth=True)[0, 0]
    lat = sarim.get_data('lat', midrange=True, midazimuth=True)[0, 0]
    inc = sarim.get_data('incidence', midrange=True, midazimuth=True)[0, 0]
    heading = sarim.get_info('platform_heading')
    north = 90 + heading # north dir in image
    im_num = sarim.get_info('image_number')
    dist = (np.array(ssr.shape)-1)*spacing_m/1000.
    # vmin = ssr[ind].min()
    # vmax = ssr[ind].max()
    vmin = scoreatpercentile(ssr[25:-25, 25:-25], 0.1)
    if vmax is None:
        vmax = scoreatpercentile(ssr[25:-25, 25:-25], 99.9)
    if sarim.get_info('pass') == 'Descending':
        ssr = ssr[::-1, ::-1]
        north = north+180
    # Make figure
    dpi = 100
    imsize = np.array(ssr.shape[::-1])
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
    plt.imshow(ssr, origin='lower', interpolation='nearest', vmin=vmin,
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
    plt.colorbar(label='sea surface roughness', cax=cbax, format='%.1e')
    meanstr = r'$\mu$=%.1f' % ssr.mean()
    cbax.text(0.5, -0.025, meanstr, ha='center', va='top')
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
    python_logo = os.path.join(LOGO_PATH, 'python-powered-w-70x28.png')
    logo = imread(python_logo)
    plt.figimage(logo, 5, figsize[1]-logo.shape[0]-5)
    odl_logo = os.path.join(LOGO_PATH, 'oceandatalab-85x32.png')
    logo = imread(odl_logo)
    plt.figimage(logo, 5, 5)
    # Save as PNG
    outbase = os.path.splitext(os.path.basename(infile))[0]+'-roughness.png'
    outfile = os.path.join(outdir, outbase)
    plt.savefig(outfile, dpi=dpi)#, bbox_inches='tight')
    plt.close()
    cmd = 'convert '+outfile+' -colors 256 '+outfile
    os.system(cmd)

if __name__ == "__main__":

    #infiles = ["/home/cercache/project/mpc-sentinel1/data/sample_data/wv_sl2/L2/ASA_WV_SL2__1_SH_20051128T002716_20051128T002717_019580_000000_F661.SAFE/measurement/asa-is2-slc-hh-20051127t235946-20051127t235947-019580-000000-001.tiff"]
    #infiles = ["/home/cercache/project/mpc-sentinel1/ipf/v220/S1IPF_V00220/working/S1S_CLNS_WV_60_second_test/S1A_WV_SL2__1_SH_20120101T043418_20120101T043421_001772_000001_0C3C.SAFE/measurement/s1a-wv1-slc-hh-20120101t043305-20120101t043308-001772-000001-001.tiff"]
    #outdir = "/local/home/gilles/tmp/pngtest"

    inp = sys.argv[1]
    if os.path.isdir(inp):
        infiles = glob.glob(os.path.join(inp, '*.tiff'))
    else:
        infiles = [inp]
    outdir = sys.argv[2]

    if len(sys.argv) > 3:
        vmax = sys.argv[3]
    else:
        vmax = None

    for infile in infiles:
        print infile
        make_png(infile, outdir, vmax=vmax)
