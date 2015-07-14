#!/usr/bin/env python
# coding=utf-8
"""
"""


from sar.utils.factory import sarimage
from sar.transform.sarimage2sarxspec import sarimage2sarxspec
import matplotlib.pyplot as plt
import os
import sys
import glob


def make_png(infile, outdir):
    """
    """
    sarim = sarimage(infile)
    extent = sarim.extent_max()
    sarxspec = sarimage2sarxspec(sarim, extent)
    fig = sarxspec.display_data()
    outbase = os.path.splitext(os.path.basename(infile))[0]+'-cross_spectra.png'
    outfile = os.path.join(outdir, outbase)
    fig.savefig(outfile, dpi=100)
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

    for infile in infiles:
        print infile
        make_png(infile, outdir)
