#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 21:32:13 2012

@author: mag
"""

import pygrib
from os import system, path

def main():
    """
    Download latest Hirlam field from met.no, \
    and save with correct filename according to valid field time
    """
    # set the download folder
#    downloaddir = '/data0/SERVERS/media/hyrax/data/auxdata/model/met.no/'
    downloaddir = '/media/SOLabNFS/SERVERS/media/hyrax/data/auxdata/model/met.no/'
    tempfile = downloaddir + 'DNMI-NEurope.grb'
    tempfile = downloaddir + 'HIRLAM_10kmEurope_20120422_0000.grib'

    # Download latest file
#    URL = 'http://retro.met.no/data/maritim/DNMI-NEurope.grb'

#    wget = "wget -r -nH --cut-dirs=3 -nc "
#    system("wget -nc " + URL + " -P " + downloaddir + " -o " + downloaddir + "log_hirlam.txt")

    # Get the file date
    grbs = pygrib.open(tempfile)
    grb = grbs.message(1)
    fileDate = grb.analDate
    filetime = fileDate.strftime("%Y%m%d_%H%M")
    print "File time: ", filetime
    print "File date: ", fileDate
    # Prepare for renaming
#    outfile = downloaddir + '/HIRLAM_10kmEurope_' + filetime + '.grib'

    # Update filename
#    if not path.isfile(outfile):
#        print "Moving to:\n" , outfile
#        system("mv " + tempfile + " " + outfile)

if __name__ == "__main__":
    main()