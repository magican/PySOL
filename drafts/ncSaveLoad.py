#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 12:03:29 2012

@author: mag
"""

from netCDF4 import Dataset
import datetime
import sys, getopt

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 5, 14)
__modified__ = datetime.datetime(2012, 5, 14)
__version__  = "1.0"
__status__   = "Development"

def Usage():
    print( "Usage: " + \
            "      " )
    return 1

def ncSave(argv=None):
    # default values
    pn = '/home/mag/'
    ncName = 'test'
    varName = 'test'
    dType = 'f8'
    units = 'test'
    data = (1,2,3,4)    
    if argv is None:
        argv = sys.argv

    if argv is None:
        print ( "Please specify the path to the RS2 folder! \n" + \
                "See USAGE for more details \n")
        return Usage()
    # Parse arguments
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["pn=","ncName="])
    except getopt.GetoptError:
        print 'readRS2.py -pn <outputfile> ...'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'readRS2.py -pn <outputfile> ...'
            sys.exit()
        elif opt in ("-pn", "--pn"):
            pn = arg
        elif opt in ("-ncName", "--ncName"):
            ncName = arg
        elif opt in ("-varName", "--varName"):
            varName = arg
        elif opt in ("-dType", "--dType"):
            dType = arg
        elif opt in ("-units", "--units"):
            units = arg
        elif opt in ("-data", "--data"):
            data = arg
    
    # open a new netCDF file for writing.
    ncfile = Dataset(pn + ncName + '.nc','w')
    # create the variable
    ncData = ncfile.createVariable(varName,dType)
    # set the units attribute.
    ncData.units =  units
    # write data to variables.
    ncData[:] = data
    # close the file.
    ncfile.close()
    print '*** SUCCESS writing file'