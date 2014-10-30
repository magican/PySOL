#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 19:52:39 2012

@author: mag
"""
import datetime

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 9, 23)
__modified__ = datetime.datetime(2012, 9, 26)
__version__  = "1.0"
__status__   = "Development"

from os import system

delta = datetime.timedelta(days=1)
startDate = datetime.date( year = 2010, month = 1, day = 2 )
stopDate  = datetime.date( year = 2010, month = 12, day = 30 )
dayList = range(0, (stopDate-startDate).days + 1, delta.days)
dateList = [ stopDate - datetime.timedelta(days=x) for x in dayList ]

areaList = tuple(['gulfstream'])
for date in dateList:
    date = date.strftime("%Y,%m,%d")
    for area in areaList:
        toCall = '/home/solab/python_scripts/modis_download.py' + ' --product="all"' + \
              ' --area="' + area + '"' + ' --stopD="' + date + '"'
        system(toCall)

#areaList = tuple(['blacksea'])
##areaList = tuple(['finngulf'])
#for date in dateList:
#    date = date.strftime("%Y,%m,%d")
#    for area in areaList:
#        toCall = '/home/solab/python_scripts/modis_download.py' + ' --product="all"' + \
#              ' --area="' + area + '"' + ' --stopD="' + date + '"'
#        system(toCall)

#~ areaList = 'barents', 'kara', 'laptev', 'east_siberian', 'chukchi'
#~ for date in dateList:
    #~ date = date.strftime("%Y,%m,%d")
    #~ for area in areaList:
        #~ toCall = '/home/mag/Documents/repos/solab/PySOL/modis_download.py' + ' --downDir="/home/mag/Downloads/MODIS/"' + ' --product="MOD021KM, MYD021KM"' + \
              #~ ' --area="' + area + '"' + ' --stopD="' + date + '"'
        #~ system(toCall)

#areaList = 'barents', 'kara', 'laptev', 'east_siberian', 'chukchi'
#for date in dateList:
#    date = date.strftime("%Y,%m,%d")
#    for area in areaList:
#        toCall = '/home/solab/python_scripts/modis_download.py' + ' --product="MOD021KM, MYD021KM"' + \
#              ' --area="' + area + '"' + ' --stopD="' + date + '"'
#        system(toCall)
