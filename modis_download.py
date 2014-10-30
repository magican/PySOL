#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 12:14:02 2012

@author: mag
"""
import datetime

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2012, 9, 10)
__modified__ = datetime.datetime(2012, 9, 10)
__version__  = "1.0"
__status__   = "Development"

import urllib2
import os, sys, getopt

def downloadString(url):
    """ Returns a string with the url contents """
    filein = urllib2.urlopen(url)
    data   = filein.read()
    filein.close()
    return data

def main( argv=None ):

    dwnld = ''
    downDir = '/nfs1/store/satellite/modis/'
#    downDir = '/home/mag/Downloads/MODIS'
#    downDir = '/media/magDesk/media/SOLabNFS2/modis/chukchi/'

    # Initiate download part
    children = {}
    maxjobs = 3                 # maximum number of concurrent jobs
    jobs = []                   # current list of queued jobs
    # Default wget options to use for downloading each URL
    wget = ["wget", "-q", "-r", "-nH", "-nc", "--cut-dirs=0"]

    if argv is None:
        argv = sys.argv

    if argv is None:
        print ( "Please specify the path to the RS2 folder! \n" + \
                "See USAGE for more details \n")
        return

    # Parse arguments
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["product=","downDir=","area=","stopD="])
    except getopt.GetoptError:
        print 'readRS2.py -product <inputproduct> ...'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'readRS2.py -product <inputproduct> ...'
            sys.exit()
        elif opt in ("-product", "--product"):
            product = arg
        elif opt in ("-downDir", "--downDir"):
            downDir = arg
        elif opt in ("-area", "--area"):
            area = arg
            downDir = downDir + area + '/'
        elif opt in ("-stopD", "--stopD"):
            if len(arg)==10:
                now  = datetime.date(int(arg[0:4]), int(arg[5:7]), int(arg[8:10]))
            elif arg == 'today':
                # set the time interval for today if no now os passed in sys.argv
                now = datetime.date.today()
            elif arg == 'last': # download 10 days before last downloaded day
                dwnld = 'last'

    if area == 'blacksea':
        north = '48'
        south = '40'
        west  = '27'
        east  = '42'
    elif area == 'barents':
        north = '90'
        south = '65'
        west  = '25'
        east  = '50'
    elif area == 'chukchi':
        north = '75'
        south = '65'
        west  = '170'
        east  = '-160'
    elif area == 'finngulf':
        north = '60.8'
        south = '59.2'
        west  = '22.9'
        east  = '30.3'
    elif area == 'kara':
        north = '85'
        south = '70'
        west  = '50'
        east  = '100'
    elif area == 'laptev':
        north = '85'
        south = '70'
        west  = '100'
        east  = '150'
    elif area == 'east_siberian':
        north = '85'
        south = '70'
        west  = '150'
        east  = '170'
    elif area == 'gulfstream':
        north = '33'
        south = '23'
        west  = '-89'
        east  = '-73'

    if product == 'all':
        productList = 'MYD03', 'MOD03', 'MOD35_L2', 'MYD35_L2', \
                      'MOD021KM', 'MYD021KM', 'MOD02QKM', 'MYD02QKM'
    elif len(product) == 0:
        productList = 'MYD03', 'MOD03', 'MOD35_L2', 'MYD35_L2', \
                      'MOD021KM', 'MYD021KM', 'MOD02QKM', 'MYD02QKM'
    else:
        productList = product.split()
        if len(productList)>=2:
            productList[0] = productList[0][0:-1]
        productList = tuple(productList)

    for product in productList:
        if dwnld == 'last':
            # list the download dir to find the earliest date downloaded
            subDirs = []
            for dirname, dirnames, filenames in os.walk(downDir + 'allData/5/' + product + '/'):
                for subdirname in dirnames:
                    # sub = os.path.join(dirname, subdirname)
                    subDirs.append(subdirname)
            # set the stopdate
            year = subDirs[0]
            subDirs.sort()
            day = subDirs[0]
            now = datetime.datetime.strptime(year + ' ' + day, '%Y %j').date()

#        now = datetime.date(2012, 8, 1)
        delta = datetime.timedelta(days=1)
        startDate = str(now - delta)  # download starting  from 10 days before
        stopDate = str(now)  # today
        print "startDate=%s" % startDate
        print "stopDate =%s" % stopDate

        print product
        # Search for Day (D) MODIS file IDs
        searchForFiles = 'http://modwebsrv.modaps.eosdis.nasa.gov/axis2/services/' + \
                    'MODAPSservices/searchForFiles?product=' + product + \
                    '&collection=5&start=' + startDate + '&' + 'stop=' + stopDate + \
                    '&north=' + north + '&south=' + south + \
                    '&west=' + west + '&east=' + east + \
                    '&coordsOrTiles=coords&dayNightBoth=D'

        fileIds = downloadString(searchForFiles)

        # split the line converting to list
        fileIds = fileIds.split('<return>')[1:]
        for i in range(0,fileIds.__len__()):
            fileIds[i] = fileIds[i][:-9]

        # Convert to list for further URL request
        fileIds = ','.join(fileIds)

        fileIds = fileIds.split('</return>') # split and remove the text in the end
        fileIds = fileIds[0]


        # Search for MODIS file URLs
        getFileUrls = 'http://modwebsrv.modaps.eosdis.nasa.gov/axis2/services/' + \
                    'MODAPSservices/getFileUrls?fileIds=' + fileIds

        fileURLs = downloadString(getFileUrls)

        # split the line converting to list
        fileURLs = fileURLs.split('<return>')[1:]
        for i in range(0,fileURLs.__len__()):
            fileURLs[i] = fileURLs[i][:-9]

        # Convert to list for further URL request
        fileURLs = ','.join(fileURLs)

        fileURLs = fileURLs.split('</return>') # split and remove the text in the end
        fileURLs = fileURLs[0]
        fileURLs = fileURLs.split(',') # split the line back to list


        # Download part
        # Build a list of wget jobs, one for each URL in our input file(s).
        # Spawn a new child from jobs[] and record it in children{} using
        # its PID as a key.
        def spawn(cmd, *args):
            argv = [cmd] + list(args)
            pid = None
            try:
                pid = os.spawnlp(os.P_NOWAIT, cmd, *argv)
                children[pid] = {'pid': pid, 'cmd': argv}
            except Exception, inst:
                print "'%s': %s" % ("\x20".join(argv), str(inst))
            print "spawned pid %d of nproc=%d njobs=%d for '%s'" % \
                (pid, len(children), len(jobs), "\x20".join(argv))
            return pid

        try:
            for u in fileURLs:
                cmd = wget + [u.strip('\r\n')] + ['-P'] + [downDir]
                jobs.append(cmd)
        except IOError:
            pass

        print "%d wget jobs queued" % len(jobs)

        # Spawn at most maxjobs in parallel.
        while len(jobs) > 0 and len(children) < maxjobs:
            cmd = jobs[0]
            if spawn(*cmd):
                del jobs[0]

        print "%d jobs spawned" % len(children)

        # Watch for dying children and keep spawning new jobs while
        # we have them, in an effort to keep <= maxjobs children active.
        while len(jobs) > 0 or len(children):
            (pid, status) = os.wait()
#            print "pid %d exited. status=%d, nproc=%d, njobs=%d, cmd=%s" % \
#                (pid, status, len(children) - 1, len(jobs), \
#                 "\x20".join(children[pid]['cmd']))
            del children[pid]
            if len(children) < maxjobs and len(jobs):
                cmd = jobs[0]
                if spawn(*cmd):
                    del jobs[0]

if __name__ == "__main__":
    main(sys.argv[1:])
