#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 14:49:57 2012

@author: mag
"""

import os
import re
import time
import sys

lifeline = re.compile(r"(\d) received")
report = ("No response","Partial Response","Alive")

print time.ctime()

for host in range(100,110):
    ip = "10.170.0."+str(host)
    pingaling = os.popen("ping -q -c2 "+ip,"r")
    print "Testing ",ip,
    sys.stdout.flush()
    while 1:
        line = pingaling.readline()
        if not line: break
    igot = re.findall(lifeline,line)
    if igot:
        print report[int(igot[0])]

print time.ctime()