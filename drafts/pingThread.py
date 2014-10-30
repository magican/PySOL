#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 14:44:32 2012

@author: mag
"""

import os
import re
import time
#import sys
from threading import Thread

class testit(Thread):
    def __init__ (self,ip):
        Thread.__init__(self)
        self.ip = ip
        self.status = -1
    def run(self):
        pingaling = os.popen("ping -q -c2 "+self.ip,"r")
        while 1:
            line = pingaling.readline()
            if not line: break
        igot = re.findall(testit.lifeline,line)
        if igot:
            self.status = int(igot[0])

testit.lifeline = re.compile(r"(\d) received")
report = ("No response","Partial Response","Alive")

print time.ctime()

pinglist = []

for host in range(100,110):
    ip = "10.170.0."+str(host)
    current = testit(ip)
    pinglist.append(current)
    current.start()

for pingle in pinglist:
    pingle.join()
    print "Status from ",pingle.ip,"is",report[pingle.status]

print time.ctime()