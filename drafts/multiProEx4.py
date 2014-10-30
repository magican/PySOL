#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 11:49:25 2012

@author: mag
"""

from numpy import random
import multiprocessing
import time

x = 4
y = 6

def f1(x,y):
    for i in xrange(1,400):
        j = random.normal(1,1,(600,600))
        x = x+y*j
#    print 'F1:', x
    
def f2(x,y):
        for i in xrange(1,400):
            j = random.normal(1,1,(600,600))
            x = x-y*j
#    print 'F2:', y

def main():
    print "Starting main program"
    hilo1 = multiprocessing.Process(target=f1, args=(x,y))
    hilo2 = multiprocessing.Process(target=f2, args=(x,y))
    print "Launching threads"
    hilo1.start()
    hilo2.start()
    hilo1.join()
    hilo2.join()
    print "Both threads finished"
    print "Ending program"
    print "X:",x,"Y:",y

# Important: modules should not execute
# code when you import them
if __name__ == '__main__':
    main()