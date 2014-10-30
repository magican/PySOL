#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 20:58:19 2012

@author: mag
"""

from multiprocessing import Process, Queue

def multiply(a,b,que=[]): #add a argument to function for assigning a queue
    que.put(a*b) #we're putting return value into queue

if __name__ == '__main__':
    queue1 = Queue() #create a queue object
    p = Process(target= multiply, args= (5,4,queue1)) #we're setting 3rd argument to queue1
    p.start()
    S = queue1.get() #and we're getting return value: 20
    print("S = %i") %S
    p.join()
    print("ok.")