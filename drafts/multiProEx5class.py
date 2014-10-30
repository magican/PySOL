#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 18:18:05 2012

@author: mag
"""

from multiprocessing import Pool

def f(x):
    return x*x

class calculate(object):
    def run(self):
        p = Pool()
        return p.map(f, [1,2,3])

cl = calculate()
print cl.run()