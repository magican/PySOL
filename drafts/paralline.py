# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 21:09:55 2012

@author: mag
"""

from sympy import Point, Line, Segment, Symbol
from numpy import zeros, float64

def subs_point(l, val):
    """Take an arbitrary point and make it a fixed point."""
    t = Symbol('t', real=True)
    ap = l.arbitrary_point()
    return Point(ap[0].subs(t, val), ap[1].subs(t, val))

def paralline(nl=2, x0=0, y0=100, x1=33, y1=66):
    
    p1, p2 = Point(x0, y0), Point(x1, y1)
    s1 = Segment(p1, p2)
    
    # Create a new Line perpendicular to this linear entity which passes through the point p1.
    l1 = s1.perpendicular_line(p1)
    l2 = s1.perpendicular_line(p2)
    
    p1 in l1
    p2 in l2
    s1.is_perpendicular(l2)
    s1.is_perpendicular(l1)
    
    #p11 = subs_point(l1, s1.length)
    # find coords of parallel nl segments from each side of the transect
    x11, y11 = zeros(2*nl+1), zeros(2*nl+1)
    x22, y22 = zeros(2*nl+1), zeros(2*nl+1)
    j=0
    for i in range(-nl,nl+1):
        p11 = subs_point(l1, 1*i/s1.length) # divide unit segment on its length
        x111, y111 = p11.args
        x11[j], y11[j] = float64(x111), float64(y111)
        p22 = subs_point(l2, 1*i/s1.length) # divide unit segment on its length
        x222, y222 = p22.args
        x22[j], y22[j] = float64(x222), float64(y222)
        j+=1
    
#    # Checking that segments are parallel and same length
#    s2 = Segment(p11,p22)
#    s2.is_parallel(s1)
#    s1.length - s2.length
#    plt.plot([x0, x1], [y0, y1])
#    plt.plot([x11, x22],[y11, y22], 'k')

    return x11, y11, x22, y22