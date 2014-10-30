#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 18:03:53 2012

@author: mag
"""


from multiprocessing import Pool, Queue, Process
from time import time

from numpy import ones

import cmod
reload(cmod)

from cmod import rcs2wind, rcs2windPar


w = rcs2wind(sar=-0.387*ones((300,100)), cmdv=4, windir=0, theta=20*ones((300,100)))

w2 = rcs2windPar(sar=-0.387*ones((300,600)), cmdv=4, windir=0, theta=20*ones((300,600)), nprocs=4)

def factorize_naive(n):
    """ A naive factorization method. Take integer 'n', return list of
        factors.
    """
    if n < 2:
        return []
    factors = []
    p = 2
    while True:
        if n == 1:
            return factors
        r = n % p
        if r == 0:
            factors.append(p)
            n = n / p
        elif p * p >= n:
            factors.append(n)
            return factors
        elif p > 2:
            # Advance in steps of 2 over odd numbers
            p += 2
        else:
            # If p == 2, get to 3
            p += 1
    assert False, "unreachable"

def mp_factorizer(nums, nprocs):
    def worker(nums, out_q):
        """ The worker function, invoked in a process. 'nums' is a
            list of numbers to factor. The results are placed in
            a dictionary that's pushed to a queue.
        """
        outdict = {}
        for n in nums:
            outdict[n] = factorize_naive(n)
        out_q.put(outdict)
    # Each process will get 'chunksize' nums and a queue to put his out
    # dict into
    out_q = Queue()
    chunksize = int(ceil(len(nums) / float(nprocs)))
    procs = []
    for i in range(nprocs):
        p = Process(
                target=worker,
                args=(nums[chunksize * i:chunksize * (i + 1)],
                      out_q))
        procs.append(p)
        p.start()
    # Collect all results into a single result dict. We know how many dicts
    # with results to expect.
    resultdict = []
    for i in range(nprocs):
        resultdict.append(out_q.get())
    # Wait for all worker processes to finish
    for p in procs:
        p.join()
    return resultdict

    nms = arange(1,10000,1)
    res = mp_factorizer(nms,12)











# time.clock for cputime works incorrectly in the example,
# it doesn't take into account other threads
N = 1000
K = 300
def cb(r): #optional: callback function
    pass
    print r

def CostlyFunction(z):
    r = 0
    for k in xrange(1, K+2):
        r += z ** (1 / k**1.5)
    return r

currtimeSstart = time()
for i in xrange(N):
    CostlyFunction(i)

currtimeSstop = time()

currtimePstart = time()
po = Pool()
for i in xrange(N):
    res = po.apply_async(CostlyFunction,(i,),callback=cb)

po.close()
po.join()
currtimePstop = time()

print 'serial: time elapsed:', currtimeSstop - currtimeSstart
print '2: parallel: time elapsed:', currtimePstop - currtimePstart


po = Pool()
res = po.map(CostlyFunction,(i,))
po.close()
po.join()










def processing(thread_num, iq):
    while True:
        n = iq.get()
        r = 0
        for k in xrange(1, K+2):
            r += z ** (1 / k**1.5)
        iq.task_done()
        print "%d - %d" % (thread_num, n)

num_process = 12
in_queue = JoinableQueue()
t_start_fill = time()
for i in xrange(N):
    in_queue.put(i)

t_start_calc = time()
for i in xrange(num_process):
    worker = Process(target=processing, args=(i, in_queue))
    worker.daemon = True
    worker.start()

in_queue.join()

t_end = time()

print "*"*20
print t_end - t_start_fill
print t_end - t_start_calc
