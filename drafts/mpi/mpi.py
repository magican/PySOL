#!/bin/env python

import pypar as p


def root():
	return p.rank() == 0

		
def broadcast( obj ):
	if root():
		for i in xrange(p.size()-1):
			p.send(obj,i+1)
		return obj
	else:
		return p.receive(0)


def broadcast_vec( vec, i ):
	myid = p.rank()
	if myid == i:
		for j in xrange(p.size()):
			if j != myid:
				p.send(vec[i],j)
	else:
		vec[i] = p.receive(i)
		
		
def scatter( vec ):
	if root():
		for i in xrange(p.size()-1):
			p.send(vec[i+1],i+1)
		return vec[0]
	else:
		return p.receive(0)

		
def gather( obj ):
	if root():
		result = [ None for i in xrange(p.size()) ]
		result[0] = obj
		for i in xrange(p.size()-1):
			result[i+1] = p.receive(i+1)
		return result
	else:
		p.send(obj,0)
	
		
def all_gather( obj ):
	myid = p.rank()
	nproc = p.size()
	result = [ None for i in xrange(nproc) ]
	result[myid] = obj
	for i in xrange(nproc):
		broadcast_vec(result,i)
	return result


def p_sum( vec ):
	nproc = p.size()
	pcs = [ None for i in xrange(nproc) ]
	if root():
		for i in xrange(nproc):
			pcs[i] = vec[i::nproc]
	temp = scatter(pcs)
	temp = sum(temp)
	pcs = gather(temp)
	if root():
		return sum(pcs)

		
def p_sum_all( vec ):
	nproc = p.size()
	pcs = [ None for i in xrange(nproc) ]
	if root():
		for i in xrange(nproc):
			pcs[i] = vec[i::nproc]
	temp = scatter(pcs)
	temp = sum(temp)
	pcs = all_gather(temp)
	return sum(pcs)


def p_dot( a, b ):
	nproc = p.size()
	va = [ None for i in xrange(nproc) ]
	vb = [ None for i in xrange(nproc) ]
	if root():
		for i in xrange(nproc):
			va[i] = a[i::nproc]
			vb[i] = b[i::nproc]
	ta = scatter(va)
	tb = scatter(vb)
	pv = [ ta[i]*tb[i] for i in xrange(len(ta)) ]
	ps = sum(pv)
	rv = gather(ps)
	if root():
		return sum(rv)

		
def p_dot_all( a, b ):
	nproc = p.size()
	va = [ None for i in xrange(nproc) ]
	vb = [ None for i in xrange(nproc) ]
	if root():
		for i in xrange(nproc):
			va[i] = a[i::nproc]
			vb[i] = b[i::nproc]
	ta = scatter(va)
	tb = scatter(vb)
	pv = [ ta[i]*tb[i] for i in xrange(len(ta)) ]
	ps = sum(pv)
	rv = all_gather(ps)
	return sum(rv)


def p_mv( m, v ):
	m = broadcast(m)
	v = broadcast(v)
	n = len(v)
	result = [ None for i in xrange(n) ]
	for i in xrange(n):
		result[i] = p_dot( m[i], v )
	return result


def empty_matrix( m, n ):
	res = [ [None for i in xrange(n)] for i in xrange(m) ]
	return res
	

def eye_matrix( m ):
	res = [ [ 0.0 for i in xrange(m) ] for i in xrange(m) ]
	for i in xrange(m):
		res[i][i] = 1.0
	return res
	
	
if __name__ == '__main__':
	
	if root():
		start = p.time()

	if False:	
		data = p.rank()
		data = broadcast(data)
		print data

	if False:
		vec = range(p.size())
		data = scatter(vec)
		print data, p.rank()
	
	if False:
		data = p.rank()
		vec = gather(data)
		if root():
			print vec
	
	if False:		
		data = p.rank()
		vec = all_gather(data)
		import time
		time.sleep(data*2+1)
		print p.rank(), vec, data

	if False:	
		v = [ 1 for i in xrange(10) ]
		res = p_sum_all(v)
		import time
		time.sleep(p.rank()*2+1)
		print p.rank(), res

	if True:		
		v = [ 2 for i in xrange(10000000) ]
		res = p_dot_all(v,v)
		#import time
		#time.sleep(p.rank()*2+1)
		print p.rank(), res

	if False:
		s = 0
		for i in xrange(100):
			r = p.rank()
			r = broadcast(r)
			s += (r + 1)
			p.barrier()
		print "%d %d" % ( p.rank(), s )

	if False:
		m = None
		v = None
		if root():
			m = eye_matrix(3000)
			v = range(3000)
		r = p_mv(m,v)
		if root():
			print r

	if root():
		end = p.time()
		total = end - start
		print "total time: %.14f" % total
			
	p.finalize()
		
		
		
