import datetime

#~ originally from https://github.com/josephmeiring/LeeFilter

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2013, 5, 1)
__modified__ = datetime.datetime(2013, 5, 1)
__version__  = "1.0"
__status__   = "Development"

import numpy as np
cimport numpy as np
from scipy.ndimage.filters import uniform_filter as boxcar
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t
import cython
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def lee(np.ndarray[double, ndim=2] Array, int N = 7, float sig = 5, mode="reflect"):
    """
        N = size of the filter box
        Sig = Estimate of Variance
        mode = mode for convolution, defaut is reflect
    """

    #width of the window
    cdef  int Delta = (N - 1) / 2     
    cdef  int n_row = Array.shape[0]    
    cdef  int n_col = Array.shape[1]

    cdef  float total

    cdef np.ndarray[DTYPE_t, ndim=2] mean = boxcar(Array, (N,N), mode=mode)
    cdef np.ndarray[DTYPE_t, ndim=2] z = np.zeros((n_row, n_col), dtype=DTYPE)

    cdef int nr, nc, v, w

    for nc in xrange(Delta, n_col - Delta):  
        for nr in xrange(Delta, n_row - Delta):
			# Without types on indices v and w, this is esesntially a pure python
            # loop, and takes ~20 sec
            total = 0.0
            for v in xrange(-Delta, Delta+1):
                for w in xrange(-Delta, Delta+1):
                    #print nc, nr, v, w, total
                    total += (Array[nr + v, nc +  w ] - mean[nr,nc])**2

            z[nr, nc] = total

    z = z / (N**2 - 1)
    
    #Upon starting the next equation,  Z = Var(Z). Upon exit, Z = Var(X) 
    #of equation 19 of Lee, Optical Engineering 25(5), 636-643 (May 1986)
    #
    #VAR_X = (VAR_Z + Mean^2 )/(Sigma^2 +1) - Mean^2   (19)
    #
    #Here we constrain to >= 0, because Var(x) can't be negative:
    cdef np.ndarray[DTYPE_t, ndim=2] var_x=(z + mean**2) /(sig**2 + 1.0) - mean**2

    #return value from equation 21,22 of Lee, 1986.
    #K = ( VAR_X/(mean^2 * Sigma^2 + VAR_X) )          (22)
    #Filtered_Image = Mean + K * ( Input_Image - Mean) (21)
  
    #cdef np.ndarray[DTYPE_t, ndim=2] out_array = mean + (Array - mean) * ( var_x/(mean_squared*sig**2 + var_x) ) 

    return mean + (Array - mean) * ( var_x/(mean**2 * sig**2 + var_x) )
