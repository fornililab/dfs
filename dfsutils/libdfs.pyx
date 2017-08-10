import numpy as np
cimport numpy as np
from libc.math cimport sqrt

def _get_drmsd(np.ndarray[np.float64_t, ndim=3] ref, np.ndarray[np.float64_t, ndim=3] ens):

	cdef np.ndarray[np.float64_t, ndim=1] drmsds = np.zeros([ref.shape[0]], dtype=np.float64)

	cdef np.float64_t eln = float(ens.shape[1])*(float(ens.shape[1]-1))/2.0
	cdef np.float64_t cumsum
	cdef np.int_t i, j, frame,

	for frame in range(ref.shape[0]):
		cumsum = 0.0
		for i in range(ref.shape[1]):
			for j in range(i):

				cumsum += (sqrt(
								(ref[frame,i,0]-ref[frame,j,0])**2 +
								(ref[frame,i,1]-ref[frame,j,1])**2 +
								(ref[frame,i,2]-ref[frame,j,2])**2) 
						 - sqrt(
								(ens[frame,i,0]-ens[frame,j,0])**2 +
								(ens[frame,i,1]-ens[frame,j,1])**2 +
								(ens[frame,i,2]-ens[frame,j,2])**2))**2

		drmsds[frame] = sqrt(cumsum/eln)

	return drmsds