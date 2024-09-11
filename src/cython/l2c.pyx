cimport l2c

#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

import numpy as np
cimport numpy as np

np.import_array()

cpdef enum:
    L2C
    C2L
    BOTH

cdef class Leg2Cheb:
    """Class for Legendre/Chebyshev transforms

    A fast multipole algorithm::

        M Mortensen "A faster multipole Legendre-Chebyshev transform", to be 
        published in SIAM Journal on Scientific Computing. See 
        https://github.com/spectralDNS/Legendre-to-Chebyshev/blob/main/FMM_paper.pdf

    Parameters
    ----------
    input_array : Numpy array of floats 
    output_array : Numpy array of floats
    M : int, optional
        Rank of hierarchical matrices
    maxs : int, optional
        Max size of smallest hierarchical matrix
    direction : int, optional
        0 - Legendre to Chebyshev
        1 - Chebyshev to Legendre
        2 - Assemble for both directions
    verbose : int, optional
        Verbosity level

    Note
    ----
    Input arrays are overwritten
    """
    cdef:
        l2c.fmm_plan* plan
        size_t N
        size_t direction
        size_t lagrange
        size_t precompute
        size_t verbose
        np.ndarray _input_array
        np.ndarray _output_array

    def __cinit__(self, input_array : np.ndarray, output_array : np.ndarray,
                  M : size_t=18, maxs : size_t=36, direction : size_t=BOTH,
                  lagrange : size_t=0, precompute : size_t=0, verbose : size_t=1):
        self.N = input_array.shape[0]
        self.direction = direction
        self.precompute = precompute
        self.verbose = verbose
        self.plan = <l2c.fmm_plan*>l2c.create_fmm(self.N, maxs, M, direction, lagrange, precompute, verbose)
        self._input_array = input_array
        self._output_array = output_array

    def __call__(self, input_array=None, output_array=None, direction=2):
        """
        Signature::

            __call__(input_array=None, output_array=None, direction=0, **kw)

        Compute transform and return output array

        Parameters
        ----------
        input_array : array, optional
            If not provided, then use internally stored array
        output_array : array, optional
            If not provided, then use internally stored array
        direction : int
            0 - Legendre to Chebyshev
            1 - Chebyshev to Legendre

        """
        assert direction in (0, 1), "Please specify direction of transform!"

        if input_array is not None:
            self._input_array[:] = input_array

        self._output_array[...] = 0
        l2c.execute(<double*>np.PyArray_DATA(self._input_array),
                    <double*>np.PyArray_DATA(self._output_array),
                    self.plan, direction, 1)

        if output_array is not None:
            output_array[...] = self._output_array[:]
            return output_array

        return self._output_array

    @property
    def input_array(self):
        return self._input_array

    @property
    def output_array(self):
        return self._output_array

    @property
    def verbose(self):
        return self.verbose

    def __del__(self):
        free_fmm(<l2c.fmm_plan>self.plan)

cpdef np.ndarray leg2cheb(input_array : np.ndarray, output_array : np.ndarray):
    """
    Compute Legendre to Chebyshev transform and return output array

    This function uses a direct method and should not be used for
    arrays much larger than 10,000.

    Parameters
    ----------
    input_array : array
    output_array : array

    """
    cdef l2c.direct_plan* plan = <l2c.direct_plan*>l2c.create_direct(input_array.shape[0], L2C, 0)
    l2c.direct(<double*>np.PyArray_DATA(input_array),
               <double*>np.PyArray_DATA(output_array),
               plan, L2C, 1)
    free_direct(<l2c.direct_plan>plan)
    return output_array

cpdef np.ndarray cheb2leg(input_array : np.ndarray, output_array : np.ndarray):
    """
    Compute Chebyshev to Legendre transform and return output array

    This function uses a direct method and should not be used for
    arrays much larger than 10,000.

    Parameters
    ----------
    input_array : array
    output_array : array

    """
    cdef l2c.direct_plan* plan = <l2c.direct_plan*>l2c.create_direct(input_array.shape[0], C2L, 0)
    l2c.direct(<double*>np.PyArray_DATA(input_array),
               <double*>np.PyArray_DATA(output_array),
               plan, C2L, 1)
    free_direct(<l2c.direct_plan>plan)
    return output_array

def Lambdav(np.ndarray x):
    """Return

    .. math::

            \Lambda(x) = \frac{\Gamma(x+\frac{1}{2})}{\Gamma(x+1)}

    Parameters
    ----------
    x : array of floats
        array can have any dimension and works elementwise. Returned
        array have same shape as input array.
    """
    cdef:
        size_t N = np.PyArray_Ravel(x, np.NPY_CORDER).shape[0]
        np.ndarray[double, ndim=1] a = np.empty(N)
    l2c.Lambdavec(<double*>np.PyArray_DATA(x), <double*>np.PyArray_DATA(a), N)
    return np.PyArray_Reshape(a, np.shape(x))
