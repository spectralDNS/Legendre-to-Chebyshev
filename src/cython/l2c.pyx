cimport l2c

#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

import numpy as np
cimport numpy as np

np.import_array()

cdef class Leg2Cheb:
    cdef:
        void *plan
        int N
        int use_direct
        int direction
        int verbose
        np.ndarray _input_array
        np.ndarray _output_array

    def __cinit__(self, input_array : np.ndarray, output_array : np.ndarray,
                  maxs : int=36, direction : int=2, use_direct : int=500, verbose : int=1):
        self.N = input_array.shape[0]
        self.use_direct = use_direct
        self.direction = direction
        self.verbose = verbose
        if self.N < use_direct:
            self.plan = <l2c.direct_plan*>l2c.create_direct(self.N, direction)
        else:
            self.plan = <l2c.fmm_plan*>l2c.create_fmm(self.N, maxs, direction, verbose)
        self._input_array = input_array
        self._output_array = output_array

    def __call__(self, input_array=None, output_array=None, direction=None):
        if input_array is not None:
            self._input_array[:] = input_array

        if direction is None:
            direction = self.direction

        if direction == 2:
            raise RuntimeError("Pleas specify direction of transform.")

        self._output_array[...] = 0

        if self.N > self.use_direct:
            l2c.execute(<double*>np.PyArray_DATA(self._input_array),
                        <double*>np.PyArray_DATA(self._output_array),
                        <l2c.fmm_plan*>self.plan, direction)
        else:
            l2c.direct(<double*>np.PyArray_DATA(self._input_array),
                        <double*>np.PyArray_DATA(self._output_array),
                        <l2c.direct_plan*>self.plan, direction)

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

    def __dealloc__(self):
        if self.N > self.use_direct:
            free_fmm(<l2c.fmm_plan>self.plan)
        else:
            free_direct(<l2c.direct_plan>self.plan)

cpdef np.ndarray leg2cheb(input_array : np.ndarray, output_array : np.ndarray):
    cdef l2c.direct_plan* plan = <l2c.direct_plan*>l2c.create_direct(input_array.shape[0], 0)
    l2c.direct(<double*>np.PyArray_DATA(input_array),
               <double*>np.PyArray_DATA(output_array),
               plan, 0)
    free_direct(<l2c.direct_plan>plan)
    return output_array

cpdef np.ndarray cheb2leg(input_array : np.ndarray, output_array : np.ndarray):
    cdef l2c.direct_plan* plan = <l2c.direct_plan*>l2c.create_direct(input_array.shape[0], 1)
    l2c.direct(<double*>np.PyArray_DATA(input_array),
               <double*>np.PyArray_DATA(output_array),
               plan, 1)
    free_direct(<l2c.direct_plan>plan)
    return output_array

def Lambda(np.ndarray x):
    """Return

    .. math::

            \Lambda(x) = \frac{\Gamma(x+\frac{1}{2})}{\Gamma(x+1)}

    Parameters
    ----------
    x : array of floats
        array can have any dimension
    """
    cdef:
        int N = np.PyArray_Ravel(x, np.NPY_CORDER).shape[0]
        np.ndarray[double, ndim=1] a = np.empty(N)
    l2c.__Lambda(<double*>np.PyArray_DATA(x), <double*>np.PyArray_DATA(a), N)
    return np.PyArray_Reshape(a, np.shape(x))
