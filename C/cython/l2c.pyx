cimport l2c

#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

import numpy as np
cimport numpy as np

np.import_array()

cdef class Leg2Cheb:
    cdef void *plan
    cdef np.ndarray _input_array
    cdef np.ndarray _output_array

    def __cinit__(self, input_array : np.ndarray, output_array : np.ndarray,
                  maxs : int, direction : int):
        cdef:
            int N = input_array.shape[0]

        self.plan = <l2c.fmm_plan*>l2c.create_fmm(N, maxs, direction)
        self._input_array = input_array
        self._output_array = output_array

    def __call__(self, input_array=None, output_array=None, direction=None):
        if input_array is not None:
            self._input_array[:] = input_array

        if direction is None:
            direction = (<l2c.fmm_plan*>self.plan).direction

        if direction == 2:
            raise RuntimeError("Class configured for transport in both directions. Pleas specify")

        self._output_array[...] = 0

        l2c.execute(<double*>np.PyArray_DATA(self._input_array),
                    <double*>np.PyArray_DATA(self._output_array),
                    <l2c.fmm_plan*>self.plan, direction)

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
        free_fmm(<l2c.fmm_plan>self.plan)
