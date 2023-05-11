from ctypes import cdll, POINTER, c_double, c_size_t, c_int, c_void_p, Structure
from ctypes.util import find_library

C = cdll.LoadLibrary(find_library("libleg2cheb"))

doublep = POINTER(c_double)
doublepp = POINTER(POINTER(c_double))

class direct_plan(Structure):
    __fields__ = [("direction", c_size_t),
                  ("N", c_size_t),
                  ("a", doublep),
                  ("dn", doublep),
                  ("an", doublep)]

class fmm_plan(Structure):
    __fields__ = [("direction", c_size_t),
                  ("M", c_size_t),
                  ("N", c_size_t),
                  ("Nn", c_size_t),
                  ("L", c_size_t),
                  ("s", c_size_t),
                  ("A", doublepp),
                  ("T", doublep),
                  ("TT", doublep),
                  ("Th", doublep),
                  ("ThT", doublep),
                  ("dplan", POINTER(direct_plan))]

create_fmm = C.create_fmm
create_fmm.restype = POINTER(fmm_plan)
create_fmm.argtypes = [c_size_t, c_size_t, c_size_t, c_size_t]
execute = C.execute
execute.argtypes = [c_void_p, c_void_p, POINTER(fmm_plan), c_size_t]
create_direct = C.create_direct
create_direct.restype = POINTER(direct_plan)
create_direct.argtypes = [c_size_t, c_size_t]
direct = C.direct
direct.argtypes = [c_void_p, c_void_p, POINTER(direct_plan), c_size_t]
Lambda = C._Lambda
Lambda.restype = c_double

class Leg2Cheb:
    """Class for Legendre/Chebyshev transforms

    A fast multipole algorithm similar to::

        B. K. Alpert and V. Rokhlin, A fast algorithm for the evaluation of
        Legendre expansions, 389 SIAM Journal on Scientific and Statistical
        Computing, 12 (1991), pp. 158–179, https://doi.390org/10.1137/0912009.391

    Parameters
    ----------
    input_array : Numpy array of floats
    output_array : Numpy array of floats
    maxs : int
        Max size of smallest hierarchical matrix
    direction : int
        0 - Legendre to Chebyshev
        1 - Chebyshev to Legendre
        2 - Assemble for both directions
    verbose : int
        Verbosity level
    """
    def __init__(self, input_array, output_array, maxs : int=36,
                 direction : int=2, verbose : int=1):
        self.N = input_array.shape[0]
        self.plan = create_fmm(self.N, maxs, direction, verbose)
        self._input_array = input_array
        self._output_array = output_array
        self._input_ctypes = input_array.ctypes
        self._output_ctypes = output_array.ctypes

    def __call__(self, input_array=None, output_array=None, direction=0):
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
        assert direction in (0, 1)
        if input_array is not None:
            self._input_array[...] = input_array
        self._output_array[...] = 0
        execute(self._input_ctypes, self._output_ctypes, self.plan, direction)
        if output_array is not None:
            output_array[...] = self._output_array
            return output_array
        return self._output_array

if __name__ == '__main__':
    # Test by transforming from Legendre to Chebyshev and back to Legendre
    import numpy as np
    N = 1000
    u = np.ones(N)
    v = np.zeros_like(u)
    L2C = Leg2Cheb(u.copy(), v.copy())
    L2C(u, v, direction=0)
    L2C(v, u, direction=1)
    error_f = np.linalg.norm(u-1, np.inf)
    assert error_f < 1e-8
    print(f'Errornorm ctypes wrap  : {error_f:2.6e}')