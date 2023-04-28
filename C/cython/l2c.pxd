cdef extern from "../leg2cheb.c":
    pass

cdef extern from "../leg2cheb.h":

    ctypedef struct fmm_plan:
        size_t direction
        size_t N
        size_t Nn
        size_t L
        size_t s
        size_t M
        double* A
        double* a
        double* T
        double* TT
        double* Th
        double* ThT

    cdef fmm_plan create_fmm(size_t, int, int)
    cdef void free_fmm(fmm_plan)
    cdef size_t execute(const double*, double*, fmm_plan*, int)
