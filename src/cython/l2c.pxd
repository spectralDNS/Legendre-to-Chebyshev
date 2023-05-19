cdef extern from "../C/leg2cheb.c":
    pass

cdef extern from "../C/leg2cheb.h":

    ctypedef struct direct_plan:
        size_t direction
        size_t N
        double* a
        double* dn
        double* an

    ctypedef struct fmm_plan:
        size_t direction
        size_t M
        size_t N
        size_t Nn
        size_t L
        size_t s
        double** A
        double* T
        double* TT
        double* Th
        double* ThT
        double *ia
        double *oa
        double **wk
        double **ck
        direct_plan* dplan

    ctypedef struct fmm_plan_2d:
        fmm_plan *fmmplan0
        fmm_plan *fmmplan1
        size_t N0
        size_t N1
        int axis

    cdef void free_fmm(fmm_plan)
    cdef void free_direct(direct_plan)
    cdef fmm_plan create_fmm(size_t, size_t, size_t, size_t, size_t)
    cdef fmm_plan_2d create_fmm_2d(size_t, size_t, int, size_t, size_t, size_t, size_t)
    cdef direct_plan create_direct(size_t, size_t)
    cdef size_t execute(const double*, double*, fmm_plan*, size_t, const size_t)
    cdef double Lambda(const double)
    cdef size_t get_number_of_blocks(const size_t)
    cdef size_t get_h(const size_t, const size_t)
    cdef size_t get_number_of_submatrices(const size_t, const size_t, const size_t)
    cdef void get_ij(size_t*, const size_t, const size_t, const size_t, const size_t)
    cdef size_t direct(const double*, double*, direct_plan*, size_t, size_t)
    cpdef double _Lambda(const double z)
    cdef void __Lambda(const double* z, double* w, size_t N)
