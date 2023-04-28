#ifndef LEG2CHEB_H
#define LEG2CHEB_H

#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <cblas.h>
#include <omp.h>
#include "fftw3.h"


typedef struct
{
    unsigned direction;
    size_t M;
    size_t N;
    size_t Nn;
    size_t L;
    size_t s;
    double** A;
    double* a;
    double* dn;
    double* an;
    double* T;
    double* TT;
    double* Th;
    double* ThT;
} fmm_plan;

void free_fmm(fmm_plan* plan);
fmm_plan* create_fmm(size_t N, size_t maxs, unsigned direction);
size_t execute(const double* input_array, double* output_array, fmm_plan* fmmplan, unsigned direction);
double tdiff_sec(struct timeval t0, struct timeval t1);
double Lambda(const double z);
size_t get_number_of_blocks(const size_t level);
size_t get_h(const size_t level, const size_t L);
size_t get_number_of_submatrices(const size_t N, const size_t s, const size_t L);
void get_ij(size_t* ij, const size_t level, const size_t block, const size_t s, const size_t L);
size_t direct(const double* u, double* b, fmm_plan* plan, unsigned direction);

#endif
