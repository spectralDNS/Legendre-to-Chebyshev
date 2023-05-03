#ifndef LEG2CHEB_H
#define LEG2CHEB_H

#include "fftw3.h"
#include <assert.h>
#include <cblas.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

enum { L2C = 0, C2L = 1, BOTH = 2 };

typedef struct {
  unsigned direction;
  size_t N;
  double *a;
  double *dn;
  double *an;
} direct_plan;

typedef struct {
  unsigned direction;
  size_t M;
  size_t N;
  size_t Nn;
  size_t L;
  size_t s;
  double **A;
  double *T;
  double *TT;
  double *Th;
  double *ThT;
  direct_plan *dplan;
} fmm_plan;

void free_fmm(fmm_plan *plan);
void free_direct(direct_plan *dplan);
fmm_plan *create_fmm(size_t N, size_t maxs, unsigned direction, size_t v);
direct_plan *create_direct(size_t N, unsigned direction);
size_t direct(const double *input_array, double *output_array,
              direct_plan *dplan, unsigned direction);
size_t execute(const double *input_array, double *output_array,
               fmm_plan *fmmplan, unsigned direction);
double tdiff_sec(struct timeval t0, struct timeval t1);
double norm(const double *a, const size_t N);
double _Lambda(const double z);
void __Lambda(const double* z, double* w, size_t N);
size_t get_number_of_blocks(const size_t level);
size_t get_h(const size_t level, const size_t L);
size_t get_number_of_submatrices(const size_t N, const size_t s,
                                 const size_t L);
void get_ij(size_t *ij, const size_t level, const size_t block, const size_t s,
            const size_t L);

#endif
