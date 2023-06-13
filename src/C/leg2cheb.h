#ifndef LEG2CHEB_H
#define LEG2CHEB_H

#include "fftw3.h"
#include <assert.h>
#include <unistd.h>
#include <cblas.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
//#include <Accelerate/Accelerate.h>
#ifdef OMP
#include <omp.h>
#endif
enum { L2C = 0, C2L = 1, BOTH = 2 };

typedef struct {
  size_t direction;
  size_t N;
  double *a;
  double *dn;
  double *an;
} direct_plan;

typedef struct {
  size_t direction;
  size_t M;
  size_t N;
  size_t Nn;
  size_t L;
  size_t s;
  double **A;
  double *T;
  double *Th;
  double *ThT;
  double *ia;
  double *oa;
  double *work;
  double **wk;
  double **ck;
  direct_plan *dplan;
} fmm_plan;

typedef struct {
  fmm_plan *fmmplan0;
  fmm_plan *fmmplan1;
  size_t N0;
  size_t N1;
  int axis;
} fmm_plan_2d;

void free_fmm(fmm_plan *plan);
void free_fmm_2d(fmm_plan_2d *plan);
void free_direct(direct_plan *dplan);
fmm_plan_2d *create_fmm_2d(size_t N0, size_t N1, int axis, size_t maxs, size_t M, size_t direction, size_t v);
fmm_plan *create_fmm(size_t N, size_t maxs, size_t M, size_t direction, size_t v);
direct_plan *create_direct(size_t N, size_t direction);
size_t direct(const double *input_array, double *output_array,
              direct_plan *dplan, size_t direction, size_t stride);
size_t execute2D(const double *input_array, double *output_array,
                 fmm_plan_2d *fmmplan2d, size_t direction);
size_t execute(const double *input_array, double *output_array,
               fmm_plan *fmmplan, size_t direction, const size_t stride);
double tdiff_sec(struct timeval t0, struct timeval t1);
double _Lambda(const double z);
void __Lambda(const double *z, double *w, size_t N);
size_t get_number_of_blocks(const size_t level);
size_t get_h(const size_t level, const size_t L);
size_t get_number_of_submatrices(const size_t N, const size_t s,
                                 const size_t L);
size_t get_total_number_of_blocks(const size_t N, const size_t s,
                                  const size_t L);
void get_ij(size_t *ij, const size_t level, const size_t block, const size_t s,
            const size_t L);
// Tests
void test_forward_backward(size_t N, size_t maxs, size_t M, double m, bool random, size_t verbose);
void test_speed(size_t N, size_t maxs, size_t repeat, size_t direction,
                size_t M, size_t verbose);
void test_direct(size_t N, size_t verbose);
void test_2_sizes(size_t N, size_t maxs, size_t verbose);
void test_forward_2d(size_t N0, size_t N1, size_t maxs, size_t verbose, size_t direction);
void test_forward_backward_2d(size_t N0, size_t N1, size_t maxs, size_t verbose);
void test_directM(size_t N, size_t repeat, size_t verbose);
void test_direct_speed(size_t N, size_t repeat, size_t direction, size_t verbose);

#endif
