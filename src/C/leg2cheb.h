#ifndef LEG2CHEB_H
#define LEG2CHEB_H

#include "fftw3.h"
#include <assert.h>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif
#include <stdio.h>
#include <cblas.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#ifdef OMP
#include <omp.h>
#endif

enum { L2C = 0, C2L = 1, BOTH = 2 };

#ifdef CLOCK_UPTIME_RAW
#define tic clock_gettime_nsec_np(CLOCK_UPTIME_RAW)
#define dtics(a, b) (double)(b - a) / 1.0E9
#define toc(a) (double)(tic - a) / 1.0E9
#else
#define tic clock()
#define dtics(a, b) (double)(b - a) / (double)CLOCKS_PER_SEC
#define toc(a) (double)(tic - a) / (double)CLOCKS_PER_SEC
#endif

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
  size_t lagrange;
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
fmm_plan_2d *create_fmm_2d(size_t N0, size_t N1, int axis, size_t maxs,
                           size_t M, size_t direction, const size_t lagrange, size_t v);
fmm_plan *create_fmm(const size_t N, const size_t maxs, const size_t M, const size_t direction,
                     const size_t lagrange, const size_t v);
direct_plan *create_direct(size_t N, size_t direction);
size_t direct(const double *input_array, double *output_array,
              direct_plan *dplan, size_t direction, size_t stride);
size_t directM(const double *input_array, double *output_array,
               fmm_plan *fmmplan, const size_t strides);
size_t execute2D(const double *input_array, double *output_array,
                 fmm_plan_2d *fmmplan2d, size_t direction);
size_t execute(const double *input_array, double *output_array,
               fmm_plan *fmmplan, size_t direction, const size_t stride);
void matvectri(const double *A, const double *x, double *b, double *w,
                const size_t m, const bool upper);
double _Lambda(const double z);
double _LambdaE(const double z);
void __Lambda(const double *z, double *w, size_t N);
size_t get_number_of_blocks(const size_t level);
size_t get_total_number_of_blocks(const size_t L);
size_t get_h(const size_t level, const size_t L);
size_t get_number_of_submatrices(const size_t level);
size_t get_total_number_of_submatrices(const size_t L);
void get_ij(size_t *ij, const size_t level, const size_t block, const size_t s,
            const size_t L);
void dct(double *input, double *output);
void dctH(double *input, double *output, size_t st, double *z, double *zpm);
void dct2(double *input, double *output);
void dctH2(double *input, double *output);

#endif
