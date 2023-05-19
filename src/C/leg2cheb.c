#include "leg2cheb.h"

long lmin(long a, long b) { return (a > b) ? b : a; }

size_t min(size_t a, size_t b) { return (a > b) ? b : a; }

size_t max(size_t a, size_t b) { return (a > b) ? a : b; }

double fmax(double a, double b) { return (a > b) ? a : b; }

double chebval(const double x, const double *c, size_t M) {
  const double x2 = 2 * x;
  double c0 = c[M - 2];
  double c1 = c[M - 1];
  for (size_t i = 3; i < M + 1; i++) {
    double tmp = c0;
    c0 = c[M - i] - c1;
    c1 = tmp + c1 * x2;
  }
  return c0 + c1 * x;
}

void __Lambda(const double *z, double *w, size_t N) {
  for (size_t i = 0; i < N; i++)
    w[i] = _Lambda(z[i]);
}

double _Lambda(const double z) {
  double a0[8] = {9.9688251374224490e-01,  -3.1149502185860763e-03,
                  2.5548605043159494e-06,  1.8781800445057731e-08,
                  -4.0437919461099256e-11, -9.0003384202201278e-13,
                  3.1032782098712292e-15,  1.0511830721865363e-16};

  if (z < 20)
    return exp(lgamma(z + 0.5) - lgamma(z + 1));
  return chebval(-1.0 + 40.0 / z, a0, 7) / sqrt(z); // Note - using only 7
}

double mxy(const double *x, const double *y) {
  return _Lambda(((*y) - (*x)) / 2) * _Lambda(((*y) + (*x)) / 2);
}

double lxy(const double *x, const double *y) {
  return -1.0 / (((*x) + (*y) + 1) * ((*y) - (*x))) *
         _Lambda(((*y) - (*x) - 2) / 2) * _Lambda(((*y) + (*x) - 1) / 2);
}

double tdiff_sec(struct timeval t0, struct timeval t1) {
  return (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) / 1000000.0;
}

size_t get_number_of_blocks(const size_t level) {
  return pow(2, level + 1) - 1;
}

size_t get_h(const size_t level, const size_t L) {
  return pow(2, L - level - 1);
}

size_t get_number_of_submatrices(const size_t N, const size_t s,
                                 const size_t L) {
  return 3 * N / (2 * s) - 3 * (L + 2);
}

void get_ij(size_t *ij, const size_t level, const size_t block, const size_t s,
            const size_t L) {
  size_t h = get_h(level, L);
  ij[0] = 2 * block * s * h;
  ij[1] = ij[0] + 2 * s * h;
}

size_t direct(const double *u, double *b, direct_plan *dplan,
              size_t direction, size_t strides) {
  size_t flops = 0;
  const size_t N = dplan->N;
  for (size_t i = 0; i < N; i++)
    b[i*strides] = 0.0;
  flops += N;
  if (direction == L2C) {
    const double *a = dplan->a;
    for (size_t n = 0; n < N / 2 + N % 2; n++) {
      const double *ap = &a[n];
      const double *cp = &u[2 * n * strides];
      const double a0 = ap[0];
      for (size_t i = 0; i < N - 2 * n; i++) {
        b[i*strides] += a0 * ap[i] * (*cp);
        cp += strides;
      }
      flops += 3 * (N - 2 * n);
    }
    b[0] /= 2;
    for (size_t i = 0; i < N; i++) {
      b[i*strides] *= M_2_PI;
    }
    flops += N;
  } else {
    double *vn = (double *)malloc(N * sizeof(double));
    const double *an = dplan->an;
    const double *dn = dplan->dn;
    vn[0] = u[0];
    for (size_t i = 1; i < N; i++) {
      vn[i] = u[i*strides] * i;
    }
    for (size_t n = 0; n < N; n++)
      b[n*strides] = sqrt(M_PI) * an[n] * vn[n];
    for (size_t n = 2; n < N; n = n + 2) {
      const double *ap = &an[n / 2];
      const double *vp = &vn[n];
      for (size_t i = 0; i < N - n; i++)
        b[i*strides] -= dn[n / 2] * ap[i] * vp[i];
      flops += 3 * (N - 2 * n);
    }
    for (size_t i = 0; i < N; i++)
      b[i*strides] *= (i + 0.5);
    flops += N;
    free(vn);
  }
  return flops;
}

size_t directM(const double *input_array, double *output_array,
                      fmm_plan *fmmplan, const size_t strides) {
  size_t s = fmmplan->s;
  size_t N = fmmplan->N;
  const double *a = fmmplan->dplan->a;
  size_t flops = 0;
  for (size_t n = 0; n < s; n++) {
    const double *ap = &a[n];
    const double *cp = &input_array[2 * n];
    const double a0 = ap[0];
    double *op = &output_array[0];
    if (strides == 1) {
      for (size_t i = 0; i < N - 2 * n; i++) {
        (*op++) += a0 * (*ap++) * (*cp++);
      }
    }
    else {
      for (size_t i = 0; i < N - 2 * n; i++) {
        (*op) += a0 * (*ap++) * (*cp++);
        op += strides;
      }
    }
    flops += (N - 2 * n) * 3;
  }

  size_t h = 2 * s;
  size_t n0 = 2 * s;
  for (size_t block = 0; block < get_number_of_blocks(fmmplan->L); block++) {
    size_t i0 = block * h;
    size_t j0 = n0 + i0;
    const long Nm = lmin(N - j0, h);
    for (size_t n = 0; n < h; n = n + 2) {
      double *vp = &output_array[i0*strides];
      const size_t n1 = (n + n0) / 2;
      const double *ap = &a[i0 + n1];
      const double a0 = a[n1];
      const double *up = &input_array[j0 + n];
      if ((long)(Nm - n) < 0)
        break;
      if (strides == 1) {
        for (size_t i = 0; i < (size_t)(Nm - n); i++) {
          (*vp++) += a0 * (*ap++) * (*up++);
        }
      }
      else {
        for (size_t i = 0; i < (size_t)(Nm - n); i++) {
          (*vp) += a0 * (*ap++) * (*up++);
          vp += strides;
        }
      }
      flops += 3 * (Nm - n);
    }
  }
  double *op = &output_array[0];
  if (strides == 1) {
    for (size_t i = 0; i < N; i++)
      (*op++) *= M_2_PI;
  }
  else {
    for (size_t i = 0; i < N; i++) {
      (*op) *= M_2_PI;
      op += strides;
    }
  }
  output_array[0] *= 0.5;
  return flops;
}

size_t directL(const double *input_array, double *output_array, fmm_plan *fmmplan, size_t strides) {
  size_t s = fmmplan->s;
  size_t N = fmmplan->N;
  double *an = fmmplan->dplan->an;
  double *dn = fmmplan->dplan->dn;
  double SPI = sqrt(M_PI);
  size_t flops = 0;
  double *op = &output_array[0];
  double *ia = &input_array[0];
  double *ap = &an[0];
  if (strides == 1) {
    for (size_t i = 0; i < N; i++)
      (*op++) += SPI * (*ap++) * (*ia++);
  }
  else {
    for (size_t i = 0; i < N; i++) {
      (*op) += SPI * (*ap++) * (*ia++);
      op += strides;
    }
  }

  for (size_t n = 1; n < s; n++) {
    ap = &an[n];
    ia = &input_array[2 * n];
    op = &output_array[0];
    if (strides == 1) {
      for (size_t i = 0; i < N - 2 * n; i++) {
        (*op++) -= dn[n] * (*ap++) * (*ia++);
      }
    }
    else {
      for (size_t i = 0; i < N - 2 * n; i++) {
        (*op) -= dn[n] * (*ap++) * (*ia++);
        op += strides;
      }
    }
    flops += (N - 2 * n) * 3;
  }

  size_t h = 2 * s;
  size_t n0 = 2 * s;
  for (size_t block = 0; block < get_number_of_blocks(fmmplan->L); block++) {
    size_t i0 = block * h;
    size_t j0 = n0 + i0;
    const long Nm = lmin(N - j0, h);
    for (size_t n = 0; n < h; n = n + 2) {
      const size_t n1 = (n + n0) / 2;
      double *vp = &output_array[i0*strides];
      ap = &an[i0 + n1];
      ia = &input_array[j0 + n];
      if ((long)(Nm - n) < 0)
        break;
      if (strides == 1) {
        for (size_t i = 0; i < (size_t)(Nm - n); i++) {
          (*vp++) -= dn[n1] * (*ap++) * (*ia++);
        }
      }
      else {
        for (size_t i = 0; i < (size_t)(Nm - n); i++) {
          (*vp) -= dn[n1] * (*ap++) * (*ia++);
          vp += strides;
        }
      }
      flops += 3 * (Nm - n);
    }
  }
  op = &output_array[0];
  if (strides == 1) {
    for (size_t i = 0; i < N; i++) {
      (*op++) *= (i + 0.5);
    }
  }
  else {
    for (size_t i = 0; i < N; i++) {
      (*op) *= (i + 0.5);
      op += strides;
    }
  }
  return flops;
}

void matvectri(const double *A, const double *x, double *b, const size_t m,
               const size_t n, const size_t lda, const bool upper) {
  if (upper == false) {
    for (size_t i = 0; i < m; i++) {
      double s = 0.0;
      const double *xp = &x[0];
      const double *ap = &A[i * lda];
      for (size_t j = 0; j < i + 1; j++)
        s += (*ap++) * (*xp++);
      (*b++) += s;
    }
  } else {
    for (size_t i = 0; i < m; i++) {
      double s = 0.0;
      const double *xp = &x[i];
      const double *ap = &A[i * lda + i];
      for (size_t j = i; j < n; j++)
        s += (*ap++) * (*xp++);
      (*b++) += s;
    }
  }
}

void matvectriZ(const double *A, const double *x, double *b, const size_t m,
                const size_t n, const size_t lda, const bool upper) {

  if (upper == false) {
    double *zp = (double *)malloc(m * sizeof(double));
    double *zm = (double *)malloc(m * sizeof(double));
    const double *xp = &x[lda];
    for (size_t i = 0; i < m; i = i + 2) {
      zp[i] = x[i] + xp[i];
      zm[i] = x[i] - xp[i];
    }
    for (size_t i = 1; i < m; i = i + 2) {
      zp[i] = x[i] - xp[i];
      zm[i] = x[i] + xp[i];
    }
    for (size_t i = 0; i < m; i++) {
      double s = 0.0;
      const double *z = i % 2 ? &zm[0] : &zp[0];
      const double *ap = &A[i * lda];
      for (size_t j = 0; j < i + 1; j++)
        s += (*ap++) * (*z++);
      (*b++) += s;
    }
    free(zp);
    free(zm);
  } else {
    double *bp = &b[lda];
    for (size_t i = 0; i < m; i++) {
      const double *ap = &A[i * lda + i];
      const double *xp = &x[i];
      for (size_t j = i; j < n - 1; j = j + 2) {
        double se = (*ap++) * (*xp++);
        double so = (*ap++) * (*xp++);
        (*b) += se + so;
        (*bp) += se - so;
      }
      double s = (*ap++) * (*xp++);
      (*b) += s;
      (*bp) += s;
      b++;
      bp++;
    }
  }
}

void matvec(const double *A, const double *x, double *b, const size_t m,
            const size_t n, const size_t transpose) {
  if (transpose == 0) {
    for (size_t i = 0; i < m; i++) {
      double s = 0.0;
      const double *xp = &x[0];
      for (size_t j = 0; j < n; j++)
        s += (*A++) * (*xp++);
      (*b++) += s;
    }
  } else {
    for (size_t i = 0; i < m; i++) {
      double s = (*x++);
      double *bp = &b[0];
      for (size_t j = 0; j < n; j++)
        (*bp++) += s * (*A++);
    }
  }
}

void vandermonde(double *T, const size_t h, const size_t M) {
  double *x = (double *)malloc(2 * h * sizeof(double));
  double *x2 = (double *)malloc(2 * h * sizeof(double));
  double *Tm = (double *)malloc(2 * h * M * sizeof(double));

  for (size_t i = 0; i < 2 * h; i++) {
    x[i] = -1 + (double)(i) / ((double)h);
    x2[i] = 2 * x[i];
    Tm[i * M] = (double)(1);
    Tm[i * M + 1] = x[i];
  }

  for (size_t m = 2; m < M; m++) {
    for (size_t i = 0; i < 2 * h; i++) {
      Tm[i * M + m] = Tm[i * M + m - 1] * x2[i] - Tm[i * M + m - 2];
    }
  }

  for (size_t i = 0; i < h; i++) {
    for (size_t m = 0; m < M; m++) {
      T[i * M + m] = Tm[2 * i * M + m];             // even
      T[(h + i) * M + m] = Tm[(2 * i + 1) * M + m]; // odd
    }
  }
  free(x);
  free(x2);
  free(Tm);
}

double sum(const double *a, const size_t N) {
  double x = 0;
  for (size_t i = 0; i < N; i++)
    x += a[i];
  return x;
}

void free_direct(direct_plan *plan) {
  if (plan->a != NULL) {
    free(plan->a);
    plan->a = NULL;
  }
  if (plan->an != NULL) {
    free(plan->an);
    plan->an = NULL;
  }
  if (plan->dn != NULL) {
    free(plan->dn);
    plan->dn = NULL;
  }
}

void free_fmm_2d(fmm_plan_2d *plan) {
  if (plan->fmmplan0 != NULL) {
    free_fmm(plan->fmmplan0);
    plan->fmmplan0 = NULL;
  }
  if (plan->fmmplan1 != NULL) {
    free_fmm(plan->fmmplan1);
    plan->fmmplan1 = NULL;
  }
}

void free_fmm(fmm_plan *plan) {
  if (plan->A[0] != NULL) {
    free(plan->A[0]);
    plan->A[0] = NULL;
  }
  if (plan->A[1] != NULL) {
    free(plan->A[1]);
    plan->A[1] = NULL;
  }
  if (plan->A != NULL) {
    free(plan->A);
    plan->A = NULL;
  }
  if (plan->T != NULL) {
    free(plan->T);
    plan->T = NULL;
  }
  if (plan->TT != NULL) {
    free(plan->TT);
    plan->TT = NULL;
  }
  if (plan->Th != NULL) {
    free(plan->Th);
    plan->Th = NULL;
  }
  if (plan->ThT != NULL) {
    free(plan->ThT);
    plan->ThT = NULL;
  }
  if (plan->ia != NULL) {
    free(plan->ia);
    plan->ia = NULL;
  }
  if (plan->oa != NULL) {
    free(plan->ia);
    plan->ia = NULL;
  }
  for (size_t level = 0; level < plan->L; level++) {
    if (plan->wk[level] != NULL) {
      free(plan->wk[level]);
    }
    if (plan->ck[level] != NULL) {
      free(plan->ck[level]);
    }
  }
  if (plan->wk != NULL) {
    free(plan->wk);
    plan->wk = NULL;
  }
  if (plan->ck != NULL) {
    free(plan->ck);
    plan->ck = NULL;
  }
  if (plan->dplan != NULL) {
    free_direct(plan->dplan);
    plan->dplan = NULL;
  }
}

direct_plan *create_direct(size_t N, size_t direction) {
  static direct_plan dplan = {0};
  if (direction == L2C | direction == BOTH) {
    double *a = (double *)malloc(N * sizeof(double));
    for (size_t i = 0; i < N; i++)
      a[i] = _Lambda(i);
    dplan.a = a;
  }
  if (direction == C2L | direction == BOTH) {
    double *dn = (double *)malloc(N / 2 * sizeof(double));
    double *an = (double *)malloc(N * sizeof(double));
    dn[0] = 0;
    an[0] = M_2_SQRTPI;
    for (size_t i = 1; i < N; i++)
      an[i] = 1 / (2 * _Lambda(i) * i * (i + 0.5));
    for (size_t i = 1; i < N / 2; i++)
      dn[i] = _Lambda((double)(2 * i - 2) / 2.0) / (2 * i);
    dplan.an = an;
    dplan.dn = dn;
  }
  dplan.direction = direction;
  dplan.N = N;
  return &dplan;
}

fmm_plan *create_fmm(size_t N, size_t maxs, size_t M, size_t direction, size_t v) {
  fftw_plan plan, plan1d;
  static fmm_plan fmmplan = {0};
  size_t Nn;
  size_t L = 1;
  size_t s;
  size_t ij[2];
  struct timeval t1, t2;
  size_t directions[2];
  size_t num_directions = 2;
  switch (direction) {
  case L2C:
    directions[0] = 0;
    num_directions = 1;
    break;
  case C2L:
    directions[0] = 1;
    num_directions = 1;
    break;
  default:
    directions[0] = 0;
    directions[1] = 1;
    break;
  }
  double **A = (double **)calloc(2, sizeof(double *));

  fmmplan.dplan = create_direct(N, direction);

  L = ceil(log2((double)N / (double)maxs)) - 2;
  if (L < 1) {
    if (v > 1)
      printf("Levels < 1. Using only direct method\n");
    fmmplan.Nn = N;
    fmmplan.A = A;
    return &fmmplan;
  }

  s = ceil((double)N / (double)pow(2, L + 2));
  Nn = s * pow(2, L + 2);
  fmmplan.M = M;
  fmmplan.L = L;
  fmmplan.N = N;
  fmmplan.Nn = Nn;
  fmmplan.s = s;
  if (v > 1) {
    printf("N %lu\n", N);
    printf("Num levels %lu\n", L);
    printf("Num submatrices %lu\n", get_number_of_submatrices(Nn, s, L));
    printf("Given max s %lu \n", maxs);
    printf("Computed s %lu \n", s);
    printf("Computed N %lu\n", Nn);
  }
  gettimeofday(&t1, 0);
  if (direction == BOTH) {
    A[0] = (double *)malloc(get_number_of_submatrices(Nn, s, L) * M * M *
                            sizeof(double));
    A[1] = (double *)malloc(get_number_of_submatrices(Nn, s, L) * M * M *
                            sizeof(double));
  } else {
    A[direction] = (double *)malloc(get_number_of_submatrices(Nn, s, L) * M *
                                    M * sizeof(double));
  }

  double *xj = (double *)malloc(M * sizeof(double));
  for (size_t i = 0; i < M; i++)
    xj[i] = cos((i + 0.5) * M_PI / M);

  double *fun = (double *)fftw_malloc(M * M * sizeof(double));
  double *fun_hat = (double *)fftw_malloc(M * M * sizeof(double));
  plan1d = fftw_plan_r2r_1d(M, fun, fun_hat, FFTW_REDFT10, FFTW_ESTIMATE);
  plan = fftw_plan_r2r_2d(M, M, fun, fun_hat, FFTW_REDFT10, FFTW_REDFT10,
                          FFTW_ESTIMATE);

  for (size_t di = 0; di < num_directions; di++) {
    size_t dir = directions[di];
    double *ap = A[dir];
    for (size_t level = L; level-- > 0;) // Reverse loop L-1, L-2, ..., 0
    {
      size_t h = s * get_h(level, L);
      for (size_t block = 0; block < get_number_of_blocks(level); block++) {
        get_ij(ij, level, block, s, L);
        for (size_t q = 0; q < 2; q++) {
          for (size_t p = 0; p < q + 1; p++) {
            for (size_t i = 0; i < M; i++) {
              double x = (2 * (ij[0] + p * h) + (xj[i] + 1) * h);
              for (size_t j = 0; j < M; j++) {
                double y = (2 * (ij[1] + (q + 0) * h) + (xj[j] + 1) * h);
                fun[i * M + j] = dir == 0 ? mxy(&x, &y) : lxy(&x, &y);
              }
            }
            fftw_execute(plan);
            for (size_t i = 0; i < M; i++) {
              for (size_t j = 0; j < M; j++) {
                if (i == 0 && j == 0)
                  fun_hat[i * M + j] /= (4 * M * M);
                else if (i == 0)
                  fun_hat[i * M + j] /= (2 * M * M);
                else if (j == 0)
                  fun_hat[i * M + j] /= (2 * M * M);
                else
                  fun_hat[i * M + j] /= (M * M);
              }
            }
            double *fh = &fun_hat[0];
            for (size_t i = 0; i < M * M; i++)
              *ap++ = *fh++;
          }
        }
      }
    }
  }
  fmmplan.A = A;
  double *T = (double *)malloc(2 * s * M * sizeof(double));
  double *TT = (double *)malloc(2 * s * M * sizeof(double));
  vandermonde(T, s, M);
  for (size_t i = 0; i < s; i++) {
    for (size_t m = 0; m < M; m++) {
      TT[m * s + i] = T[i * M + m];             // even
      TT[s * (M + m) + i] = T[(s + i) * M + m]; // odd
    }
  }
  fmmplan.T = T;
  fmmplan.TT = TT;

  double *Th = (double *)malloc(2 * M * M * sizeof(double));
  double *ThT = (double *)malloc(2 * M * M * sizeof(double));
  double *th = &Th[0];
  for (size_t q = 0; q < 2; q++) {
    for (size_t k = 0; k < M; k++) {
      for (size_t j = 0; j < M; j++) {
        fun[j] = cos(k * acos((xj[j] + 2 * q - 1) / 2));
        fun_hat[j] = 0.0;
      }
      fftw_execute(plan1d);
      *th++ = fun_hat[0] / M / 2;
      for (size_t j = 1; j < M; j++)
        *th++ = fun_hat[j] / M;
    }
  }
  for (size_t q = 0; q < 2; q++) {
    th = &Th[q * M * M];
    double *tht = &ThT[q * M * M];
    for (size_t i = 0; i < M; i++) {
      for (size_t j = 0; j < M; j++) {
        tht[i * M + j] = th[j * M + i];
      }
    }
  }

  double *ia = (double *)malloc(Nn / 2 * sizeof(double));
  double *oa = (double *)malloc(Nn / 2 * sizeof(double));
  double **wk = (double **)malloc(L * sizeof(double *));
  double **ck = (double **)malloc(L * sizeof(double *));
  for (size_t level = 0; level < L; level++) {
    size_t b = get_number_of_blocks(level);
    wk[level] = (double *)malloc(b * 2 * M * sizeof(double));
    ck[level] = (double *)malloc(b * 2 * M * sizeof(double));
  }
  fmmplan.ia = ia;
  fmmplan.oa = oa;
  fmmplan.wk = wk;
  fmmplan.ck = ck;

  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1d);
  free(xj);
  fmmplan.Th = Th;
  fmmplan.ThT = ThT;
  gettimeofday(&t2, 0);
  if (v > 1)
    printf("Initialization %2.4e s\n", tdiff_sec(t1, t2));
  return &fmmplan;
}

fmm_plan_2d *create_fmm_2d(size_t N0, size_t N1, int axis, size_t maxs, size_t M, size_t direction, size_t v) {
  static fmm_plan_2d fmmplan2d = {NULL, NULL, 0, 0, -1};
  if ( axis == 0 ) {
    fmmplan2d.fmmplan0 = create_fmm(N0, maxs, M, direction, v);
  }
  else if ( axis == 1 ) {
    fmmplan2d.fmmplan1 = create_fmm(N1, maxs, M, direction, v);
  }
  else if ( axis == -1 ) {
    fmmplan2d.fmmplan0 = create_fmm(N0, maxs, M, direction, v);
    if ( N0 == N1 ) {
      fmmplan2d.fmmplan1 = fmmplan2d.fmmplan0;
    }
    else {
     fmmplan2d.fmmplan1 = create_fmm(N1, maxs, M, direction, v);
    }
  }
  fmmplan2d.N0 = N0;
  fmmplan2d.N1 = N1;
  fmmplan2d.axis = axis;
  return &fmmplan2d;
}

size_t execute2D(const double *input_array, double *output_array,
                 fmm_plan_2d *fmmplan2d, size_t direction) {
  if ( fmmplan2d->axis == 0 || fmmplan2d->axis == -1) {
    for (size_t i = 0; i < fmmplan2d->N1; i++) {
      execute(&input_array[i], &output_array[i], fmmplan2d->fmmplan0, direction, fmmplan2d->N1);
    }
  }
  else if (fmmplan2d->axis == 1 || fmmplan2d->axis == -1) {
    for (size_t i = 0; i < fmmplan2d->N0; i++) {
      size_t N1 = fmmplan2d->N1;
      execute(&input_array[i*N1], &output_array[i*N1], fmmplan2d->fmmplan1, direction, 1);
    }
  }
}

size_t execute(const double *input_array, double *output_array,
               fmm_plan *fmmplan, size_t direction, const size_t stride) {
  size_t Nn = fmmplan->Nn;
  size_t N = fmmplan->N;
  size_t L = fmmplan->L;
  size_t s = fmmplan->s;
  size_t M = fmmplan->M;
  double *TT = fmmplan->TT;
  double *T = fmmplan->T;
  double *Th = fmmplan->Th;
  double *ThT = fmmplan->ThT;
  double *A = fmmplan->A[direction];
  size_t flops = 0;
  assert(direction == C2L | direction == L2C);

  if (TT == NULL) {
    flops = direct(input_array, output_array, fmmplan->dplan, direction, stride);
    return flops;
  }

  double *ia = fmmplan->ia;
  double *oa = fmmplan->oa;
  double **wk = fmmplan->wk;
  double **ck = fmmplan->ck;
  double *input;
  if (direction == C2L) { // Need to modify input array, so make copy
    input = (double *)malloc(N * sizeof(double));
    input[0] = input_array[0];
    input[1] = input_array[stride];
    double *w0 = &input[2];
    size_t ii = 2*stride;
    for (size_t i = 2; i < N; i++) {
      (*w0++) = input_array[ii] * i;
      ii += stride;
    }
  }

  for (size_t odd = 0; odd < 2; odd++) {
    for (size_t i = 0; i < Nn / 2; i++) {
      oa[i] = 0;
    }
    for (size_t level = 0; level < L; level++) {
      for (size_t i = 0; i < get_number_of_blocks(level) * 2 * M; i++) {
        wk[level][i] = 0;
        ck[level][i] = 0;
      }
    }
    const double *ap;
    switch (direction) {
    case L2C:
      ap = &input_array[odd*stride];
      break;

    case C2L:
      ap = &input[odd];
      break;
    }

    size_t rest = N % 2;
    double *iap = &ia[0];
    if (stride == 1) {
      for (size_t i = 0; i < N / 2 + rest * (1 - odd); i++) {
        *iap++ = *ap;
        ap += 2;
      }
    }
    else {
      for (size_t i = 0; i < N / 2 + rest * (1 - odd); i++) {
        *iap++ = *ap;
        ap += 2*stride;
      }
    }

    for (size_t i = N / 2 + rest * (1 - odd); i < Nn / 2; i++) {
      *iap++ = 0;
    }

    size_t ik = 0;
    size_t MM = M * M;
    for (size_t level = L; level-- > 0;) {
      if (level == L - 1) {
        size_t K = get_number_of_blocks(L - 1) * 2;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, K, M, s, 1.0,
                    &ia[2 * s], s, &T[odd * s * M], M, 0.0, wk[level], M);
        flops += 2 * K * s * M;
      }

      for (size_t block = 0; block < get_number_of_blocks(level); block++) {
        size_t Nd = block * 2 * M;
        double *c0 = &ck[level][Nd];
        double *wq = &wk[level][Nd];
#ifdef TRI
        if (level > 0 && block > 0) {
          int b0 = (block - 1) / 2;
          int q0 = (block - 1) % 2;
          matvectriZ(&Th[0], &wk[level][Nd], &wk[level - 1][(b0 * 2 + q0) * M],
                     M, M, M, false);
          flops += MM; //+2*M;
        }
#endif
        for (size_t q = 0; q < 2; q++) {
#ifndef TRI
          if (level > 0 && block > 0) {
            int b0 = (block - 1) / 2;
            int q0 = (block - 1) % 2;
            matvectri(&Th[q * MM], &wq[q * M],
                      &wk[level - 1][(b0 * 2 + q0) * M], M, M, M, false);
            flops += MM;
          }
#endif
          for (size_t p = 0; p < q + 1; p++) {
            cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1.0, &A[ik * MM], M,
                        &wq[q * M], 1, 1.0, &c0[p * M], 1);
            flops += 2 * MM;
            ik++;
          }
        }
      }
    }
    for (size_t level = 0; level < L - 1; level++) {
      double *c0 = ck[level];
      double *c1 = ck[level + 1];
      size_t j1 = 0;
      for (size_t block = 0; block < get_number_of_blocks(level + 1) - 1;
           block++) {
#ifdef TRI
        matvectriZ(&ThT[0], &c0[block * M], &c1[block * 2 * M], M, M, M, true);
        flops += 1.5 * MM;
#else
        for (size_t p = 0; p < 2; p++) {
          matvectri(&ThT[p * M * M], &c0[block * M], &c1[j1 * M], M, M, M,
                    true);
          j1++;
          flops += MM;
        }
#endif
      }
    }
    size_t K = get_number_of_blocks(L - 1) * 2;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, K, s, M, 1.0,
                ck[L - 1], M, &TT[odd * s * M], s, 0.0, &oa[0], s);
    flops += 2 * K * s * M;

    double *oaa = &output_array[odd*stride];
    double *oap = &oa[0];
    if (stride == 1) {
      for (size_t i = 0; i < N / 2 + rest * (1 - odd); i++) {
        *oaa += (*oap++);
        oaa += 2;
      }
    }
    else {
      const size_t s2 = 2*stride;
      for (size_t i = 0; i < N / 2 + rest * (1 - odd); i++) {
        *oaa += (*oap++);
        oaa += s2;
      }
    }
    // flops += 3*N/2;
  }

  switch (direction) {
  case L2C:
    flops += directM(input_array, output_array, fmmplan, stride);
    break;

  case C2L:
    flops += directL(input, output_array, fmmplan, stride);
    break;
  }

  if (direction == C2L) {
    free(input);
  }
  return flops;
}

void test_foreward_backward(size_t N, size_t maxs, double m, size_t verbose) {
  fmm_plan *fmmplan = create_fmm(N, maxs, 18, 2, verbose);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // Initialize some input array
  for (size_t i = 0; i < N; i++)
    input_array[i] = 1.0 / pow(i + 1, m);
  // Leg to Cheb
  size_t flops = execute(input_array, output_array, fmmplan, L2C, 1);
  double *ia = (double *)calloc(N, sizeof(double));
  // Cheb to Leg
  flops += execute(output_array, ia, fmmplan, C2L, 1);
  // Compute L2 error norm
  double error = 0.0;
  for (size_t j = 0; j < N; j++)
    error = fmax(fabs(input_array[j] - ia[j]), error);
  printf("N %6lu L inf Error = %2.4e \n", N, error);
  printf("               Flops = %lu\n\n", flops);
#ifdef TEST
  assert(error < 1e-10);
#endif
  free(ia);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_speed(size_t N, size_t maxs, size_t repeat, size_t direction,
                size_t M, size_t verbose) {
  fmm_plan *fmmplan = create_fmm(N, maxs, M, direction, verbose);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // Initialize some input array
  for (size_t i = 0; i < N; i++)
    input_array[i] = 1.0;

  struct timeval t0, t1;
  gettimeofday(&t0, 0);
  double min_time = 1e8;
  size_t flops;
  for (size_t i = 0; i < repeat; i++) {
    struct timeval r0, r1;
    gettimeofday(&r0, 0);
    for (size_t j = 0; j < N; j++)
      output_array[j] = 0.0;
    flops = execute(input_array, output_array, fmmplan, direction, 1);
    gettimeofday(&r1, 0);
    double s1 = tdiff_sec(r0, r1);
    min_time = s1 < min_time ? s1 : min_time;
  }
  gettimeofday(&t1, 0);
  printf("Timing N %6lu avg / min = %2.4e / %2.4e flops = %lu\n", N,
         tdiff_sec(t0, t1) / repeat, min_time, flops);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_direct(size_t N, size_t verbose) {
  direct_plan *dplan = create_direct(N, 2);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // Initialize some input array
  for (size_t i = 0; i < N; i++)
    input_array[i] = 1.0;

  size_t flops;
  flops = direct(input_array, output_array, dplan, L2C, 1);
  double *ia = (double *)calloc(N, sizeof(double));
  // Cheb to Leg
  flops += direct(output_array, ia, dplan, C2L, 1);
  // Compute L2 error norm
  double error = 0.0;
  for (size_t j = 0; j < N; j++) {
    error += pow(input_array[j] - ia[j], 2);
  }
  printf("L2 Error direct = %2.4e \n", sqrt(error));
#ifdef TEST
  assert(sqrt(error) < 1e-10);
#endif
  free(input_array);
  free(output_array);
  free_direct(dplan);
}

void test_2_sizes(size_t N, size_t maxs, size_t verbose) {
  fmm_plan *fmmplan = create_fmm(N, maxs, 18, 2, verbose);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // Initialize some input array
  for (size_t i = 0; i < N; i++)
    input_array[i] = 1.0;
  // Leg to Cheb
  size_t flops = execute(input_array, output_array, fmmplan, L2C, 1);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
  fmmplan = create_fmm(2 * N, maxs, 18, 2, verbose);
  input_array = (double *)calloc(2 * N, sizeof(double));
  output_array = (double *)calloc(2 * N, sizeof(double));
  // Initialize some input array
  for (size_t i = 0; i < 2 * N; i++)
    input_array[i] = 1.0;
  // Leg to Cheb
  flops = execute(input_array, output_array, fmmplan, L2C, 1);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_foreward_2d(size_t N0, size_t N1, size_t maxs, size_t verbose) {
  if (verbose > 1) {
    printf("test_foreward_2d\n");
  }

  double *input_array = (double *)calloc(N0*N1, sizeof(double));
  double *output_array = (double *)calloc(N0*N1, sizeof(double));
  double *input_array1d = (double *)calloc(max(N0, N1), sizeof(double));
  double *output_array1d = (double *)calloc(max(N0, N1), sizeof(double));

  for (size_t axis = 0; axis < 2; axis++)
  {
    // 1D first
    size_t N = (axis == 0) ? N0 : N1;
    fmm_plan *fmmplan = create_fmm(N, maxs, 18, 2, verbose);
    for (size_t i = 0; i < N; i++) {
      input_array1d[i] = 1.0;
      output_array1d[i] = 0.0;
    }
    size_t flops = execute(input_array1d, output_array1d, fmmplan, L2C, 1);

    fmm_plan_2d *fmmplan2d = create_fmm_2d(N0, N1, axis, maxs, 18, 2, verbose);
    // Initialize some input array
    for (size_t i = 0; i < N0*N1; i++) {
      input_array[i] = 1.0;
      output_array[i] = 0.0;
    }
    flops = execute2D(input_array, output_array, fmmplan2d, L2C);

    size_t strides = (axis == 0) ? N1 : 1;
    double error = 0.0;
    for (size_t j = 0; j < N; j++)
      error = fmax(fabs(output_array[j*strides] - output_array1d[j]), error);
    printf("%2.4e %2.4e %2.4e \n", output_array1d[0], output_array1d[1], output_array1d[2]);
    printf("%2.4e %2.4e %2.4e \n", output_array[0], output_array[strides], output_array[2*strides]);
    printf("N %6lu axis %lu st %lu L inf Error = %2.4e \n", N, axis, strides, error);
    printf("               Flops = %lu\n\n", flops);
#ifdef TEST
    assert(error < 1e-10);
#endif
    free_fmm(fmmplan);
  }
  free(input_array);
  free(output_array);
  free(input_array1d);
  free(output_array1d);
}
