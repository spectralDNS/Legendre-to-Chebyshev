#include "leg2cheb.h"
#include <getopt.h>

size_t min(size_t a, size_t b) { return (a > b) ? b : a; }

size_t max(size_t a, size_t b) { return (a > b) ? a : b; }

double fmax(double a, double b) { return (a > b) ? a : b; }

void test_forward_backward(size_t N, size_t maxs, size_t M, double m,
                           bool random, size_t lagrange, size_t verbose) {
  if (verbose > 1)
    printf("test_forward_backward\n");
  fmm_plan *fmmplan = create_fmm(N, maxs, M, BOTH, lagrange, verbose);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));

  // srand48((unsigned int)time(NULL));
  srand48(1);
  // Initialize some input array
  if (random) {
    for (size_t i = 0; i < N; i++)
      input_array[i] = (2 * drand48() - 1) / pow(i + 1, m);
  } else {
    for (size_t i = 0; i < N; i++)
      input_array[i] = 0.5 / pow(i + 1, m);
  }

  // Leg to Cheb
  if (verbose > 1)
    printf("Leg2cheb\n");

  size_t flops = execute(input_array, output_array, fmmplan, L2C, 1);
  double *ia = (double *)calloc(N, sizeof(double));
  // Cheb to Leg
  flops += execute(output_array, ia, fmmplan, C2L, 1);

  // Compute L2 error norm
  double error = 0.0;
  for (size_t j = 0; j < N; j++)
    error = fmax(fabs(input_array[j] - ia[j]), error);

  printf("N %6lu L inf Error = %2.8e \n", N, error);
  printf("               Flops = %lu\n\n", flops);
#ifdef TEST
  assert(error < 1e-10);
#endif
  free(ia);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_lambda() {
  printf("Testing Lambda\n");
  const double lam[24] = {
      8.8622692545275805e-01, 6.6467019408956851e-01, 4.8465534985697706e-01,
      3.4807557771536241e-01, 2.4805479961430874e-01, 1.7608753625554593e-01,
      1.2475610011693392e-01, 8.8302073254996866e-02, 6.2469489890636082e-02,
      4.4183385549636807e-02, 3.1246185535707110e-02, 2.2095738254098853e-02,
      1.5624523170118865e-02, 1.1048374869932079e-02, 7.8124403955826070e-03,
      5.5242506546358426e-03, 3.9062425494265085e-03, 2.7621332298331754e-03,
      1.9531240686776474e-03, 1.3810676027327610e-03, 9.7656238358468511e-04,
      6.9053392484345725e-04, 4.8828123544808499e-04, 3.4526697785636500e-04};
  for (size_t i = 0; i < 24; i++) {
    double j = pow(2, i);
    double x = _Lambda(j);
    double y = _LambdaE(j);
    printf("%d %2.6e %2.6e \n", (int)j, x - lam[i], y - lam[i]);
    // assert(fabs(x-y) < 1e-10);
  }
}

void test_accuracy(size_t N, size_t maxs, size_t M, size_t lagrange,
                   size_t verbose) {
  const double a[5] = {1.7913439415860921e+00, 2.3408718378050253e+00,
                       1.8833606631960720e+00, 1.6587642683880404e+00,
                       1.4417173322290182e+00};
  const double b[5] = {5.2631578947368418e-01, 7.5187969924812026e-02,
                       1.3268465280849182e-01, 1.7768205680441512e-01,
                       2.4367824933176932e-01};

  assert(N < 1000);
  if (verbose > 1)
    printf("test_accuracy N=20 direct and N=%d FMM\n", N);
  fmm_plan *fmmplan = create_fmm(20, maxs, M, BOTH, lagrange, 1);

  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // Initialize input array
  for (size_t i = 0; i < N; i++)
    input_array[i] = 1.0;

  execute(input_array, output_array, fmmplan, L2C, 1);

  double error0 = 0;
  for (size_t i = 0; i < 5; i++) {
    error0 += pow(output_array[i] - a[i], 2);
  }
  assert(sqrt(error0) < 1e-7);

  for (size_t i = 0; i < N; i++)
    output_array[i] = 0.0;
  execute(input_array, output_array, fmmplan, C2L, 1);

  double error1 = 0;
  for (size_t i = 0; i < 5; i++) {
    error1 += pow(output_array[i] - b[i], 2);
  }
  assert(sqrt(error1) < 1e-7);

  if (verbose > 1)
    printf("Direct N=20 Error L2C: %2.6e C2L: %2.6e\n", sqrt(error0),
           sqrt(error1));

  for (size_t i = 0; i < N; i++)
    output_array[i] = 0.0;

  fmm_plan *fmmplan2 = create_fmm(N, maxs, M, BOTH, lagrange, 1);
  direct(input_array, output_array, fmmplan2->dplan, L2C, 1);

  double *out = (double *)calloc(N, sizeof(double));

  execute(input_array, out, fmmplan2, L2C, 1);
  error0 = 0.0;
  for (size_t i = 0; i < N; i++) {
    error0 += pow(output_array[i] - out[i], 2);
  }

  // assert(sqrt(error0) < 1e-7);

  direct(input_array, output_array, fmmplan2->dplan, C2L, 1);
  for (size_t i = 0; i < N; i++)
    out[i] = 0.0;
  execute(input_array, out, fmmplan2, C2L, 1);
  error1 = 0.0;
  for (size_t i = 0; i < N; i++) {
    error1 += pow(output_array[i] - out[i], 2);
  }
  // assert(sqrt(error1) < 1e-7);

  if (verbose > 1)
    printf("FMM   N=%d Error L2C: %2.6e C2L: %2.6e\n", N, sqrt(error0),
           sqrt(error1));

  // Test only planning one direction
  fmm_plan *fmmplan3 = create_fmm(N, maxs, M, L2C, lagrange, 1);
  direct(input_array, output_array, fmmplan3->dplan, L2C, 1);

  for (size_t i = 0; i < N; i++)
    out[i] = 0.0;

  execute(input_array, out, fmmplan3, L2C, 1);
  error0 = 0.0;
  for (size_t i = 0; i < N; i++) {
    error0 += pow(output_array[i] - out[i], 2);
  }
  // assert(sqrt(error0) < 1e-7);

  fmm_plan *fmmplan4 = create_fmm(N, maxs, M, C2L, lagrange, 1);
  direct(input_array, output_array, fmmplan4->dplan, C2L, 1);

  for (size_t i = 0; i < N; i++)
    out[i] = 0.0;

  execute(input_array, out, fmmplan4, C2L, 1);
  error1 = 0.0;
  for (size_t i = 0; i < N; i++) {
    error1 += pow(output_array[i] - out[i], 2);
  }
  // assert(sqrt(error1) < 1e-7);

  if (verbose > 1)
    printf("FMM   N=%d Error L2C: %2.6e C2L: %2.6e (plan one direction)\n", N,
           sqrt(error0), sqrt(error1));

  free(out);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
  free_fmm(fmmplan2);
  free_fmm(fmmplan3);
  free_fmm(fmmplan4);
}

void test_speed(size_t N, size_t maxs, size_t repeat, size_t direction,
                size_t M, size_t lagrange, size_t verbose) {
  if (verbose > 1)
    printf("test_speed %lu\n", direction);
  fmm_plan *fmmplan = create_fmm(N, maxs, M, direction, lagrange, verbose);

  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // Initialize some input array
  for (size_t i = 0; i < N; i++)
    input_array[i] = 1.0;

  uint64_t t0 = tic;
  double min_time = 1e8;
  size_t flops;
  for (size_t i = 0; i < repeat; i++) {
    for (size_t j = 0; j < N; j++)
      output_array[j] = 0.0;
    uint64_t g0 = tic;
    flops = execute(input_array, output_array, fmmplan, direction, 1);
    double s1 = toc(g0);
    min_time = s1 < min_time ? s1 : min_time;
  }
  // for (size_t i = 0; i < 5; i++)
  //{
  //   printf("%2.6e\n", output_array[i]);
  // }

  uint64_t t1 = tic;
  printf("Timing N %6lu avg / min = %2.4e / %2.4e flops = %lu\n", N,
         dtics(t0, t1) / repeat, min_time, flops);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_direct_speed(size_t N, size_t repeat, size_t direction,
                       size_t verbose) {
  if (verbose > 1)
    printf("test_direct_speed %lu\n", direction);
  direct_plan *dplan = create_direct(N, 2);

  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // Initialize some input array
  for (size_t i = 0; i < N; i++)
    input_array[i] = 1.0;

  uint64_t t0 = tic;
  double min_time = 1e8;
  size_t flops;
  for (size_t i = 0; i < repeat; i++) {
    for (size_t j = 0; j < N; j++)
      output_array[j] = 0.0;
    uint64_t r0 = tic;
    flops = direct(input_array, output_array, dplan, direction, 1);
    double s1 = toc(r0);
    min_time = s1 < min_time ? s1 : min_time;
  }
  uint64_t t1 = tic;
  printf("Timing N %6lu avg / min = %2.4e / %2.4e flops = %lu\n", N,
         dtics(t0, t1) / repeat, min_time, flops);
  free(input_array);
  free(output_array);
  free_direct(dplan);
}

void test_direct(size_t N, size_t verbose) {
  direct_plan *dplan = create_direct(N, 2);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // Initialize some input array
  // for (size_t i = 0; i < N; i++)
  //  input_array[i] = 1.0;
  for (size_t i = 0; i < N; i++)
    input_array[i] = (2 * drand48() - 1);

  size_t flops;
  flops = direct(input_array, output_array, dplan, L2C, 1);
  double *ia = (double *)calloc(N, sizeof(double));
  // Cheb to Leg
  flops += direct(output_array, ia, dplan, C2L, 1);
  // Compute L2 error norm
  double error = 0.0;
  for (size_t j = 0; j < N; j++) {
    error += pow(input_array[j] - ia[j], 2);
    // printf("%d %2.4e \n", j, input_array[j] - ia[j]);
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
  fmm_plan *fmmplan = create_fmm(N, maxs, 18, 2, 0, verbose);
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
  fmmplan = create_fmm(2 * N, maxs, 18, 2, 0, verbose);
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

void test_forward_2d(size_t N0, size_t N1, size_t maxs, size_t verbose,
                     size_t direction) {
  if (verbose > 1) {
    printf("test_forward_2d\n");
  }

  double *input_array = (double *)calloc(N0 * N1, sizeof(double));
  double *output_array = (double *)calloc(N0 * N1, sizeof(double));
  double *input_array1d = (double *)calloc(max(N0, N1), sizeof(double));
  double *output_array1d = (double *)calloc(max(N0, N1), sizeof(double));

  for (size_t axis = 0; axis < 2; axis++) {
    // 1D first
    size_t N = ((axis == 0) ? N0 : N1);
    fmm_plan *fmmplan = create_fmm(N, maxs, 18, 2, 0, verbose);
    for (size_t i = 0; i < N; i++) {
      input_array1d[i] = 1.0 + i;
      output_array1d[i] = 0.0;
    }
    size_t flops =
        execute(input_array1d, output_array1d, fmmplan, direction, 1);

    fmm_plan_2d *fmmplan2d =
        create_fmm_2d(N0, N1, axis, maxs, 18, 2, 0, verbose);
    // Initialize some input array
    for (size_t i = 0; i < N0; i++) {
      for (size_t j = 0; j < N1; j++) {
        input_array[i * N1 + j] = 1.0 + (axis == 0 ? i : j);
        output_array[i * N1 + j] = 0.0;
      }
    }
    flops = execute2D(input_array, output_array, fmmplan2d, direction);

    size_t strides = (axis == 0 ? N1 : 1);
    double error = 0.0;
    for (size_t j = 0; j < N; j++)
      error = fmax(fabs(output_array[j * strides] - output_array1d[j]), error);
    printf("%2.4e %2.4e %2.4e \n", output_array1d[0], output_array1d[1],
           output_array1d[2]);
    printf("%2.4e %2.4e %2.4e \n", output_array[0], output_array[strides],
           output_array[2 * strides]);
    printf("N %6lu axis %lu st %lu L inf Error = %2.4e \n", N, axis, strides,
           error);
    printf("               Flops = %lu\n\n", flops);
#ifdef TEST
    assert(error < 1e-10);
#endif
    free_fmm(fmmplan);
    free_fmm_2d(fmmplan2d);
  }
  free(input_array);
  free(output_array);
  free(input_array1d);
  free(output_array1d);
}

void test_forward_backward_2d(size_t N0, size_t N1, size_t maxs,
                              size_t verbose) {
  if (verbose > 1) {
    printf("test_forward_backward_2d\n");
  }

  double *input_array = (double *)calloc(N0 * N1, sizeof(double));
  double *output_array = (double *)calloc(N0 * N1, sizeof(double));
  double *output_array2 = (double *)calloc(N0 * N1, sizeof(double));

  int axis = -1;
  fmm_plan_2d *fmmplan2d = create_fmm_2d(N0, N1, axis, maxs, 18, 2, 0, verbose);
  // Initialize some input array
  for (size_t i = 0; i < N0; i++) {
    for (size_t j = 0; j < N1; j++) {
      input_array[i * N1 + j] = 1.0;
    }
  }

  size_t flops = execute2D(input_array, output_array, fmmplan2d, L2C);
  flops = execute2D(output_array, output_array2, fmmplan2d, C2L);

  double error = 0.0;
  for (size_t i = 0; i < N0; i++) {
    for (size_t j = 0; j < N1; j++) {
      error += pow(output_array2[i * N1 + j] - input_array[i * N1 + j], 2);
    }
  }
  printf("Error %2.6e\n", sqrt(error));

#ifdef TEST
  assert(sqrt(error) < 1e-8);
#endif
  free_fmm_2d(fmmplan2d);
  free(input_array);
  free(output_array);
  free(output_array2);
}

void test_directM(size_t N, size_t repeat, size_t verbose, size_t s, size_t M) {
  fmm_plan *fmmplan = create_fmm(N, s, M, 2, 0, verbose);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  size_t flops;
  double min_time = 1e8;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t r0 = tic;
    flops = directM(input_array, output_array, fmmplan, 1);
    double s1 = toc(r0);
    min_time = s1 < min_time ? s1 : min_time;
  }
  printf("directM min time %2.6e\n", min_time);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_dct(size_t N, size_t repeat) {
  double *fun = (double *)fftw_malloc(N * N * sizeof(double));
  double *fun_hat = (double *)fftw_malloc(N * N * sizeof(double));
  double *fun_hat2 = (double *)fftw_malloc(N * N * sizeof(double));
  size_t M = N * N;
  fftw_plan plan =
      fftw_plan_r2r_1d(N, fun, fun_hat, FFTW_REDFT10, FFTW_MEASURE);
  double min_time = 1e8;
  double min_time2 = 1e8;
  for (size_t i = 0; i < M; i++) {
    fun[i] = i * i;
  }
  dct2(&fun[0], &fun_hat[0]);
  //dctH2(&fun[0], &fun_hat2[0]);
  double error = 0;
  for (size_t i = 0; i < M; i++) {
    error += fabs(fun_hat2[i] - fun_hat[i]);
    // printf("%2.6e %2.6e \n", fun_hat[i], fun_hat2[i]);
  }
  printf("Error %2.6e\n", error);
  double z[9], zp[9];
  uint64_t t0 = tic;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    ///fftw_execute(plan);
    //dctH2(fun, fun_hat);
    double s1 = toc(g0);
    min_time = s1 < min_time ? s1 : min_time;
  }

  printf("Time avg fftw %2.6e\n", min_time);// toc(t0) / repeat);
  t0 = tic;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    dct2(fun, fun_hat2);
    double s1 = toc(g0);
    min_time2 = s1 < min_time2 ? s1 : min_time2;
  }
  double err = fun_hat[5] - fun_hat2[5];
  printf("Time avg dct2 %2.6e\n", min_time2); // toc(t0) / repeat);
  printf("%lu %2.6e %2.6e \n", N, min_time, err);
  fftw_free(fun);
  fftw_free(fun_hat);
  fftw_free(fun_hat2);
  fftw_destroy_plan(plan);
}

void test_matvectriZ(size_t N, size_t repeat, size_t M, size_t lagrange,
                     size_t verbose) {
  fmm_plan *fmmplan = create_fmm(N, 64, M, 2, lagrange, verbose);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  double min_time = 1e8;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t r0 = tic;
    /*//////////
    for (size_t level = fmmplan->L; level-- > 1;) {
      double *w1 = fmmplan->wk[level - 1];
      for (size_t block = 1; block < get_number_of_blocks(level); block++) {
        size_t Nd = block * 2 * M;
        double *wq = &(fmmplan->wk[level][Nd]);
        int b0 = (block - 1) / 2;
        int q0 = (block - 1) % 2;
        if (lagrange == 0) {
          matvectri(&fmmplan->Th[0], wq, &w1[(b0 * 2 + q0) * M], fmmplan->work,
                     M, false);
        } else {
          cblas_dgemv(CblasRowMajor, CblasTrans, M, M, 1, &fmmplan->Th[0], M,
    &wq[0], 1, 0, &w1[(b0 * 2 + q0) * M], 1); cblas_dgemv(CblasRowMajor,
    CblasTrans, M, M, 1, &fmmplan->Th[M * M], M, &wq[M], 1, 0, &w1[(b0 * 2 + q0)
    * M], 1);
        }
      }
    }
    /*/////////////
    for (size_t level = 0; level < fmmplan->L - 1; level++) {
      double *c0 = fmmplan->ck[level];
      double *c1 = fmmplan->ck[level + 1];
      for (size_t block = 0; block < get_number_of_blocks(level + 1) - 1;
           block++) {
        if (lagrange == 0) {
          matvectri(&fmmplan->ThT[0], &c0[block * M], &c1[block * 2 * M], NULL,
                    M, true);
        } else {
          cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &fmmplan->Th[0], M,
                      &c0[block * M], 1, 1, &c1[block * 2 * M], 1);
          cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &fmmplan->Th[M * M],
                      M, &c0[block * M], 1, 1, &c1[block * 2 * M + M], 1);
        }
      }
    }
    ///////////
    double s1 = toc(r0);
    min_time = s1 < min_time ? s1 : min_time;
  }
  printf("matvectriZ min time %2.6e\n", min_time);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_init(size_t N, size_t maxs, size_t repeat, size_t direction,
               size_t M, size_t lagrange, size_t verbose) {
  if (verbose > 1)
    printf("test_init %lu\n", direction);

  uint64_t t0 = tic;
  double min_time = 1e8;
  size_t flops;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    fmm_plan *fmmplan = create_fmm(N, maxs, M, direction, lagrange, verbose);
    double s1 = toc(g0);
    min_time = s1 < min_time ? s1 : min_time;
    free(fmmplan);
  }
  uint64_t t1 = tic;
  printf("Timing init N %6lu avg / min = %2.4e / %2.4e \n", N,
         dtics(t0, t1) / repeat, min_time);
}


int main(int argc, char *argv[]) {
  int opt;
  size_t N;
  size_t maxs = 64;
  size_t verbose = 2;
  size_t num_threads = 1;
  size_t lagrange = 0;
  size_t M = 18;
  bool R = false;
  double m = 0;
  size_t repeat = 1;
  unsigned direction;
  char *help =
      "Usage: ./l2c options\n"
      "  -h      Show this message\n"
      "  -N      Size of transform\n"
      "  -s      Estimated max size of smallest submatrix  (optional, "
      "default=64)\n"
      "  -M      Chebyshev coefficients (optional, default=18)\n"
      "  -m      Data decay coefficient (optional, default=0). \n"
      "  -r      Repeat computation this many times (for timing, default=1)\n"
      "  -v      Level of verbosity (optional, default=0)\n"
      "  -t      Number of threads if using openmp\n"
      "  -R      Use random data (optional, default=false)\n"
      "  -d      Kind of transform to run\n"
      "       0 - Test speed of Legendre to Chebyshev transform\n"
      "       1 - Test speed of Chebyshev to Legendre transform\n"
      "       2 - Test accuracy of one transform back and forth\n"
      "       3 - Test direct transform back and forth\n";

  while ((opt = getopt(argc, argv, ":N:d:s::M::m::r::v::t::R::l::h")) != -1) {
    switch (opt) {
    case 'N':
      N = atoi(optarg);
      break;
    case 's':
      maxs = atoi(optarg);
      break;
    case 'm':
      m = atof(optarg);
      break;
    case 'r':
      repeat = atoi(optarg);
      break;
    case 'd':
      direction = atoi(optarg);
      break;
    case 'v':
      verbose = atoi(optarg);
      break;
    case 't':
      num_threads = atoi(optarg);
      break;
    case 'M':
      M = atoi(optarg);
      break;
    case 'R':
      R = atoi(optarg);
      break;
    case 'l':
      lagrange = atoi(optarg);
      break;
    case 'h':
      puts(help);
      return 0;
      break;
    default:
      printf("Wrong argument, exiting...");
      exit(-1);
    }
  }
#ifdef OMP
  omp_set_num_threads(num_threads);
#endif
  switch (direction) {
  case 0:
  case 1:
    test_speed(N, maxs, repeat, direction, M, lagrange, verbose);
    break;

  case 2:
    test_forward_backward(N, maxs, M, m, R, lagrange, verbose);
    break;

  case 3:
    test_direct(N, verbose);
    break;

  case 4:
    test_forward_2d(N, N + 2, maxs, verbose, L2C);
    break;

  case 5:
    test_forward_backward_2d(N, N, maxs, verbose);
    break;

  case 6:
    test_direct_speed(N, repeat, L2C, verbose);
    break;

  case 7:
    test_dct(N, repeat);
    break;

  case 8:
    test_lambda();
    break;

  case 9:
    test_accuracy(N, maxs, M, lagrange, verbose);
    break;

  case 10:
    test_init(N, maxs, repeat, L2C, M, lagrange, verbose);
    break;

  default:
    // test_directM(N, repeat, verbose, maxs, M);
    test_matvectriZ(N, repeat, M, lagrange, verbose);
    break;
  }
}
