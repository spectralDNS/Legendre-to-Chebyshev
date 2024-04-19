#include "leg2cheb.h"
#include <getopt.h>

size_t min(size_t a, size_t b) { return (a > b) ? b : a; }

size_t max(size_t a, size_t b) { return (a > b) ? a : b; }

double fmax(double a, double b) { return (a > b) ? a : b; }

void test_forward_backward(size_t N, size_t maxs, size_t M, double m,
                           size_t random, size_t lagrange, size_t precompute, size_t verbose) {
  if (verbose > 1)
    printf("test_forward_backward\n");
  fmm_plan *fmmplan = create_fmm(N, maxs, M, BOTH, lagrange, precompute, verbose);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));

  if (verbose > 1)
    printf("Initializing\n");
  // srand48((unsigned int)time(NULL));
  srand48(1);
  // Initialize some input array
  if (random == 1) {
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

  // Compute maximum error
  double error = 0.0;
  for (size_t j = 0; j < N; j++) {
    error = fmax(fabs(input_array[j] - ia[j]), error);
  }

  double e0 = 0;
  for (size_t j = 0; j < N; j++)
    e0 = fmax(fabs(ia[j]), e0);

  double ulp = nextafter(e0, 1e8) - e0;

  printf("N %6lu L inf Error = %2.8e \n", N, error / e0);
  printf("               Flops = %lu\n", flops);
  printf("               ulp   = %2.16e\n", ulp);
  printf("          error ulp  = %2.1f\n", error / ulp);
#ifdef TEST
  assert(error < 1e-10);
#endif
  free(ia);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_accuracy(size_t N, size_t maxs, size_t M, size_t lagrange,
                   size_t precompute, size_t verbose) {
  const double a[5] = {1.7913439415860921e+00, 2.3408718378050253e+00,
                       1.8833606631960720e+00, 1.6587642683880404e+00,
                       1.4417173322290182e+00};
  const double b[5] = {5.2631578947368418e-01, 7.5187969924812026e-02,
                       1.3268465280849182e-01, 1.7768205680441512e-01,
                       2.4367824933176932e-01};

  assert(N < 1000);
  if (verbose > 1)
    printf("test_accuracy N=20 direct and N=%d FMM\n", N);
  fmm_plan *fmmplan = create_fmm(20, maxs, M, BOTH, lagrange, precompute, 1);

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

  fmm_plan *fmmplan2 = create_fmm(N, maxs, M, BOTH, lagrange, precompute, 1);
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
  fmm_plan *fmmplan3 = create_fmm(N, maxs, M, L2C, lagrange, precompute, 1);
  direct(input_array, output_array, fmmplan3->dplan, L2C, 1);

  for (size_t i = 0; i < N; i++)
    out[i] = 0.0;

  execute(input_array, out, fmmplan3, L2C, 1);
  error0 = 0.0;
  for (size_t i = 0; i < N; i++) {
    error0 += pow(output_array[i] - out[i], 2);
  }
  // assert(sqrt(error0) < 1e-7);

  fmm_plan *fmmplan4 = create_fmm(N, maxs, M, C2L, lagrange, precompute, 1);
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
                size_t M, size_t lagrange, size_t precompute, size_t verbose) {
  if (verbose > 1)
    printf("test_speed %lu\n", direction);
  fmm_plan *fmmplan = create_fmm(N, maxs, M, direction, lagrange, precompute, verbose);

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

  uint64_t t1 = tic;
  printf("Timing N %6lu avg / min = %2.4e / %2.4e flops = %lu\n", N,
         dtics(t0, t1) / repeat, min_time, flops);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_direct_speed(size_t N, size_t repeat, size_t direction,
                       size_t precompute, size_t verbose) {
  if (verbose > 1)
    printf("test_direct_speed %lu\n", direction);
  direct_plan *dplan = create_direct(N, 2, precompute);

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

void test_direct(size_t N, size_t precompute, size_t verbose) {
  direct_plan *dplan = create_direct(N, 2, precompute);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // Initialize some input array
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
  fmm_plan *fmmplan = create_fmm(N, maxs, 18, 2, 0, 0, verbose);
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
  fmmplan = create_fmm(2 * N, maxs, 18, 2, 0, 0, verbose);
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
    fmm_plan *fmmplan = create_fmm(N, maxs, 18, 2, 0, 0, verbose);
    for (size_t i = 0; i < N; i++) {
      input_array1d[i] = 1.0 + i;
      output_array1d[i] = 0.0;
    }
    size_t flops =
        execute(input_array1d, output_array1d, fmmplan, direction, 1);

    fmm_plan_2d *fmmplan2d =
        create_fmm_2d(N0, N1, axis, maxs, 18, 2, 0, 0, verbose);
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

  srand48(1);

  int axis = -1;
  fmm_plan_2d *fmmplan2d = create_fmm_2d(N0, N1, axis, maxs, 18, 2, 0, 0, verbose);
  // Initialize some input array
  for (size_t i = 0; i < N0; i++) {
    for (size_t j = 0; j < N1; j++) {
      input_array[i * N1 + j] = 2 * drand48() - 1;
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

void test_directM(size_t N, size_t repeat, size_t verbose, size_t s, size_t M, size_t precompute) {
  fmm_plan *fmmplan = create_fmm(N, s, M, 2, 0, precompute, verbose);
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

void test_dct0(size_t N, size_t repeat) {
  double *fun = (double *)fftw_malloc(N * sizeof(double));
  double *fun_hat = (double *)fftw_malloc(N * sizeof(double));

  fftw_plan plan =
      fftw_plan_r2r_1d(N, fun, fun_hat, FFTW_REDFT10, FFTW_PATIENT);

  double min_time = 1e8;
  for (size_t i = 0; i < N; i++) {
    fun[i] = i * i;
  }

  uint64_t t0 = tic;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    fftw_execute(plan);
    double s1 = toc(g0);
    min_time = s1 < min_time ? s1 : min_time;
  }
  printf("Time %2.6e %2.6e N = %d\n", min_time, toc(t0) / repeat, N);
  fftw_free(fun);
  fftw_free(fun_hat);
  fftw_destroy_plan(plan);
  return;
}

void test_dct(size_t N, size_t repeat) {
  double *fun = (double *)fftw_malloc(N * sizeof(double));
  double *fun_hat = (double *)fftw_malloc(N * sizeof(double));

  fftw_plan plan =
      fftw_plan_r2r_1d(N, fun, fun_hat, FFTW_REDFT10, FFTW_PATIENT);

  double min_time = 1e8;
  double min_time2 = 1e8;
  for (size_t i = 0; i < N; i++) {
    fun[i] = i;
  }

  uint64_t t0 = tic;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    fftw_execute(plan);
    double s1 = toc(g0);
    min_time = s1 < min_time ? s1 : min_time;
  }
  printf("Time avg fftw %2.6e %2.6e N = %d\n", min_time, toc(t0) / repeat, N);

  double add, mul, fma;
  fftw_flops(plan, &add, &mul, &fma);
  printf("N %d Flops %f\n", N, add + mul + 2 * fma);

  double fun2[18], fun3[324];
  double fun_hat2[18], fun_hat3[324];
  for (size_t i = 0; i < 18; i++) {
    fun2[i] = i * i;
    for (size_t j = 0; j < 18; j++) {
      fun3[i * 18 + j] = i + j;
    }
  }
  dct(fun2, fun_hat, 1);
  dct_radix23(fun2, fun_hat2, 1);
  double err0 = 0;
  for (size_t i = 0; i < 18; i++) {
    err0 += pow(fun_hat[i] - fun_hat2[i], 2);
  }
  printf("Err dct_fft %2.16e\n", err0);

  t0 = tic;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    dct(fun2, fun_hat2, 1);
    double s1 = toc(g0);
    min_time2 = s1 < min_time2 ? s1 : min_time2;
  }
  printf("Time avg dct     %2.6e %2.6e\n", min_time2, toc(t0) / repeat);

  t0 = tic;
  min_time2 = 1e8;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    dct_fft(fun2, fun_hat2, 1);
    double s1 = toc(g0);
    min_time2 = s1 < min_time2 ? s1 : min_time2;
  }
  printf("Time avg dct_fft %2.6e %2.6e\n", min_time2, toc(t0) / repeat);

  t0 = tic;
  min_time2 = 1e8;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    dct_radix2(fun2, fun_hat2, 1);
    double s1 = toc(g0);
    min_time2 = s1 < min_time2 ? s1 : min_time2;
  }
  printf("Time avg dct_r22 %2.6e %2.6e\n", min_time2, toc(t0) / repeat);

  t0 = tic;
  min_time2 = 1e8;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    dct_radix23(fun2, fun_hat2, 1);
    double s1 = toc(g0);
    min_time2 = s1 < min_time2 ? s1 : min_time2;
  }
  printf("Time avg dct_r23 %2.6e %2.6e\n", min_time2, toc(t0) / repeat);

  t0 = tic;
  min_time2 = 1e8;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    dct0(fun2, fun_hat2, 1);
    double s1 = toc(g0);
    min_time2 = s1 < min_time2 ? s1 : min_time2;
  }
  printf("Time avg dct0    %2.6e %2.6e\n", min_time2, toc(t0) / repeat);

  double v[18];
  for (size_t i = 0; i < 18; i++) {
    v[i] = i * i;
  }

  fftw_complex x[9];
  fftw_complex *xc = (fftw_complex *)fftw_malloc(9 * sizeof(fftw_complex));
  fftw_complex *xh = (fftw_complex *)fftw_malloc(9 * sizeof(fftw_complex));
  fftw_plan planc = fftw_plan_dft_1d(9, xc, xh, FFTW_FORWARD, FFTW_PATIENT);

  for (int i = 0; i < 9; ++i) {
    x[i] = (fftw_complex){v[2 * i], v[2 * i + 1]};
    xc[i] = (fftw_complex){v[2 * i], v[2 * i + 1]};
  }

  fft9c(x);
  fftw_execute(planc);

  double err = 0;
  for (size_t i = 0; i < 9; i++) {
    err +=
        pow(creal(xh[i]) - creal(x[i]), 2) + pow(cimag(xh[i]) - cimag(x[i]), 2);
  }
  printf("Error fft9 %2.6e %2.6e \n", err);

  t0 = tic;
  min_time2 = 1e8;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    fft9(x);
    double s1 = toc(g0);
    min_time2 = s1 < min_time2 ? s1 : min_time2;
  }
  printf("Time avg fft9  fftw_complex %2.6e %2.6e\n", min_time2,
         toc(t0) / repeat);

  t0 = tic;
  min_time2 = 1e8;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    fft9c(x);
    double s1 = toc(g0);
    min_time2 = s1 < min_time2 ? s1 : min_time2;
  }
  printf("Time avg fft9c fftw_complex %2.6e %2.6e\n", min_time2,
         toc(t0) / repeat);

  t0 = tic;
  min_time2 = 1e8;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    fftw_execute(planc);
    double s1 = toc(g0);
    min_time2 = s1 < min_time2 ? s1 : min_time2;
  }
  printf("Time avg fftw               %2.6e %2.6e\n", min_time2,
         toc(t0) / repeat);

  fftw_print_plan(planc);
  fftw_print_plan(plan);
  fftw_destroy_plan(planc);
  fftw_destroy_plan(plan);
  fftw_free(xc);
  fftw_free(xh);
  fftw_free(fun);
  fftw_free(fun_hat);
}

void test_flops_dct_fftw() {
  double add, mul, fma;
  for (size_t i = 6; i < 12; i++) {
    size_t N = pow(2, i);
    double *fun = (double *)fftw_malloc(N * sizeof(double));
    double *fun_hat = (double *)fftw_malloc(N * sizeof(double));
    fftw_plan plan =
        fftw_plan_r2r_1d(N, fun, fun_hat, FFTW_REDFT10, FFTW_MEASURE);

    fftw_flops(plan, &add, &mul, &fma);
    printf("%d %8.2f %8.2f \n", N, (add + mul + 2 * fma), 2 * N * log2(N));
  }
}

void test_matvectri(size_t N, size_t repeat, size_t M, size_t lagrange,
                     size_t precompute, size_t verbose) {
  fmm_plan *fmmplan = create_fmm(N, 64, M, 2, lagrange, precompute, verbose);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  double min_time = 1e8;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t r0 = tic;
    for (size_t level = 0; level < fmmplan->L - 1; level++) {
      double *c0 = fmmplan->ck[level];
      double *c1 = fmmplan->ck[level + 1];
      for (size_t block = 0; block < get_number_of_blocks(level + 1) - 1;
           block++) {
        if (lagrange == 0) {
          matvectri(&fmmplan->BT[0], &c0[block * M], &c1[block * 2 * M], NULL,
                    M, true);
        } else {
          cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &fmmplan->B[0], M,
                      &c0[block * M], 1, 1, &c1[block * 2 * M], 1);
          cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, &fmmplan->B[M * M],
                      M, &c0[block * M], 1, 1, &c1[block * 2 * M + M], 1);
        }
      }
    }
    double s1 = toc(r0);
    min_time = s1 < min_time ? s1 : min_time;
  }
  printf("matvectri min time %2.6e\n", min_time);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_init(size_t N, size_t maxs, size_t repeat, size_t direction, size_t M,
               size_t lagrange, size_t precompute, size_t verbose) {
  if (verbose > 1)
    printf("test_init %lu\n", direction);

  uint64_t t0 = tic;
  double min_time = 1e8;
  size_t flops;
  for (size_t i = 0; i < repeat; i++) {
    uint64_t g0 = tic;
    fmm_plan *fmmplan = create_fmm(N, maxs, M, direction, lagrange, precompute, verbose);
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
  size_t R = 0;
  double m = 0;
  size_t precompute = 0;
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
      "  -t      Number of threads if using openmp (optional, default=1)\n"
      "  -p      Precompute direct part of matrix (optional, default=0)\n"
      "  -R      Use random data (optional, default=0)\n"
      "  -d      Kind of transform to run\n"
      "       0 - Test speed of Legendre to Chebyshev transform\n"
      "       1 - Test speed of Chebyshev to Legendre transform\n"
      "       2 - Test accuracy of one transform back and forth\n"
      "       3 - Test direct transform back and forth\n";

  while ((opt = getopt(argc, argv, ":N:d:s::M::m::r::v::t::p::R::l::h")) != -1) {
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
    case 'p':
      precompute = atoi(optarg);
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
    test_speed(N, maxs, repeat, direction, M, lagrange, precompute, verbose);
    break;

  case 2:
    test_forward_backward(N, maxs, M, m, R, lagrange, precompute, verbose);
    break;

  case 3:
    test_direct(N, precompute, verbose);
    break;

  case 4:
    test_forward_2d(N, N + 2, maxs, verbose, L2C);
    break;

  case 5:
    test_forward_backward_2d(N, N, maxs, verbose);
    break;

  case 6:
    test_direct_speed(N, repeat, L2C, precompute, verbose);
    break;

  case 7:
    test_dct0(N, repeat);
    break;

  case 8:
    test_accuracy(N, maxs, M, lagrange, precompute, verbose);
    break;

  case 9:
    test_init(N, maxs, repeat, L2C, M, lagrange, precompute, verbose);
    break;

  case 10:
    test_flops_dct_fftw();
    break;

  default:
    test_matvectri(N, repeat, M, lagrange, precompute, verbose);
    break;
  }
}
