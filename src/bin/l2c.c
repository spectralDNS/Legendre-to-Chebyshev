#include "leg2cheb.h"
#include <getopt.h>

void test_forward_backward(size_t N, size_t maxs, size_t M, double m,
                           bool random, size_t verbose)
{
  if (verbose > 1)
    printf("test_forward_backward\n");
  fmm_plan *fmmplan = create_fmm(N, maxs, M, 2, verbose);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // srand48((unsigned int)time(NULL));
  srand48(1);
  // Initialize some input array
  switch (random)
  {
  case true:
  {
    for (size_t i = 0; i < N; i++)
      input_array[i] = (2 * drand48() - 1) / pow(i + 1, m);
  }
  break;
  case false:
  {
    for (size_t i = 0; i < N; i++)
      input_array[i] = 1.0 / pow(i + 1, m);
  }
  break;
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

void test_speed(size_t N, size_t maxs, size_t repeat, size_t direction,
                size_t M, size_t verbose)
{
  if (verbose > 1)
    printf("test_speed %lu\n", direction);
  fmm_plan *fmmplan = create_fmm(N, maxs, M, direction, verbose);

  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // Initialize some input array
  for (size_t i = 0; i < N; i++)
    input_array[i] = 1.0;

  uint64_t t0 = tic;
  double min_time = 1e8;
  size_t flops;
  for (size_t i = 0; i < repeat; i++)
  {
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
                       size_t verbose)
{
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
  for (size_t i = 0; i < repeat; i++)
  {
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

void test_direct(size_t N, size_t verbose)
{
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
  for (size_t j = 0; j < N; j++)
  {
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

void test_2_sizes(size_t N, size_t maxs, size_t verbose)
{
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

void test_forward_2d(size_t N0, size_t N1, size_t maxs, size_t verbose,
                     size_t direction)
{
  if (verbose > 1)
  {
    printf("test_forward_2d\n");
  }

  double *input_array = (double *)calloc(N0 * N1, sizeof(double));
  double *output_array = (double *)calloc(N0 * N1, sizeof(double));
  double *input_array1d = (double *)calloc(max(N0, N1), sizeof(double));
  double *output_array1d = (double *)calloc(max(N0, N1), sizeof(double));

  for (size_t axis = 0; axis < 2; axis++)
  {
    // 1D first
    size_t N = ((axis == 0) ? N0 : N1);
    fmm_plan *fmmplan = create_fmm(N, maxs, 18, 2, verbose);
    for (size_t i = 0; i < N; i++)
    {
      input_array1d[i] = 1.0 + i;
      output_array1d[i] = 0.0;
    }
    size_t flops =
        execute(input_array1d, output_array1d, fmmplan, direction, 1);

    fmm_plan_2d *fmmplan2d = create_fmm_2d(N0, N1, axis, maxs, 18, 2, verbose);
    // Initialize some input array
    for (size_t i = 0; i < N0; i++)
    {
      for (size_t j = 0; j < N1; j++)
      {
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
                              size_t verbose)
{
  if (verbose > 1)
  {
    printf("test_forward_backward_2d\n");
  }

  double *input_array = (double *)calloc(N0 * N1, sizeof(double));
  double *output_array = (double *)calloc(N0 * N1, sizeof(double));
  double *output_array2 = (double *)calloc(N0 * N1, sizeof(double));

  int axis = -1;
  fmm_plan_2d *fmmplan2d = create_fmm_2d(N0, N1, axis, maxs, 18, 2, verbose);
  // Initialize some input array
  for (size_t i = 0; i < N0; i++)
  {
    for (size_t j = 0; j < N1; j++)
    {
      input_array[i * N1 + j] = 1.0;
    }
  }

  size_t flops = execute2D(input_array, output_array, fmmplan2d, L2C);
  flops = execute2D(output_array, output_array2, fmmplan2d, C2L);

  double error = 0.0;
  for (size_t i = 0; i < N0; i++)
  {
    for (size_t j = 0; j < N1; j++)
    {
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

void test_directM(size_t N, size_t repeat, size_t verbose)
{
  fmm_plan *fmmplan = create_fmm(N, 64, 18, 2, verbose);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  size_t flops;
  double min_time = 1e8;
  for (size_t i = 0; i < repeat; i++)
  {
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

void test_dct(size_t N, size_t repeat)
{

  double *fun = (double *)fftw_malloc(N * sizeof(double));
  double *fun_hat = (double *)fftw_malloc(N * sizeof(double));
  fftw_plan plan1d =
      fftw_plan_r2r_1d(N, fun, fun_hat, FFTW_REDFT10, FFTW_MEASURE);
  double min_time = 1e8;
  uint64_t t0 = tic;
  for (size_t i = 0; i < repeat; i++)
  {
    uint64_t g0 = tic;
    fftw_execute(plan1d);
    double s1 = toc(g0);
    min_time = s1 < min_time ? s1 : min_time;
  }
  printf("%lu %2.6e %2.6e\n", N, min_time, toc(t0) / repeat);
  fftw_free(fun);
  fftw_free(fun_hat);
  fftw_destroy_plan(plan1d);
}

void test_matvectriZ(size_t N, size_t repeat, size_t verbose)
{
  fmm_plan *fmmplan = create_fmm(N, 36, 18, 2, verbose);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  size_t flops;
  double min_time = 1e8;
  for (size_t i = 0; i < repeat; i++)
  {
    uint64_t r0 = tic;
    for (size_t level = fmmplan->L; level-- > 1;)
    {
      double *w0 = fmmplan->wk[level];
      double *w1 = fmmplan->wk[level - 1];
      size_t M = fmmplan->M;
      for (size_t block = 1; block < get_number_of_blocks(level); block++)
      {
        size_t Nd = block * 2 * M;
        double *wq = &(fmmplan->wk[level][Nd]);
        int b0 = (block - 1) / 2;
        int q0 = (block - 1) % 2;
        matvectriZ(&fmmplan->Th[0], wq, &w1[(b0 * 2 + q0) * M], fmmplan->work,
                   M, M, M, true);
      }
    }
    double s1 = toc(r0);
    min_time = s1 < min_time ? s1 : min_time;
  }
  printf("matvectriZ min time %2.6e\n", min_time);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}


int main(int argc, char *argv[]) {
  int opt;
  size_t N;
  size_t maxs = 64;
  size_t verbose = 2;
  size_t num_threads = 1;
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
      "  -d      Kind of transform to run\n"
      "       0 - Test speed of Legendre to Chebyshev transform\n"
      "       1 - Test speed of Chebyshev to Legendre transform\n"
      "       2 - Test accuracy of one transform back and forth\n"
      "       3 - Test direct transform back and forth\n";

  while ((opt = getopt(argc, argv, ":N:d:s::M::m::r::v::t::R::h")) != -1) {
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
  case 2:
    test_forward_backward(N, maxs, M, m, R, verbose);
    break;

  case 3:
    test_direct(N, verbose);
    break;

  case 4:
    test_forward_2d(N, N+2, maxs, verbose, L2C);
    break;

  case 5:
    test_forward_backward_2d(N, N, maxs, verbose);
    break;

  case 6:
    test_direct_speed(N, repeat, L2C, verbose);
    break;

  case 0:
  case 1:
    test_speed(N, maxs, repeat, direction, M, verbose);
    break;

  case 7:
    test_dct(N, repeat);
    break;

  case 8:
    test_matvectriZ(N, repeat, verbose);
    break;

  default:
    test_directM(N, repeat, verbose);
    break;

  }
}
