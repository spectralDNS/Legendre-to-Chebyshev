#include "leg2cheb.h"

void test_foreward_backward(size_t N, size_t maxs) {
  fmm_plan *fmmplan = create_fmm(N, maxs, 2);
  double *input_array = (double *)calloc(fmmplan->Nn, sizeof(double));
  double *output_array = (double *)calloc(fmmplan->Nn, sizeof(double));
  // Initialize some input array
  for (size_t i = 0; i < N; i++)
    input_array[i] = 1.0;
  // Leg to Cheb
  size_t flops = execute(input_array, output_array, fmmplan, 0);
  double *ia = (double *)calloc(fmmplan->Nn, sizeof(double));
  // Cheb to Leg
  flops += execute(output_array, ia, fmmplan, 1);
  // Compute L2 error norm
  double error = 0.0;
  for (size_t j = 0; j < N; j++)
    error += pow(input_array[j] - ia[j], 2);
  printf("L2 Error = %2.4e \n", sqrt(error));
  printf("Flops %lu\n\n", flops);
  assert(error < 1e-10);
  free(ia);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_speed(size_t N, size_t maxs, size_t repeat, unsigned direction) {
  fmm_plan *fmmplan = create_fmm(N, maxs, direction);
  double *input_array = (double *)calloc(fmmplan->Nn, sizeof(double));
  double *output_array = (double *)calloc(fmmplan->Nn, sizeof(double));
  // Initialize some input array
  for (size_t i = 0; i < N; i++)
    input_array[i] = 1.0;

  struct timeval t0, t1;
  gettimeofday(&t0, 0);
  double min_time = 1e8;
  for (size_t i = 0; i < repeat; i++) {
    struct timeval r0, r1;
    gettimeofday(&r0, 0);
    for (size_t j = 0; j < N; j++)
      output_array[j] = 0.0;
    execute(input_array, output_array, fmmplan, direction);
    gettimeofday(&r1, 0);
    double s1 = tdiff_sec(r0, r1);
    min_time = s1 < min_time ? s1 : min_time;
  }
  gettimeofday(&t1, 0);
  printf("Timing avg / min = %2.4e / %2.4e \n", tdiff_sec(t0, t1) / repeat,
         min_time);
  free(input_array);
  free(output_array);
  free_fmm(fmmplan);
}

void test_direct(size_t N) {
  direct_plan *dplan = create_direct(N, 2);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  // Initialize some input array
  for (size_t i = 0; i < N; i++)
    input_array[i] = 1.0;

  size_t flops;
  flops = direct(input_array, output_array, dplan, 0);
  double *ia = (double *)calloc(N, sizeof(double));
  // Cheb to Leg
  flops += direct(output_array, ia, dplan, 1);
  // Compute L2 error norm
  double error = 0.0;
  for (size_t j = 0; j < N; j++) {
    error += pow(input_array[j] - ia[j], 2);
  }
  printf("L2 Error direct = %2.4e \n", sqrt(error));
  assert(sqrt(error) < 1e-10);
  free(input_array);
  free(output_array);
  free_direct(dplan);
}

int main(int argc, char *argv[]) {
  test_foreward_backward(1000, 100);
  test_foreward_backward(1000, 36);
  test_foreward_backward(10000, 100);
  test_foreward_backward(10000, 36);
  test_speed(100, 10, 100, 0);
  for (size_t i = 4; i < 12; i++)
    test_speed(36*pow(2, i+2), 36, 10, 0);
  test_direct(10);
  test_direct(200);
  test_direct(1000);
}
