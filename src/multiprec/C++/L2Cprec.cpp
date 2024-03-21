#include <assert.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <getopt.h>
#include <iostream>

extern "C" {
#include "leg2cheb.h"
}

using boost::multiprecision::number;
using boost::multiprecision::cpp_dec_float;

typedef number<cpp_dec_float<32> > cpp_dec_float_32;

template <class T> T Lambda(T z) {
  return boost::math::tgamma(z + T(0.5)) / boost::math::tgamma(z + 1);
}

template <class T> void leg2cheb(T *u, T *b, size_t N) {
  T *a = (T *)malloc(N * sizeof(T));
  T pi = boost::math::constants::pi<T>();

  for (size_t i = 0; i < N; i++)
    a[i] = Lambda(T(i));

  for (size_t n = 0; n < N; n = n + 2) {
    for (size_t i = 0; i < N - n; i++) {
      b[i] += a[n / 2] * a[n / 2 + i] * u[n + i];
    }
  }
  b[0] /= T(2);
  for (size_t i = 0; i < N; i++) {
    b[i] *= T(2 / pi);
  }
  free(a);
}

template <class T> void cheb2leg(T *u, T *b, size_t N) {
  T *a = (T *)malloc(N * sizeof(T));
  T *dn = (T *)malloc(N / 2 * sizeof(T));
  T *un = (T *)malloc(N * sizeof(T));
  T pi = boost::math::constants::pi<T>();

  dn[0] = 0;
  for (size_t i = 1; i < N / 2; i++) {
    dn[i] = Lambda(T(i - 1)) / (2 * i);
  }
  a[0] = 2 / T(boost::multiprecision::sqrt(pi));
  for (size_t i = 1; i < N; i++)
    a[i] = 1 / (2 * Lambda(T(i)) * i * (i + T(0.5)));
  un[0] = u[0];
  un[1] = u[1];
  for (size_t i = 2; i < N; i++)
    un[i] = u[i] * i;
  for (size_t n = 0; n < N; n++) {
    b[n] = boost::multiprecision::sqrt(pi) * a[n] * un[n];
  }
  for (size_t n = 2; n < N; n = n + 2) {
    for (size_t i = 0; i < N - n; i++) {
      b[i] -= dn[n / 2] * a[n / 2 + i] * un[n + i];
    }
  }
  for (size_t n = 0; n < N; n++)
    b[n] *= T(n + 0.5);

  free(a);
  free(dn);
  free(un);
}

double test_accuracy_C(size_t N, size_t m, size_t direction, size_t norm, size_t random) {

  //typedef boost::multiprecision::cpp_dec_float_100 T;
  //typedef boost::multiprecision::cpp_dec_float_50 T;
  typedef cpp_dec_float_32 T;
  srand(time(NULL));   // Initialization, should only be called once.
  T *u = (T *)malloc(N * sizeof(T));
  T *b = (T *)calloc(N, sizeof(T));
  switch (random)
  {
  case 0:
    for (size_t i = 0; i < N; i++)
      u[i] = T(1) / (boost::multiprecision::pow(T(i + 1), m));
    break;

  case 1:
    for (size_t i = 0; i < N; i++)
      u[i] =  (T(rand()) / T(RAND_MAX)) / boost::multiprecision::pow(T(i + 1), m);
      //u[i] =  (2 * T(rand()) / T(RAND_MAX) - 1) / boost::multiprecision::pow(T(i + 1), m);
    break;

  }

  switch (direction) {
  case 0:
    leg2cheb(u, b, N);
    break;

  default:
    cheb2leg(u, b, N);
    break;
  }

  fmm_plan *fmmplan = create_fmm(N, 64, 18, direction, 0, 1);
  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_array = (double *)calloc(N, sizeof(double));
  for (size_t i = 0; i < N; i++)
    input_array[i] = (double)u[i];
  size_t flops = execute(input_array, output_array, fmmplan, direction, 1);

  double error = 0;
  switch (norm)
  {
  case 0: // L2 norm
    {
      for (size_t i = 0; i < N; i++)
        error += pow(output_array[i] - (double)b[i], 2);
      error = sqrt(error);
    }
    break;

  case 1: // inf norm
    {
      for (size_t i = 0; i < N; i++)
        error = fmax(fabs(output_array[i] - (double)b[i]), error);
      double e0 = 0;
      for (size_t i = 0; i < N; i++)
        e0 = fmax(fabs(output_array[i]), e0);
      error /= e0;
    }
  default:
    break;
  }
  free(input_array);
  free(output_array);
  free(u);
  free(b);
  std::cout << error << std::endl;

  return error;
}

void test_accuracy(size_t N, size_t m) {
  //typedef boost::multiprecision::cpp_dec_float_50 T;
  typedef boost::multiprecision::cpp_dec_float_100 T;
  T *u = (T *)malloc(N * sizeof(T));
  T *b = (T *)calloc(N, sizeof(T));
  T *c = (T *)calloc(N, sizeof(T));
  for (size_t i = 0; i < N; i++) {
    u[i] = T(1) / boost::multiprecision::pow(T(i + 1), m);
  }
  leg2cheb(u, b, N);
  cheb2leg(b, c, N);
  std::cout << std::scientific
            << std::setprecision(std::numeric_limits<T>::digits10);
  for (size_t i = 0; i < 5; i++) {
    std::cout << c[i] << std::endl;
    assert(boost::multiprecision::abs(c[i] - u[i]) <
           T(1000 * std::numeric_limits<T>::epsilon()));
  }
  free(u);
  free(b);
  free(c);
}

int main(int argc, char *argv[]) {
  size_t N;
  size_t a = 0;
  size_t d = 0;
  size_t n = 0;
  size_t R = 0;
  double m = 0.0;
  int opt;
  while ((opt = getopt(argc, argv, ":N:m::a::d::n::R::")) != -1) {
    switch (opt) {
    case 'N':
      N = atoi(optarg);
      break;
    case 'm':
      m = atof(optarg);
      break;
    case 'a':
      a = atoi(optarg);
      break;
    case 'd':
      d = atoi(optarg);
      break;
    case 'n':
      n = atoi(optarg);
      break;
    case 'R':
      R = atoi(optarg);
      break;
    default:
      exit(-1);
    }
  }
  switch (a)
  {
  case 0:
    test_accuracy(N, m);
    break;

  case 1:
    test_accuracy_C(N, m, d, n, R);
    break;

  default:
    break;
  }
}