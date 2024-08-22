#include <assert.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <getopt.h>
#include <iostream>

extern "C" {
#include "leg2cheb.h"
}

using boost::multiprecision::cpp_dec_float;
using boost::multiprecision::number;

typedef number<cpp_dec_float<32>> cpp_dec_float_32;

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

void test_accuracy_C(size_t N, double m, size_t direction, size_t norm,
                     size_t random) {

  // typedef boost::multiprecision::cpp_dec_float_100 T;
  // typedef boost::multiprecision::cpp_dec_float_50 T;
  typedef cpp_dec_float_32 T;

  //srand(time(NULL));   // Initialization, should only be called once.
  srand(1);
  
  T *u = (T *)malloc(N * sizeof(T));
  T *b = (T *)calloc(N, sizeof(T));
  switch (random) {
  case 0:
    for (size_t i = 0; i < N; i++)
      u[i] = T(1) / (boost::multiprecision::pow(T(i + 1), m));
    break;

  case 1:
    for (size_t i = 0; i < N; i++)
      u[i] =  (T(rand()) / T(RAND_MAX)) / boost::multiprecision::pow(T(i + 1), m);
      //u[i] = (2 * T(rand()) / T(RAND_MAX) - 1) /
      //       boost::multiprecision::pow(T(i + 1), m);
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

  double *input_array = (double *)calloc(N, sizeof(double));
  double *output_arrayM = (double *)calloc(N, sizeof(double));
  double *output_arrayN = (double *)calloc(N, sizeof(double));

  fmm_plan *fmmplanM = create_fmm(N, 32, 18, direction, 0, 0, 1);
  fmm_plan *fmmplanN = create_fmm(N, 32, 18, direction, 1, 0, 1);

  for (size_t i = 0; i < N; i++)
    input_array[i] = (double)u[i];

  size_t flopsM = execute(input_array, output_arrayM, fmmplanM, direction, 1);
  size_t flopsN = execute(input_array, output_arrayN, fmmplanN, direction, 1);

  double errorM = 0;
  double errorN = 0;
  double max_output;
  double ulp;
  max_output = 0;
  for (size_t i = 0; i < N; i++) {
    max_output = fmax((double)b[i], max_output);
  }
  //std::cout << max_output << std::endl;
  ulp = nextafter(max_output, 1e8) - max_output;
  switch (norm) {
  case 0: // L2 norm
  {
    for (size_t i = 0; i < N; i++) {
      errorM += pow(output_arrayM[i] - (double)b[i], 2);
      errorN += pow(output_arrayN[i] - (double)b[i], 2);
    }
    errorM = sqrt(errorM);
    errorN = sqrt(errorN);
  } break;

  case 1: // inf norm
  {

    for (size_t i = 0; i < N; i++) {
      double e0 = fabs(output_arrayM[i] - (double)b[i]);
      double e1 = fabs(output_arrayN[i] - (double)b[i]);
      if (e0 > errorM) {
        errorM = e0;
        //std::cout << "Max for i=" <<  i << std::endl;;
      }
      if (e1 > errorN) {
        errorN = e1;
        //std::cout << "Max for i=" <<  i << std::endl;;
      }

      //errorM = fmax(fabs(output_arrayM[i] - (double)b[i]), errorM);
      //errorN = fmax(fabs(output_arrayN[i] - (double)b[i]), errorN);
      // std::cout << std::setprecision(16) << errorM << " " <<
      // fabs(output_arrayM[i] - (double)b[i]) << " " << (double)b[i] <<  " " <<
      // ulp << " " << errorM/ulp << std::endl;
    }
    errorM /= max_output;
    errorN /= max_output;
  }
  default:
    break;
  }
  //for (size_t i = 0; i < 20; i++)
  //{
  //  std::cout << (double)b[i] << std::endl;
  //}

  free(input_array);
  free(output_arrayM);
  free(output_arrayN);
  free(u);
  free(b);

  std::cout << std::setprecision(16) << errorM * max_output << " "
            << errorM * max_output / ulp << " " << errorN * max_output << " "
            << errorN * max_output / ulp << " " << ulp << std::endl;
  // std::cout << nextafter(1.0, 1e8)-1.0 << " " << nextafter(2.0, 1e8) -2.0 <<
  // " " << nextafter(4.0, 1e8) - 4.0 << std::endl;
  return;
}

void test_accuracy(size_t N, size_t m) {
  // typedef boost::multiprecision::cpp_dec_float_50 T;
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
  switch (a) {
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