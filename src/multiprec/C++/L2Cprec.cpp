#include <getopt.h>
#include <assert.h>
#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>

using boost::multiprecision::number;

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
  T *dn = (T *)malloc(N/2 * sizeof(T));
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
    b[n] = boost::multiprecision::sqrt(pi) * a[n] *
           un[n];
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

int main(int argc, char* argv[]) {
  size_t N;
  double m = 0.0;
  int opt;

    while ((opt = getopt(argc, argv, ":N:m::")) != -1)
    {
       switch (opt)
       {
        case 'N':
          N = atoi(optarg);
          break;
        case 'm':
          m = atof(optarg);
          break;
        default:
          exit(-1);
       }
    }

  //typedef boost::multiprecision::cpp_dec_float_50 T;
  typedef boost::multiprecision::cpp_dec_float_100 T;

  T *u = (T *)malloc(N * sizeof(T));
  T *b = (T *)calloc(N, sizeof(T));
  T *c = (T *)calloc(N, sizeof(T));
  for (size_t i = 0; i < N; i++) {
    u[i] = T(1)/boost::multiprecision::pow(T(i+1), m);
  }
  leg2cheb(u, b, N);
  cheb2leg(b, c, N);
  std::cout << std::scientific
            << std::setprecision(std::numeric_limits<T>::digits10);
  for (size_t i = 0; i < 5; i++) {
    std::cout << c[i] << std::endl;
    assert(boost::multiprecision::abs(c[i]-u[i]) < T(1000*std::numeric_limits<T>::epsilon()));
  }
  free(u);
  free(b);
  free(c);
}