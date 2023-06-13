#include "leg2cheb.h"
#include <getopt.h>

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
    test_forward_2d(N, N, maxs, verbose, L2C);
    //test_forward_2d(N, 2*N, maxs, verbose, C2L);
    break;

  case 5:
    test_forward_backward_2d(N, N, maxs, verbose);
    break;

  case 6:
    test_direct_speed(N, repeat, direction, verbose);
    break;

  case 0:
  case 1:
    test_speed(N, maxs, repeat, direction, M, verbose);
    break;

  case 7:
    test_dct(N, repeat);
    break;

  default:
    test_directM(N, repeat, verbose);
    break;

  }
}
