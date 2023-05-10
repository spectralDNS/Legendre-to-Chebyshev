#include <getopt.h>
#include "leg2cheb.h"

int main(int argc, char *argv[])
{
    int opt;
    size_t N;
    size_t maxs = 36;
    size_t verbose = 2;
    double m = 0;
    size_t repeat = 1;
    unsigned direction;
    char *help = "Usage: ./l2c options\n"
                 "  -h      Show this message\n"
                 "  -N      Size of transform\n"
                 "  -s      Estimated max size of smallest submatrix  (optional, default=36)\n"
                 "  -m      Data decay coefficient (optional, default=0). \n"
                 "  -r      Repeat computation this many times (for timing, default=1)\n"
                 "  -v      Level of verbosity (optional, default=0)\n"
                 "  -d      Kind of transform to run\n"
                 "       0 - Test speed of Legendre to Chebyshev transform\n"
                 "       1 - Test speed of Chebyshev to Legendre transform\n"
                 "       2 - Test accuracy of one transform back and forth\n"
                 "       3 - Test direct transform back and forth\n";


    while ((opt = getopt(argc, argv, ":N:d:s::m::r::v::h")) != -1)
    {
       switch (opt)
       {
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
    omp_set_num_threads(1);
#endif
    switch (direction)
    {
    case 2:
      test_foreward_backward(N, maxs, m, verbose);
      break;

    case 3:
      test_direct(N, verbose);
      break;

    default:
      test_speed(N, maxs, repeat, direction, verbose);
      break;
    }
}
