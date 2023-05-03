#include <getopt.h>
#include <string.h>
#include "leg2cheb.h"

int main(int argc, char *argv[])
{
    int opt;
    size_t N;
    size_t maxs = 36;
    size_t v = 0;
    double m = 0;
    size_t repeat = 1;
    unsigned direction = 0;
    char *help = "Usage: ./l2c options\n"
                 "  -h      Show this message\n"
                 "  -N      Size of transform\n"
                 "  -s      Estimated max size of smallest submatrix  (optional, default=36)\n"
                 "  -m      Data decay coefficient. \n"
                 "  -r      Repeat computation this many times (for timing, default=1)\n"
                 "  -v      Level of verbosity\n"
                 "  -d      Direction of transform. 0 for Legendre to Chebyshev,\n"
                 "          1 for Chebyshev to Legendre and 2 use both and check accuracy\n";

    while ((opt = getopt(argc, argv, "N:d:s::m::r::v:h")) != -1)
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
        case 'r':
          repeat = atoi(optarg);
          break;
        case 'd':
          direction = atoi(optarg);
          break;
        case 'h':
          puts(help);
          return 0;
          break;
        default:
          printf("Exiting...");
          exit(-1);
       }
    }
    fmm_plan* fmmplan = create_fmm(N, maxs, direction, v);
    double* input_array=(double *)calloc(fmmplan->Nn, sizeof(double));
    double* output_array=(double *)calloc(fmmplan->Nn, sizeof(double));
    for (size_t i = 0; i < N; i++)
      input_array[i] = 1.0/pow(i, m);

    struct timeval t0, t1;
    double min_time = 1e8;
    size_t flops;
    gettimeofday(&t0, 0);
    if (direction == 2)
    {
        double error = 0.0;
        gettimeofday(&t0, 0);
        for (size_t j = 0; j < N; j++)
            output_array[j] = 0.0;
        flops = execute(input_array, output_array, fmmplan, 0);
        double* ia=(double *)calloc(fmmplan->Nn, sizeof(double));
        flops = execute(output_array, ia, fmmplan, 1);
        gettimeofday(&t1, 0);
        double s1 = tdiff_sec(t0, t1);
        for (size_t j = 0; j < N; j++)
            error += pow(input_array[j]-ia[j], 2);
        printf("L2 Error back and forth %2.4e in %2.4e s\n", sqrt(error), s1);
        free(ia);
    }
    else
    {
        for (size_t i = 0; i < repeat; i++)
        {
            struct timeval r0, r1;
            gettimeofday(&r0, 0);
            for (size_t j = 0; j < N; j++)
                output_array[j] = 0.0;
            flops = execute(input_array, output_array, fmmplan, direction);
            gettimeofday(&r1, 0);
            double s1 = tdiff_sec(r0, r1);
            min_time = s1 < min_time ? s1 : min_time;
        }
        size_t eflops = (4.5*fmmplan->s + 4*fmmplan->M + (6+2.5)*fmmplan->M*fmmplan->M/fmmplan->s )*N;
        gettimeofday(&t1, 0);
        printf("Timing %2.6e s average %2.6e min of %lu times %lu N %lu S %lu F\n", tdiff_sec(t0, t1)/repeat, min_time, repeat, fmmplan->Nn, get_number_of_submatrices(fmmplan->N, fmmplan->s, fmmplan->L), flops);
        printf("Estimated flops %lu\n", eflops);
        printf("Norm %2.8e\n", norm(output_array, N));

        if (N < 10001)
        {
            struct timeval r0, r1;
            gettimeofday(&r0, 0);
            size_t fd;
            double* od=(double *)calloc(N, sizeof(double));
            for (size_t i = 0; i < repeat; i++)
            {
                fd = direct(input_array, od, fmmplan->dplan, direction);
            }
            gettimeofday(&r1, 0);
            double s1 = tdiff_sec(r0, r1);
            printf("Norm %2.8e time %2.4e F %lu \n", norm(od, N), s1/repeat, fd);

            double err = 0.0;
            double errinf = 0.0;
            double omax = 0.0;
            for (size_t i = 0; i < N; i++)
            {
                double e0 = output_array[i]-od[i];
                err += pow(e0, 2);
                errinf = fmax(errinf, e0);
                omax = fmax(omax, output_array[i]);
            }
            printf("Error     %2.6e\n", sqrt(err));
            printf("Error inf %2.8e\n", errinf/omax);
            free(od);
        }
    }
    free(input_array);
    free(output_array);
    free_fmm(fmmplan);
}