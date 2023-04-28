#include "leg2cheb.h"

long lmin(long a, long b){
    return (a > b) ? b : a;
}

size_t min(size_t a, size_t b){
    return (a > b) ? b : a;
}

size_t max(size_t a, size_t b){
    return (a > b) ? a : b;
}

double fmax(double a, double b){
    return (a > b) ? a : b;
}

double norm(const double* a, const size_t N)
{
    double x = 0;
    for (size_t i = 0; i < N; i++)
    {
        x += a[i]*a[i];
    }
    return sqrt(x);
}

double norm_inf(const double* a, const size_t N)
{
    double x = 0;
    for (size_t i = 0; i < N; i++)
    {
        x = fmax(x, a[i]);
    }
    return x;
}

double chebval(const double x, const double* c, size_t M)
{
    double x2 = 2*x;
    double c0 = c[M-2];
    double c1 = c[M-1];
    for (size_t i = 3; i < M+1; i++)
    {
        double tmp = c0;
        c0 = c[M-i] - c1;
        c1 = tmp + c1*x2;
    }
    return c0 + c1*x;
}

double Lambda(const double z)
{
    double a0[8] = {9.9688251374224490e-01,
                    -3.1149502185860763e-03,
                    2.5548605043159494e-06,
                    1.8781800445057731e-08,
                    -4.0437919461099256e-11,
                    -9.0003384202201278e-13,
                    3.1032782098712292e-15,
                    1.0511830721865363e-16};

    if (z < 20)
        return exp(lgamma(z+0.5) - lgamma(z+1));
    return chebval(-1.0+40.0/z, a0, 7)/sqrt(z); // Note - using only 7
}

double mxy(const double *x, const double *y){
    return Lambda(((*y)-(*x))/2)*Lambda(((*y)+(*x))/2);
}

double lxy(const double *x, const double *y){
    return -1.0/(((*x)+(*y)+1)*((*y)-(*x)))*Lambda(((*y)-(*x)-2)/2)*Lambda(((*y)+(*x)-1)/2);
}

double tdiff_sec(struct timeval t0, struct timeval t1)
{
    return (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) / 1000000.0;
}

size_t get_number_of_blocks(const size_t level)
{
    return pow(2, level+1)-1;
}

size_t get_h(const size_t level, const size_t L)
{
    return pow(2, L-level-1);
}

size_t get_number_of_submatrices(const size_t N, const size_t s, const size_t L)
{
    return 3*N/(2*s) - 3*(L+2);
}

void get_ij(size_t* ij, const size_t level, const size_t block, const size_t s, const size_t L)
{
    ij[0] = 2*block*s*get_h(level, L);
    ij[1] = ij[0]+2*s*get_h(level, L);
}

size_t direct(const double* u, double* b, fmm_plan* plan, unsigned direction)
{
    const size_t N = plan->N;
    size_t flops = 0;
    for (size_t i = 0; i < N; i++)
        b[i] = 0.0;
    flops += N;

    if (direction == 0)
    {
        const double* a = plan->a;

        for (size_t n = 0; n < N/2+N%2; n++)
        {
            const double* ap = &a[n];
            const double* cp = &u[2*n];
            const double a0 = ap[0];
            for (size_t i = 0; i < N-2*n; i++)
                b[i] += a0*ap[i]*cp[i];
            flops += 3*(N-2*n);
        }
        b[0] /= 2;
        for (size_t i = 0; i < N; i++)
            b[i] *= M_2_PI;
        flops += N;
    }
    else
    {
        double* vn = (double*)malloc(N*sizeof(double));
        const double* an = plan->an;
        const double* dn = plan->dn;
        vn[0] = u[0];
        for (size_t i = 1; i < N; i++)
            vn[i] = u[i]*i;

        for (size_t n = 0; n < N; n++)
            b[n] = sqrt(M_PI)*an[n]*vn[n];

        for (size_t n = 2; n < N; n=n+2)
        {
            const double* ap = &an[n/2];
            const double* vp = &vn[n];
            for (size_t i = 0; i < N-n; i++)
                b[i] -= dn[n/2]*ap[i]*vp[i];
            flops += 3*(N-2*n);
        }
        for (size_t i = 0; i < N; i++)
            b[i] *= (i+0.5);
        flops += N;
        free(vn);
    }
    return flops;
}

size_t directM(double* input_array, double* output_array, fmm_plan* fmmplan)
{
    size_t s = fmmplan->s;
    size_t N = fmmplan->N;
    double* a = fmmplan->a;

    size_t flops = 0;
    for (size_t n = 0; n < s; n++)
    {
        const double* ap = &a[n];
        const double* cp = &input_array[2*n];
        const double a0 = ap[0];
        double* op = &output_array[0];
        for (size_t i = 0; i < N-2*n; i++)
            (*op++) += a0*(*ap++)*(*cp++);
        flops += (N-2*n)*3;
    }

    size_t h = 2*s;
    size_t n0 = 2*s;
    for (size_t block = 0; block < get_number_of_blocks(fmmplan->L); block++)
    {
        size_t i0 = block*h;
        size_t j0 = n0+i0;
        double* vp = &output_array[i0];
        const long Nm = lmin(N-j0, h);
        for (size_t n = 0; n < h; n=n+2)
        {
            const size_t n1 = (n+n0)/2;
            const double* ap = &a[i0+n1];
            const double a0 = a[n1];
            const double* up = &input_array[j0+n];
            //printf("Nm %ld %ld %ld   %ld %ld\n", Nm, n, Nm-n, N-j0, h);
            if ((long)(Nm-n) < 0)
                break;
            for (size_t i = 0; i < (size_t)(Nm-n); i++)
                vp[i] += a0*(*ap++)*(*up++);
            flops += 3*(Nm-n);
        }
    }
    double* op = &output_array[0];
    for (size_t i = 0; i < N; i++)
        (*op++) *= M_2_PI;
    output_array[0] *= 0.5;
    return flops;
}

size_t directL(double* input_array, double* output_array, fmm_plan* fmmplan)
{
    size_t s = fmmplan->s;
    size_t N = fmmplan->Nn;
    double* an = fmmplan->an;
    double* dn = fmmplan->dn;
    double SPI = sqrt(M_PI);

    size_t flops = 0;
    for (size_t i = 0; i < N; i++)
        output_array[i] += SPI*an[i]*input_array[i];

    for (size_t n = 1; n < s; n++)
    {
        const double* ap = &an[n];
        const double* cp = &input_array[2*n];
        double* op = &output_array[0];
        for (size_t i = 0; i < N-2*n; i++)
            (*op++) -= dn[n]*(*ap++)*(*cp++);
        flops += (N-2*n)*3;
    }

    size_t h = 2*s;
    size_t n0 = 2*s;
    for (size_t block = 0; block < get_number_of_blocks(fmmplan->L); block++)
    {
        size_t i0 = block*h;
        size_t j0 = n0+i0;
        double* vp = &output_array[i0];
        const long Nm = lmin(N-j0, h);
        for (size_t n = 0; n < h; n=n+2)
        {
            const size_t n1 = (n+n0)/2;
            const double* ap = &an[i0+n1];
            const double* up = &input_array[j0+n];
            if ((long)(Nm-n) < 0)
                break;
            for (size_t i = 0; i < (size_t)(Nm-n); i++)
                vp[i] -= dn[n1]*(*ap++)*(*up++);
            flops += 3*(Nm-n);
        }
    }
    double* op = &output_array[0];
    for (size_t i = 0; i < N; i++)
        (*op++) *= (i+0.5);
    return flops;
}

void matvectri(const double* A, const double* x, double* b, const size_t m, const size_t n, const size_t lda, const size_t upper)
{
    if (upper == 0)
    {
        for (size_t i = 0; i < m; i++)
        {
            double s = 0.0;
            const double* xp = &x[0];
            const double* ap = &A[i*lda];
            for (size_t j = 0; j < i+1; j++)
                s += (*ap++)*(*xp++);
            (*b++) += s;
        }
    }
    else
    {
        for (size_t i = 0; i < m; i++)
        {
            double s = 0.0;
            const double* xp = &x[i];
            const double* ap = &A[i*lda+i];
            for (size_t j = i; j < n; j++)
                s += (*ap++)*(*xp++);
            (*b++) += s;
        }
    }
}

void matvectriZ(const double* A, const double* x, double* b, const size_t m, const size_t n, const size_t lda, const size_t upper)
{

    if (upper == 0)
    {
        double* zp = (double*)malloc(m*sizeof(double));
        double* zm = (double*)malloc(m*sizeof(double));
        const double *xp = &x[lda];
        for (size_t i = 0; i < m; i=i+2)
        {
            zp[i] = x[i]+xp[i];
            zm[i] = x[i]-xp[i];
        }
        for (size_t i = 1; i < m; i=i+2)
        {
            zp[i] = x[i]-xp[i];
            zm[i] = x[i]+xp[i];
        }
        for (size_t i = 0; i < m; i++)
        {
            double s = 0.0;
            const double* z = i % 2 == 0 ? &zp[0] : &zm[0];
            const double* ap = &A[i*lda];
            for (size_t j = 0; j < i+1; j++)
                s += (*ap++)*(*z++);
            (*b++) += s;
        }
        free(zp);
        free(zm);
    }
    else
    {
        double *bp = &b[lda];
        for (size_t i = 0; i < m; i++)
        {
            const double* ap = &A[i*lda+i];
            const double* xp = &x[i];
            for (size_t j = i; j < n-1; j=j+2)
            {
                double se = (*ap++)*(*xp++);
                double so = (*ap++)*(*xp++);
                (*b) += se+so;
                (*bp) += se-so;
            }
            double s = (*ap++)*(*xp++);
            (*b) += s;
            (*bp) += s;
            b++;
            bp++;
        }
    }
}

void matvec(const double* A, const double* x, double* b, const size_t m, const size_t n, const size_t transpose)
{
    if (transpose == 0)
    {
        //cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, n, x, 1, 1.0, b, 1);
        for (size_t i = 0; i < m; i++)
        {
            double s = 0.0;
            const double* xp = &x[0];
            for (size_t j = 0; j < n; j++)
                s += (*A++)*(*xp++);
            (*b++) += s;
        }
    }
    else
    {
        //cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, A, n, x, 1, 0.0, b, 1);
        for (size_t i = 0; i < m; i++)
        {
            double s = (*x++);
            double* bp = &b[0];
            for (size_t j = 0; j < n; j++)
                (*bp++) += s*(*A++);
        }
    }
}

void vandermonde(double* T, const size_t h, const size_t M)
{
    double* x = (double*)malloc(2*h*sizeof(double));
    double* x2 = (double*)malloc(2*h*sizeof(double));
    double* Tm = (double*)malloc(2*h*M*sizeof(double));

    for (size_t i = 0; i < 2*h; i++){
        x[i] = -1+(double)(i)/((double)h);
        x2[i] = 2*x[i];
        Tm[i*M] = (double)(1);
        Tm[i*M+1] = x[i];
    }

    for (size_t m = 2; m < M; m++)
    {
        for (size_t i = 0; i < 2*h; i++)
        {
            Tm[i*M+m] = Tm[i*M+m-1]*x2[i] - Tm[i*M+m-2];
        }
    }

    for (size_t i = 0; i < h; i++)
    {
        for (size_t m = 0; m < M; m++)
        {
            T[i*M+m] = Tm[2*i*M+m];          // even
            T[(h+i)*M+m] = Tm[(2*i+1)*M+m];  // odd
        }
    }
    free(x);
    free(x2);
    free(Tm);
}

double sum(const double* a, const size_t N)
{
    double x = 0;
    for (size_t i = 0; i < N; i++)
    {
        x += a[i];
    }
    return x;
}

void free_fmm(fmm_plan* plan)
{
    free(plan->A[0]);
    free(plan->A[1]);
    free(plan->A);
    free(plan->a);
    free(plan->an);
    free(plan->dn);
    free(plan->T);
    free(plan->TT);
    free(plan->Th);
    free(plan->ThT);
}

fmm_plan* create_fmm(size_t N, size_t maxs, unsigned direction)
{
    size_t M = 18;
    fftw_plan plan, plan1d;
    static fmm_plan fmmplan;
    size_t Nn;
    size_t L = 1;
    size_t s;
    size_t ij[2];
    struct timeval t1, t2;
    unsigned* directions=(unsigned*)malloc(2*sizeof(unsigned));
    unsigned num_directions = 2;
    switch (direction)
    {
    case 0:
        directions[0] = 0;
        num_directions = 1;
        directions = realloc(directions, sizeof(unsigned));
        break;
    case 1:
        directions[0] = 1;
        num_directions = 1;
        directions = realloc(directions, sizeof(unsigned));
        break;
    default:
        directions[0] = 0;
        directions[1] = 1;
        break;
    }

    L = ceil(log2((double)N/(double)maxs))-2;
    if (L < 1)
    {
        printf("Levels < 1. Try smaller maxs. Exiting...\n");
        exit(-1);
    }

    s = ceil((double)N/(double)pow(2, L+2));
    Nn = s*pow(2, L+2);

    fmmplan.M = M;
    fmmplan.L = L;
    fmmplan.N = N;
    fmmplan.Nn = Nn;
    fmmplan.s = s;

    printf("N %lu\n", N);
    printf("Num levels %lu\n", L);
    printf("Num submatrices %lu\n", get_number_of_submatrices(Nn, s, L));
    printf("Given max s %lu \n", maxs);
    printf("Computed s %lu \n", s);
    printf("Computed N %lu\n", Nn);

    gettimeofday(&t1, 0);

    double** A = (double**)malloc(2*sizeof(double*));
    if (direction == 2)
    {
        A[0] = (double *)malloc(get_number_of_submatrices(Nn, s, L)*M*M*sizeof(double));
        A[1] = (double *)malloc(get_number_of_submatrices(Nn, s, L)*M*M*sizeof(double));
    }
    else
    {
        A[direction] = (double *)malloc(get_number_of_submatrices(Nn, s, L)*M*M*sizeof(double));
    }

    double* a = (double*)malloc(N*sizeof(double));
    double* dn = (double*)malloc(N/2*sizeof(double));
    double* an = (double*)malloc(N*sizeof(double));
    double* xj = (double*)malloc(M*sizeof(double));
    for (size_t i = 0; i < N; i++)
        a[i] =  Lambda(i);
    dn[0] = 0;
    an[0] = 2/sqrt(M_PI);
    for (size_t i = 1; i < N; i++)
        an[i] = 1/(2*Lambda(i)*i*(i+0.5));
    for (size_t i = 1; i < N/2; i++)
        dn[i] = Lambda((double)(2*i-2)/2.0)/(2*i);
    for (size_t i = 0; i < M; i++)
        xj[i] = cos((i+0.5)*M_PI/M);

    double* fun = (double*)fftw_malloc(M*M*sizeof(double));
    double* fun_hat = (double*)fftw_malloc(M*M*sizeof(double));
    plan1d = fftw_plan_r2r_1d(M, fun, fun_hat, FFTW_REDFT10, FFTW_ESTIMATE);
    plan = fftw_plan_r2r_2d(M, M, fun, fun_hat, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
    double* xj0 = (double *)malloc(M*sizeof(double));
    for (size_t i = 0; i < M; i++)
        xj0[i] = cos((i+0.5)*M_PI/M);

    for (size_t di = 0; di < num_directions; di++)
    {
        size_t dir = directions[di];
        double* ap = A[dir];
        for (size_t level = L; level --> 0 ;) // Reverse loop L-1, L-2, ..., 0
        {
            size_t h = s*get_h(level, L);
            for (size_t block = 0; block < get_number_of_blocks(level); block++)
            {
                get_ij(ij, level, block, s, L);
                for (size_t q = 0; q < 2; q++)
                {
                    for (size_t p = 0; p < q+1; p++)
                    {
                        for (size_t i = 0; i < M; i++)
                        {
                            double x = (2*(ij[0]+p*h)+(xj0[i]+1)*h);
                            for (size_t j = 0; j < M; j++)
                            {
                                double y = (2*(ij[1]+(q+0)*h)+(xj0[j]+1)*h);
                                fun[i*M+j] = dir == 0 ? mxy(&x, &y) : lxy(&x, &y);
                            }
                        }
                        fftw_execute(plan);
                        for (size_t i = 0; i < M; i++){
                            for (size_t j = 0; j < M; j++){
                                if (i == 0 && j == 0)
                                    fun_hat[i*M+j] /= (4*M*M);
                                else if (i == 0)
                                    fun_hat[i*M+j] /= (2*M*M);
                                else if (j == 0)
                                    fun_hat[i*M+j] /= (2*M*M);
                                else
                                    fun_hat[i*M+j] /= (M*M);
                            }
                        }
                        double* fh = &fun_hat[0];
                        for (size_t i = 0; i < M*M; i++)
                            *ap++ = *fh++;
                    }
                }
            }
        }
    }

    fmmplan.A = A;
    fmmplan.a = a;
    fmmplan.dn = dn;
    fmmplan.an = an;

    double* T = (double*)malloc(2*s*M*sizeof(double));
    double* TT = (double*)malloc(2*s*M*sizeof(double));
    vandermonde(T, s, M);
    for (size_t i = 0; i < s; i++)
    {
        for (size_t m = 0; m < M; m++)
        {
            TT[m*s+i] = T[i*M+m];          // even
            TT[s*(M+m)+i] = T[(s+i)*M+m];  // odd
        }
    }
    fmmplan.T = T;
    fmmplan.TT = TT;

    double* Th = (double *)malloc(2*M*M*sizeof(double));
    double* ThT = (double *)malloc(2*M*M*sizeof(double));
    double* th = &Th[0];
    for (size_t q = 0; q < 2; q++)
    {
        for (size_t k = 0; k < M; k++)
        {
            for (size_t j = 0; j < M; j++)
            {
                fun[j] = cos(k*acos((xj0[j]+2*q-1)/2));
                fun_hat[j] = 0.0;
            }
            fftw_execute(plan1d);
            *th++ = fun_hat[0]/M/2;
            for (size_t j = 1; j < M; j++)
                *th++ = fun_hat[j]/M;
        }
    }
    for (size_t q = 0; q < 2; q++)
    {
        th = &Th[q*M*M];
        double* tht = &ThT[q*M*M];
        for (size_t i = 0; i < M; i++)
        {
            for (size_t j = 0; j < M; j++)
            {
                tht[i*M+j] = th[j*M+i];
            }
        }
    }

    fftw_destroy_plan(plan);
    fftw_destroy_plan(plan1d);
    free(xj);
    free(xj0);
    fmmplan.Th = Th;
    fmmplan.ThT = ThT;
    gettimeofday(&t2, 0);
    printf("Initialization %2.4e s\n", tdiff_sec(t1, t2));
    return &fmmplan;
}

size_t execute(const double* input_array, double* output_array, fmm_plan* fmmplan, unsigned direction)
{
    size_t Nn = fmmplan->Nn;
    size_t N = fmmplan->N;
    size_t L = fmmplan->L;
    size_t s = fmmplan->s;
    size_t M = fmmplan->M;
    double* TT = fmmplan->TT;
    double* T = fmmplan->T;
    double* Th = fmmplan->Th;
    double* ThT = fmmplan->ThT;
    double* A = fmmplan->A[direction];
    double* ia = (double *)malloc(Nn/2*sizeof(double));
    double* oa = (double *)calloc(Nn/2, sizeof(double));
    double** wk = (double **)malloc(L*sizeof(double*));
    double** ck = (double **)malloc(L*sizeof(double*));
    double* input = (double*)malloc(Nn*sizeof(double));
    memcpy(input, input_array, N*sizeof(double));
    for (size_t i = N; i < Nn; i++)
        input[i] = 0.0;

    if (direction == 1)
    {
        double* w0 = &input[2];
        for (size_t i = 2; i < N; i++)
            (*w0++) = input_array[i]*i;
    }

    for (size_t level = 0; level < L; level++)
    {
        size_t b = get_number_of_blocks(level);
        wk[level] = (double*)calloc(b*2*M, sizeof(double));
        ck[level] = (double*)calloc(b*2*M, sizeof(double));
    }
    size_t flops = 0;
    for (size_t odd = 0; odd < 2; odd++)
    {
        if (odd == 1)
        {
            for (size_t i = 0; i < Nn/2; i++){
                ia[i] = 0;
                oa[i] = 0;
            }
            for (size_t level = 0; level < L; level++)
            {
                double* wl = wk[level];
                double* cl = ck[level];
                for (size_t i = 0; i < get_number_of_blocks(level)*2*M; i++)
                {
                    *wl++ = 0.0;
                    *cl++ = 0.0;
                }
            }
        }
        const double* ap = &input[odd];
        unsigned int rest = N % 2;
        double* iap = &ia[0];
        for (size_t i = 0; i < Nn/2; i++)
        {
            *iap++ = *ap;
            ap += 2;
        }
        size_t Nc = 0;
        size_t ik = 0;
        for (size_t level = L; level --> 0;)
        {
            if (level == L-1){
                size_t K = get_number_of_blocks(L-1)*2;
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, K, M, s, 1.0, &ia[2*s], s, &T[odd*s*M], M, 0.0, wk[level], M);
                flops += 2*K*s*M;
            }

            for (size_t block = 0; block < get_number_of_blocks(level); block++)
            {
                size_t Nd = block*2*M;
                double* c0 = &ck[level][Nd];
                double* wq = &wk[level][Nd];
#ifdef TRI
                if (level > 0 && block > 0) // optimization for constant block size 2
                {
                    int b0 = (block-1)/2;
                    int q0 = (block-1)%2;
                    matvectriZ(&Th[0], &wk[level][Nd], &wk[level-1][(b0*2+q0)*M], M, M, M, 0);
                    flops += M*M; //+2*M;
                }
#endif

                for (size_t q = 0; q < 2; q++)
                {
                    //if (level == L-1)
                    //{
                    //    size_t ij[2];
                    //    get_ij(ij, level, block, s, L);
                    //    matvec(&TT[odd*s*M0], &ia[ij[1]+q*s], &wq[q*M], M, s, 0);
                    //    //cblas_dgemv(CblasRowMajor, CblasNoTrans, M, s, 1.0, &TT[odd*s*M], s, &ia[ij[1]+q*s], 1, 0.0, &wq[q*M], 1);
                    //    flops += 2*s*M;
                    //}
#ifndef TRI
                    if (level > 0 && block > 0)
                    {
                        int b0 = (block-1)/2;
                        int q0 = (block-1)%2;
                        matvectri(&Th[q*M*M], &wq[q*M], &wk[level-1][(b0*2+q0)*M], M, M, M, 0);
                        flops += M*M;
                    }
#endif
                    for (size_t p = 0; p < q+1; p++)
                    {
                        //matvec(&A[Nc], &wq[q*M], &c0[p*M], M, M, 0);
                        cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1.0, &A[Nc], M, &wq[q*M], 1, 1.0, &c0[p*M], 1);
                        flops += 2*M*M;
                        Nc += M*M;
                        ik++;
                    }
                }
            }
        }
        for (size_t level = 0; level < L-1; level++)
        {
            double* c0 = ck[level];
            double* c1 = ck[level+1];
            for (size_t block = 0; block < get_number_of_blocks(level+1)-1; block++)
            {
#ifdef TRI
                matvectriZ(&ThT[0], &c0[block*M], &c1[block*2*M], M, M, M, 1);
                flops += 1.5*M*M;
#else
                for (size_t p = 0; p < 2; p++)
                {
                    matvectri(&ThT[p*M*M], &c0[block*M], &c1[j1*M], M, M, M, 1);
                    j1++;
                    flops += M*M;
                }
#endif
            }
        }
        //for (size_t block = 0; block < get_number_of_blocks(L-1, D, L); block++)
        //{
        //    size_t ij[2];
        //    size_t level = L-1;
        //    get_ij(ij, level, block, diagonals, s, D, L);
        //    double* c0 = &ck[level][block*D[level]*Mmax];
        //    const double* Tp = &T[odd*s*M0];
        //    for (size_t p = 0; p < D[level]; p++)
        //    {
        //        //cblas_dgemv(CblasRowMajor, CblasNoTrans, s, M0, 1.0, Tp, M0, &c0[p*Mmax], 1, 0.0, &oa[ij[0]+p*s], 1);
        //        matvec(Tp, &c0[p*Mmax], &oa[ij[0]+p*s], s, M0, 0);
        //        flops += 2*s*M0;
        //    }
        //}
        size_t K = get_number_of_blocks(L-1)*2;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, K, s, M, 1.0, ck[L-1], M, &TT[odd*s*M], s, 0.0, &oa[0], s);
        flops += 2*K*s*M;

        double* oaa = &output_array[odd];
        double* oap = &oa[0];
        for (size_t i = 0; i < N/2 + rest*(1-odd); i++)
        {
            *oaa += (*oap++);
            oaa += 2;
        }
        //flops += 3*N/2;
    }

    switch (direction)
    {
    case 0:
        flops += directM(input, output_array, fmmplan);
        break;

    default:
        flops += directL(input, output_array, fmmplan);
        break;
    }

    free(oa);
    free(ia);
    for (size_t level = 0; level < L; level++)
    {
        free(wk[level]);
        free(ck[level]);
    }
    free(wk);
    free(ck);
    free(input);
    return flops;
}
