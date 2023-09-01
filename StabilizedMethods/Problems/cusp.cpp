#include "../problems.h"


static const double sigma = 1. / 144.;


static void fcusp(const unsigned* n, const double* x, const double* y, double* fy)
{
    const unsigned N = *n / 3;
    const double D = N * N * sigma;
   
    double v = (y[0] - 0.7) * (y[0] - 1.3) / (0.1 + (y[0] - 0.7) * (y[0] - 1.3));
    fy[0] = -1.e4 * (pow(y[0], 3) + y[1] * y[0] + y[2]) + D * (y[3 * N - 3] - 2 * y[0] + y[3]);
    fy[1] = y[2] + 0.07 * v + D * (y[3 * N - 2] - 2 * y[1] + y[4]);
    fy[2] = (1 - y[1] * y[1]) * y[2] - y[1] - 0.4 * y[0] + 0.035 * v + D * (y[3 * N - 1] - 2 * y[2] + y[5]);
    for (unsigned i = 2; i < N; i++)
    {
        v = (y[3 * i - 3] - 0.7) * (y[3 * i - 3] - 1.3) / (0.1 + (y[3 * i - 3] - 0.7) * (y[3 * i - 3] - 1.3));
        fy[3 * i - 3] = -1.e4 * (pow(y[3 * i - 3], 3) + y[3 * i - 2] * y[3 * i - 3] + 
            y[3 * i - 1]) + D * (y[3 * N - 6] - 2 * y[3 * i - 3] + y[3 * i]);
        fy[3 * i - 2] = y[3 * i - 1] + 0.07 * v + D * (y[3 * N - 5] - 2 * y[3 * i - 2] + y[3 * i + 1]);
        fy[3 * i - 1] = (1 - y[3 * i - 2] * y[3 * i - 2]) * y[3 * i - 1] - y[3 * i - 2] -
            0.4 * y[3 * i - 3] + 0.035 * v + D * (y[3 * N - 4] - 2 * y[3 * i - 1] + y[3 * i + 2]);
    }
    v = (y[3 * N - 3] - 0.7) * (y[3 * N - 3] - 1.3) / (0.1 + (y[3 * N - 3] - 0.7) * (y[3 * N - 3] - 1.3));
    fy[3 * N - 3] = -1.e4 * (pow(y[3 * N - 3], 3) + y[3 * N - 2] * y[3 * N - 3] + 
        y[3 * N - 1]) + D * (y[3 * N - 6] - 2 * y[3 * N - 3] + y[0]);
    fy[3 * N - 2] = y[3 * N - 1] + 0.07 * v + D * (y[3 * N - 5] - 2 * y[3 * N - 2] + y[1]);
    fy[3 * N - 1] = (1 - y[3 * N - 2] * y[3 * N - 2]) * y[3 * N - 1] - y[3 * N - 2] - 
        0.4 * y[3 * N - 3] + 0.035 * v + D * (y[3 * N - 4] - 2 * y[3 * N - 1] + y[2]);
}


void get_cusp(
    unsigned* n,
    FcnEqDiff* fcn,
    double* x,
    double* h,
    double* xend,
    double** y)
{
    const unsigned N = 32;
    *n = 3 * N;

    *fcn = fcusp;

    *x = 0.0; *xend = 1.1;

    *h = 1e-8;

    *y = (double*)malloc(*n * sizeof(double));
    for (unsigned j = 0; j < N; j++) {
        (*y)[3 * j] = 0;
        (*y)[3 * j + 1] = -2 * cos(2 * j * M_PI / N);
        (*y)[3 * j + 2] = 2 * sin(2 * j * M_PI / N);
    }
}