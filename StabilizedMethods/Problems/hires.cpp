#include "../problems.h"


static void fhires(const unsigned* n, const double* x, const double* y, double* fy)
{
    fy[0] = -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007;
    fy[1] = 1.71 * y[0] - 8.75 * y[1];
    fy[2] = -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4];
    fy[3] = 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3];
    fy[4] = -1.745 * y[4] + 0.43 * y[5] + 0.43 * y[6];
    fy[5] = -280. * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6];
    fy[6] = 280. * y[5] * y[7] - 1.81 * y[6];
    fy[7] = -280. * y[5] * y[7] + 1.81 * y[6];
}


void get_hires(
    unsigned* n,
    FcnEqDiff* fcn,
    double* x,
    double* h,
    double* xend,
    double** y)
{
    *n = 8;

    *fcn = fhires;

    *x = 0.0; *xend = 321.8122;

    *h = 1e-5;
    
    *y = (double*)malloc(*n * sizeof(double));
    (*y)[0] = 1.0; (*y)[1] = 0.0; (*y)[2] = 0.0; (*y)[3] = 0.0;
    (*y)[4] = 0.0; (*y)[5] = 0.0; (*y)[6] = 0.0; (*y)[7] = 0.0057;
}