#include "../problems.h"


static void frober(const unsigned* n, const double* x, const double* y, double* fy)
{
    fy[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
    fy[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1];
    fy[2] = 3.0e7 * y[1] * y[1];
}


void get_rober(
    unsigned* n,
    FcnEqDiff* fcn,
    Rho* rho,
    double* x,
    double* h,
    double* xend,
    double** y)
{
    *n = 3;

    *fcn = frober;

    *x = 0.0; *xend = 1.0e6;

    *h = 1e-4;

    *y = (double*)malloc(*n * sizeof(double));
    (*y)[0] = 1; (*y)[1] = 0; (*y)[2] = 0;
}