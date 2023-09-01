#include "../problems.h"


static const double mu = 0.0003;


static void fburgers(const unsigned* n, const double* x, const double* y, double* fy)
{
    double h = 1 / (double)(*n + 1);
    double h_sqr = h * h;
    double tmp = 0.25 / h;
    fy[0] = (mu / h_sqr) * (y[1] - 2. * y[0]) - tmp * (y[1] * y[1]);
    for (unsigned i = 1; i < *n - 1; i++)
    {
        fy[i] = (mu / h_sqr) * (y[i + 1] - 2 * y[i] + y[i - 1]) - tmp * (y[i + 1] * y[i + 1] - y[i - 1] * y[i - 1]);
    }
    fy[*n - 1] = (mu / h_sqr) * (-2. * y[*n - 1] + y[*n - 2]) + tmp * (y[*n - 2] * y[*n - 2]);
}


void get_burgers(
    unsigned* n,
    FcnEqDiff* fcn,
    double* x, 
    double* h,
    double* xend,
    double** y)
{
    *n = 500;

    *fcn = fburgers;

    *x = 0.0; *xend = 2.5;

    *h = 1e-4;

    *y = (double*)malloc(*n * sizeof(double));
    for (unsigned i = 0; i < *n; i++)
    {
        double temp = (double)(i + 1) / (double)(*n + 1);
        (*y)[i] = 1.5 * temp * (1. - temp) * (1. - temp);
    }
}