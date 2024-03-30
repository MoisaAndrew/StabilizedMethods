#include "../problems.h"


static void fheat(const unsigned* n, const double* x, const double* y, double* fy)
{
    const double h = 1. / (*n + 1);

    double xx = h;
    fy[0] = (y[1] - 2 * y[0]) / (h * h) - exp(-(*x)) * (2 + xx + xx * xx);
    for (unsigned i = 1; i < *n - 1; i++)
    {
        xx = (double)(i + 1) / (*n + 1);
        fy[i] = (y[i + 1] - 2 * y[i] + y[i - 1]) / (h * h) - exp(-(*x)) * (2 + xx + xx * xx);
    }
    xx = 1 - h;
    fy[*n - 1] = (2 * exp(-(*x)) - 2 * y[*n - 1] + y[*n - 2]) / (h * h) - exp(-(*x)) * (2 + xx + xx * xx);
}


/// <summary>
/// This problem can be used to test a method's convergence.
/// Just remove an adaptive step size selection from the method's code.
/// </summary>
void get_heat(
    unsigned* n,
    FcnEqDiff* fcn,
    double* x,
    double* h,
    double* xend,
    double** y)
{
    *n = 159;

    *fcn = fheat;

    *x = 0.0; *xend = 1.0;

    *h = 1 / 20.;//1e-5;

    *y = (double*)malloc(*n * sizeof(double));
    for (unsigned i = 0; i < *n; i++)
    {
        double x = (double)(i + 1) / (*n + 1);
        (*y)[i] = x * (x + 1);
    }
}