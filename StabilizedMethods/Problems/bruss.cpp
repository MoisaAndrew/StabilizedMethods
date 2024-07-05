#include "bruss.h"


/* right-hand-side function */
static void fbruss(const unsigned* n, const double* x, const double* y, double* fy)
{
    const unsigned N = *n / 2;
    const double dx_sqr = 1.0 / ((N + 1) * (N + 1));

    fy[0] = 1.0 + y[0] * y[0] * y[1] - 4.0 * y[0] + (alpha / dx_sqr) * (1.0 - 2 * y[0] + y[2]);
    fy[1] = 3.0 * y[0] - y[0] * y[0] * y[1] + (alpha / dx_sqr) * (3.0 - 2 * y[1] + y[3]);
    for (unsigned i = 2; i < N; i++)
    {
        fy[2 * i - 2] = 1.0 + y[2 * i - 2] * y[2 * i - 2] * y[2 * i - 1] - 4.0 * y[2 * i - 2] +
            (alpha / dx_sqr) * (y[2 * i - 4] - 2 * y[2 * i - 2] + y[2 * i]);
        fy[2 * i - 1] = 3.0 * y[2 * i - 2] - y[2 * i - 2] * y[2 * i - 2] * y[2 * i - 1] +
            (alpha / dx_sqr) * (y[2 * i - 3] - 2 * y[2 * i - 1] + y[2 * i + 1]);
    }
    fy[2 * N - 2] = 1.0 + y[2 * N - 2] * y[2 * N - 2] * y[2 * N - 1] - 4.0 * y[2 * N - 2] +
        (alpha / dx_sqr) * (y[2 * N - 4] - 2 * y[2 * N - 2] + 1.0);
    fy[2 * N - 1] = 3.0 * y[2 * N - 2] - y[2 * N - 2] * y[2 * N - 2] * y[2 * N - 1] +
        (alpha / dx_sqr) * (y[2 * N - 3] - 2 * y[2 * N - 1] + 3.0);
}


void get_bruss(ProblemParams** params, FcnEqDiff* fcn, Rho* rho)
{
    *params = new BrussParams();

    *fcn = fbruss;
}