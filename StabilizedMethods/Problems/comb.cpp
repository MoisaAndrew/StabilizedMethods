#include "../problems.h"


static const double d = 1.0;
static const double R = 5.0;
static const double alpha = 1.0;
static const double delta = 20.0;


/* right-hand-side function */
static void fcomb(const unsigned* n, const double* x, const double* y, double* fy)
{
    const unsigned N = (unsigned)(0.5 + sqrt(*n));
    const double h2 = 1. / ((1 + N) * (1 + N));
    unsigned i, j, jN;

    fy[0] = d * (y[1] - 2 * y[0] + y[N] - 2 * y[0]) / h2 +
        R * (1 + alpha - y[0]) * exp(delta * (1 - 1 / y[0])) / (alpha * delta);
    for (i = 1; i < N - 1; i++)
    {
        fy[i] = d * (y[i + 1] - 2 * y[i] + y[i - 1] + y[i + N] - 2 * y[i]) / h2 +
            R * (1 + alpha - y[i]) * exp(delta * (1 - 1 / y[i])) / (alpha * delta);
    }
    fy[N - 1] = d * (1 - 2 * y[N - 1] + y[N - 2] + y[2 * N - 1] - 2 * y[N - 1]) / h2 +
        R * (1 + alpha - y[N - 1]) * exp(delta * (1 - 1 / y[N - 1])) / (alpha * delta);

    for (j = 1; j < N - 1; j++)
    {
        jN = j * N;
        fy[jN] = d * (y[jN + 1] - 2 * y[jN] + y[jN + N] - 2 * y[jN] + y[jN - N]) / h2 +
            R * (1 + alpha - y[jN]) * exp(delta * (1 - 1 / y[jN])) / (alpha * delta);
        for (i = 1; i < N - 1; i++)
        {
            fy[jN + i] = d * (y[jN + i + 1] - 2 * y[jN + i] + y[jN + i - 1] + y[jN + i + N] - 2 * y[jN + i] + y[jN + i - N]) / h2 +
                R * (1 + alpha - y[jN + i]) * exp(delta * (1 - 1 / y[jN + i])) / (alpha * delta);
        }
        fy[jN + N - 1] = d * (1 - 2 * y[jN + N - 1] + y[jN + N - 2] + y[jN + 2 * N - 1] - 2 * y[jN + N - 1] + y[jN - 1]) / h2 +
            R * (1 + alpha - y[jN + N - 1]) * exp(delta * (1 - 1 / y[jN + N - 1])) / (alpha * delta);
    }

    jN = (N - 1) * N;
    fy[jN] = d * (y[jN + 1] - 2 * y[jN] + 1 - 2 * y[jN] + y[jN - N]) / h2 +
        R * (1 + alpha - y[jN]) * exp(delta * (1 - 1 / y[jN])) / (alpha * delta);
    for (i = 1; i < N - 1; i++)
    {
        fy[jN + i] = d * (y[jN + i + 1] - 2 * y[jN + i] + y[jN + i - 1] + 1 - 2 * y[jN + i] + y[jN + i - N]) / h2 +
            R * (1 + alpha - y[jN + i]) * exp(delta * (1 - 1 / y[jN + i])) / (alpha * delta);
    }
    fy[jN + N - 1] = d * (1 - 2 * y[jN + N - 1] + y[jN + N - 2] + 1 - 2 * y[jN + N - 1] + y[jN - 1]) / h2 +
        R * (1 + alpha - y[jN + N - 1]) * exp(delta * (1 - 1 / y[jN + N - 1])) / (alpha * delta);
}


void get_comb(
    unsigned* n,
    FcnEqDiff* fcn,
    Rho* rho,
    double* x,
    double* h,
    double* xend,
    double** y)
{
    const unsigned N = 80;
    *n = N * N;

    *fcn = fcomb;

    *x = 0.0; *xend = 0.32;

    *h = 1e-7;

    *y = (double*)malloc(*n * sizeof(double));
    for (unsigned i = 0; i < *n; i++)
    {
        (*y)[i] = 1.0;
    }
}