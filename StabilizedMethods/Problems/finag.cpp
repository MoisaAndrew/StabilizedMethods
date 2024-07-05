#include "finag.h"


static void ffinag(const unsigned* n, const double* x, const double* y, double* fy)
{
    const unsigned N = *n / 2;
    unsigned i;
    double fv;
    const double dis = pow(0.01 * N, 2), rhodn = 30. / N;

    fv = y[0] * (y[0] - alpha) * (y[0] - 1);
    fy[0] = dis * (rhodn - y[0] + y[2]) - fv - y[1];
    for (i = 2; i < *n - 3; i += 2)
    {
        fv = y[i] * (y[i] - alpha) * (y[i] - 1);
        fy[i] = dis * (y[i - 2] - 2 * y[i] + y[i + 2]) - fv - y[i + 1];
    }
    fv = y[*n - 2] * (y[*n - 2] - alpha) * (y[*n - 2] - 1);
    fy[*n - 2] = dis * (y[*n - 4] - y[*n - 2]) - fv - y[*n - 1];

    for (i = 1; i < *n; i += 2)
    {
        fy[i] = eta * (y[i - 1] + beta * y[i]);
    }
}


void get_finag(ProblemParams** params, FcnEqDiff* fcn, Rho* rho)
{
    *params = new FinagParams();

    *fcn = ffinag;
}