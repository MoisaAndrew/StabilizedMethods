#include "rober.h"


static void frober(const unsigned* n, const double* x, const double* y, double* fy)
{
    fy[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
    fy[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1];
    fy[2] = 3.0e7 * y[1] * y[1];
}


void get_rober(ProblemParams** params, FcnEqDiff* fcn, Rho* rho)
{
    *params = new RoberParams();

    *fcn = frober;
}