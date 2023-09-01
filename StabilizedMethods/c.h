#include "methods_common.h"

extern "C" 
{
    int rock2c(const unsigned n, const double x, const double xend, double* h, double* y,
        const FcnEqDiff f, const Rho rho, const SolTrait solout,
        const double* atol, const double* rtol,
        unsigned iwork[12]);


    int rkcc(const unsigned n, const double x, const double xend, double* y,
        const FcnEqDiff f, const Rho rho, const SolTrait solout,
        const double* atol, const double rtol,
        unsigned iwork[10]);


    int tsrkc2(const unsigned n,
        const double x0, const double x1, const double xend,
        double* h, double* y0, double* y1,
        const FcnEqDiff f, const Rho rho, const SolTrait solout,
        const double* atol, const double* rtol,
        unsigned iwork[12]);
}