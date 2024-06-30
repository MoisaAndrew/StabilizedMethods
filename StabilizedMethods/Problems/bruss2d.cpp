#include "../problems.h"

#include <math.h>
#include <stdlib.h>

static const double alpha = 0.1;

/* right-hand-side function */
static void fbruss2d(const unsigned* n, const double* x, const double* y, double* fy)
{
    unsigned ns_sqr = (*n) / 2;
    double ans = sqrt(ns_sqr);
    unsigned ns = (unsigned)(ans + 0.1);
    unsigned tmp = ns_sqr - ns;

    double uleft, uright, uup, ulow;
    double vleft, vright, vup, vlow;
    double bet;

    if (*x >= 1.1)
    {
        bet = 5.0;
    }
    else
    {
        bet = 0.0;
    }

    for (unsigned i = 1; i <= ns_sqr; i++)
    {
        if (i % ns == 1)
        {
            uleft = y[i + ns - 2];
            vleft = y[ns_sqr + i + ns - 2];
        }
        else
        {
            uleft = y[i - 2];
            vleft = y[ns_sqr + i - 2];
        }

        if (i % ns == 0)
        {
            uright = y[i - ns];
            vright = y[ns_sqr + i - ns];
        }
        else
        {
            uright = y[i];
            vright = y[ns_sqr + i];
        }

        if (i <= ns)
        {
            ulow = y[i + tmp - 1];
            vlow = y[ns_sqr + i + tmp - 1];
        }
        else
        {
            ulow = y[i - ns - 1];
            vlow = y[ns_sqr + i - ns - 1];
        }

        if (i > tmp)
        {
            uup = y[i - tmp - 1];
            vup = y[ns_sqr + i - tmp - 1];
        }
        else
        {
            uup = y[i + ns - 1];
            vup = y[ns_sqr + i + ns - 1];
        }

        double uij = y[i - 1];
        double vij = y[i + ns_sqr - 1];
        fy[i - 1] = 1.0 + uij * uij * vij - 4.4 * uij + alpha * ns_sqr * (uleft + uright + ulow + uup - 4.0 * uij);
        fy[i + ns_sqr - 1] = 3.4 * uij - uij * uij * vij + alpha * ns_sqr * (vleft + vright + vlow + vup - 4.0 * vij);

        double iy = (i - 1) / ns + 1.;
        double ix = (i - (iy - 1.) * ans);
        double yy = iy / ans;
        double xx = ix / ans;

        if (((xx - 0.3) * (xx - 0.3) + (yy - 0.6) * (yy - 0.6)) <= 0.01) {
            fy[i - 1] += bet;
        }
    }
}


void get_bruss2d(
    unsigned* n,
    FcnEqDiff* fcn,
    Rho* rho,
    double* x,
    double* h,
    double* xend,
    double** y)
{
    const unsigned N = 128;
    *n = 2 * N * N;

    *fcn = fbruss2d;

    *x = 0.0; *xend = 11.5;

    *h = 1e-7;

    const double ans = (double)N;
    const unsigned N_sqr = N * N;

    *y = (double*)malloc(*n * sizeof(double));
    for (unsigned j = 0; j < N; j++) {
        const double yy = j / ans;
        for (unsigned i = 0; i < N; i++) {
            (*y)[j * N + i] = 22.0 * yy * pow(1.0 - yy, 1.5);
        }
    }
    for (unsigned i = 0; i < N; i++) {
        const double xx = i / ans;
        for (unsigned j = 0; j < N; j++) {
            (*y)[N_sqr + j * N + i] = 27.0 * xx * pow(1.0 - xx, 1.5);
        }
    }
}