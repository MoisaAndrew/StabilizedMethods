#include "raddiff.h"


/* right-hand-side function */
static void fraddiff(const unsigned* n, const double* x, const double* y, double* fy)
{
    const unsigned N = (unsigned)(0.5 + 0.5 * *n);
    const double h = 1. / N;

    double sigmai, sigmaim1, sigmaip1;
    double D1i, D1im1, D1ip1;
    double Evirt;

    unsigned i = 0;
    sigmai = 1 / pow(y[i + 1], 3);
    sigmaim1 = sigmai;
    sigmaip1 = 1 / pow(y[i + 3], 3);
    Evirt = (1 - (0.125 - 1 / (6 * sigmaim1 * h)) * y[i]) / (0.125 + 1 / (6 * sigmaim1 * h));
    D1i = 1 / (3 * sigmai + fabs(y[i + 2] - Evirt) / (2 * h * y[i]));
    D1im1 = 1 / (3 * sigmaim1 + fabs(y[i] - Evirt) / (h * Evirt));
    D1ip1 = 1 / (3 * sigmaip1 + fabs(y[i + 4] - y[i]) / (2 * h * y[i + 2]));
    fy[i] = ((D1i + D1ip1) * (y[i + 2] - y[i]) - (D1i + D1im1) * (y[i] - Evirt)) / (2 * h * h) 
        + sigmai * (pow(y[i + 1], 4) - y[i]);
    fy[i + 1] = k * (pow(y[i + 1], 2.5) + pow(y[i + 3], 2.5)) * (y[i + 3] - y[i + 1]) / (2 * h * h) 
        - sigmai * (pow(y[i + 1], 4) - y[i]);

    i += 2;
    sigmai = sigmaip1;
    sigmaip1 = 1 / pow(y[i + 3], 3);
    D1i = 1 / (3 * sigmai + fabs(y[i + 2] - y[i - 2]) / (2 * h * y[i]));
    D1im1 = 1 / (3 * sigmaim1 + fabs(y[i] - Evirt) / (2 * h * y[i - 2]));
    D1ip1 = 1 / (3 * sigmaip1 + fabs(y[i + 4] - y[i]) / (2 * h * y[i + 2]));
    fy[i] = ((D1i + D1ip1) * (y[i + 2] - y[i]) - (D1i + D1im1) * (y[i] - y[i - 2])) / (2 * h * h)
        + sigmai * (pow(y[i + 1], 4) - y[i]);
    fy[i + 1] = k * (
        (pow(y[i + 1], 2.5) + pow(y[i + 3], 2.5)) * (y[i + 3] - y[i + 1]) 
        - (pow(y[i + 1], 2.5) + pow(y[i - 1], 2.5)) * (y[i + 1] - y[i - 1])
        ) / (2 * h * h) - sigmai * (pow(y[i + 1], 4) - y[i]);

    for (i = 4; i < 2 * N - 4; i += 2)
    {
        sigmaim1 = sigmai;
        sigmai = sigmaip1;
        sigmaip1 = 1 / pow(y[i + 3], 3);
        D1i = 1 / (3 * sigmai + fabs(y[i + 2] - y[i - 2]) / (2 * h * y[i]));
        D1im1 = 1 / (3 * sigmaim1 + fabs(y[i] - y[i - 4]) / (2 * h * y[i - 2]));
        D1ip1 = 1 / (3 * sigmaip1 + fabs(y[i + 4] - y[i]) / (2 * h * y[i + 2]));
        fy[i] = ((D1i + D1ip1) * (y[i + 2] - y[i]) - (D1i + D1im1) * (y[i] - y[i - 2])) / (2 * h * h)
            + sigmai * (pow(y[i + 1], 4) - y[i]);
        fy[i + 1] = k * (
            (pow(y[i + 1], 2.5) + pow(y[i + 3], 2.5)) * (y[i + 3] - y[i + 1])
            - (pow(y[i + 1], 2.5) + pow(y[i - 1], 2.5)) * (y[i + 1] - y[i - 1])
            ) / (2 * h * h) - sigmai * (pow(y[i + 1], 4) - y[i]);
    }

    sigmaim1 = sigmai;
    sigmai = sigmaip1;
    sigmaip1 = 1 / pow(y[i + 3], 3);
    Evirt = -(0.125 - 1 / (6 * sigmaip1 * h)) * y[i + 2] / (0.125 + 1 / (6 * sigmaip1 * h));
    D1i = 1 / (3 * sigmai + fabs(y[i + 2] - y[i - 2]) / (2 * h * y[i]));
    D1im1 = 1 / (3 * sigmaim1 + fabs(y[i] - y[i - 4]) / (2 * h * y[i - 2]));
    D1ip1 = 1 / (3 * sigmaip1 + fabs(Evirt - y[i]) / (2 * h * y[i + 2]));
    fy[i] = ((D1i + D1ip1) * (y[i + 2] - y[i]) - (D1i + D1im1) * (y[i] - y[i - 2])) / (2 * h * h)
        + sigmai * (pow(y[i + 1], 4) - y[i]);
    fy[i + 1] = k * (
        (pow(y[i + 1], 2.5) + pow(y[i + 3], 2.5)) * (y[i + 3] - y[i + 1])
        - (pow(y[i + 1], 2.5) + pow(y[i - 1], 2.5)) * (y[i + 1] - y[i - 1])
        ) / (2 * h * h) - sigmai * (pow(y[i + 1], 4) - y[i]);

    i += 2;
    sigmaim1 = sigmai;
    sigmai = sigmaip1;
    D1i = 1 / (3 * sigmai + fabs(Evirt - y[i - 2]) / (2 * h * y[i]));
    D1im1 = 1 / (3 * sigmaim1 + fabs(y[i] - y[i - 4]) / (2 * h * y[i - 2]));
    D1ip1 = 1 / (3 * sigmaip1 + fabs(Evirt - y[i]) / (h * Evirt));
    fy[i] = ((D1i + D1ip1) * (Evirt - y[i]) - (D1i + D1im1) * (y[i] - y[i - 2])) / (2 * h * h)
        + sigmai * (pow(y[i + 1], 4) - y[i]);
    fy[i + 1] = -k * (pow(y[i + 1], 2.5) + pow(y[i - 1], 2.5)) * (y[i + 1] - y[i - 1])
        / (2 * h * h) - sigmai * (pow(y[i + 1], 4) - y[i]);
}


static void rhoraddiff(const unsigned* n, const double* x, const double* y, double* eigmax)
{
    const unsigned N = (unsigned)(0.5 + 0.5 * *n);
    const double h = 1. / N;

    *eigmax = 8 / (h * h) + 6000;
}


void get_raddiff(ProblemParams** params, FcnEqDiff* fcn, Rho* rho)
{
    *params = new RadDiffParams();

    *fcn = fraddiff;
    *rho = rhoraddiff;
}