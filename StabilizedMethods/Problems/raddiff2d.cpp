#include "raddiff2d.h"


//                ijp2
//        im1jp1  ijp1  ip1jp1
//  im2j   im1j    ij    ip1j   ip2j
//        im1jm1  ijm1  ip1jm1
//                ijm2
static void fij(
    const double h,
    const double Eij, const double Tij, 
    double Eim1j, double Tim1j, 
    double Eip1j, double Tip1j,
    double Eijm1, double Tijm1, 
    double Eijp1, double Tijp1,
    double Eim2j, double Eip2j, double Eijm2, double Eijp2,
    double Eim1jm1, double Eim1jp1, double Eip1jm1, double Eip1jp1,
    double* fE, double* fT)
{
    double sigmaij = 1 / pow(Tij, 3);
    double sigmaim1j, sigmaip1j, sigmaijm1, sigmaijp1;

    if (isnan(Eim1j))
    {
        sigmaim1j = sigmaij;
        sigmaip1j = 1 / pow(Tip1j, 3);
        Eim1j = (1 - (0.125 - 1 / (6 * sigmaij * h)) * Eij) / (0.125 + 1 / (6 * sigmaij * h));
        Tim1j = Tij;
        if (!isnan(Eijm1))
        {
            Eim1jm1 = (1 - (0.125 - pow(Tijm1, 3) / (6 * h)) * Eijm1) / (0.125 + pow(Tijm1, 3) / (6 * h));
        }
        if (!isnan(Eijp1))
        {
            Eim1jp1 = (1 - (0.125 - pow(Tijp1, 3) / (6 * h)) * Eijp1) / (0.125 + pow(Tijp1, 3) / (6 * h));
        }
    }
    else
    {
        sigmaim1j = 1 / pow(Tim1j, 3);
        if (isnan(Eip1j))
        {
            sigmaip1j = sigmaij;
            Eip1j = -(0.125 - 1 / (6 * sigmaij * h)) * Eij / (0.125 + 1 / (6 * sigmaij * h));
            Tip1j = Tij;
            if (!isnan(Eijm1))
            {
                Eip1jm1 =  -(0.125 - pow(Tijm1, 3) / (6 * h)) * Eijm1 / (0.125 + pow(Tijm1, 3) / (6 * h));
            }
            if (!isnan(Eijp1))
            {
                Eip1jp1 = -(0.125 - pow(Tijp1, 3) / (6 * h)) * Eijp1 / (0.125 + pow(Tijp1, 3) / (6 * h));
            }
        }
        else
        {
            sigmaip1j = 1 / pow(Tip1j, 3);
            if (isnan(Eim2j))
            {
                Eim2j = (1 - (0.125 - 1 / (6 * sigmaim1j * h)) * Eim1j) / (0.125 + 1 / (6 * sigmaim1j * h));
            }
            else if (isnan(Eip2j))
            {
                Eip2j = -(0.125 - 1 / (6 * sigmaip1j * h)) * Eip1j / (0.125 + 1 / (6 * sigmaip1j * h));
            }
        }
    }
    
    if (isnan(Eijm1))
    {
        sigmaijm1 = sigmaij;
        sigmaijp1 = 1 / pow(Tijp1, 3);
        Eijm1 = Eij;
        Tijm1 = Tij;
        if (!isnan(Eim2j))
        {
            Eim1jm1 = Eim1j;
        }
        if (!isnan(Eip2j))
        {
            Eip1jm1 = Eip1j;
        }
    }
    else
    {
        sigmaijm1 = 1 / pow(Tijm1, 3);
        if (isnan(Eijp1))
        {
            sigmaijp1 = sigmaij;
            Eijp1 = Eij;
            Tijp1 = Tij;
            if (!isnan(Eim2j))
            {
                Eim1jp1 = Eim1j;
            }
            if (!isnan(Eip2j))
            {
                Eip1jp1 = Eip1j;
            }
        }
        else
        {
            sigmaijp1 = 1 / pow(Tijp1, 3);
            if (isnan(Eijm2))
            {
                Eijm2 = Eijm1;
            }
            else if(isnan(Eijp2))
            {
                Eijp2 = Eijp1;
            }
        }
    }
    
    double D1ij = 1 / (3 * sigmaij + sqrt(pow(Eip1j - Eim1j, 2) + pow(Eijp1 - Eijm1, 2)) / (2 * h * Eij));
    double D1im1j, D1ip1j, D1ijm1, D1ijp1;
    
    if (isnan(Eim2j))
    {
        if (isnan(Eim1jm1))
        {
            D1im1j = 1 / (3 * sigmaim1j + sqrt(pow(Eij - Eim1j, 2) + pow(Eim1jp1 - Eim1j, 2)) / (h * Eim1j));
        }
        else if (isnan(Eim1jp1))
        {
            D1im1j = 1 / (3 * sigmaim1j + sqrt(pow(Eij - Eim1j, 2) + pow(Eim1j - Eim1jm1, 2)) / (h * Eim1j));
        }
        else
        {
            D1im1j = 1 / (3 * sigmaim1j + sqrt(pow(Eij - Eim1j, 2) + pow((Eim1jp1 - Eim1jm1) / 2, 2)) / (h * Eim1j));
        }
    }
    else
    {
        D1im1j = 1 / (3 * sigmaim1j + sqrt(pow(Eij - Eim2j, 2) + pow(Eim1jp1 - Eim1jm1, 2)) / (2 * h * Eim1j));
    }
    
    if (isnan(Eip2j))
    {
        if (isnan(Eip1jm1))
        {
            D1ip1j = 1 / (3 * sigmaip1j + sqrt(pow(Eip1j - Eij, 2) + pow(Eip1jp1 - Eip1j, 2)) / (h * Eip1j));
        }
        else if (isnan(Eip1jp1))
        {
            D1ip1j = 1 / (3 * sigmaip1j + sqrt(pow(Eip1j - Eij, 2) + pow(Eip1j - Eip1jm1, 2)) / (h * Eip1j));
        }
        else
        {
            D1ip1j = 1 / (3 * sigmaip1j + sqrt(pow(Eip1j - Eij, 2) + pow((Eip1jp1 - Eip1jm1) / 2, 2)) / (h * Eip1j));
        }
    }
    else
    {
        D1ip1j = 1 / (3 * sigmaip1j + sqrt(pow(Eip2j - Eij, 2) + pow(Eip1jp1 - Eip1jm1, 2)) / (2 * h * Eip1j));
    }
    
    if (isnan(Eijm2))
    {
        if (isnan(Eim1jm1))
        {
            D1ijm1 = 1 / (3 * sigmaijm1 + sqrt(pow(Eip1jm1 - Eij, 2) + pow(Eij - Eijm1, 2)) / (h * Eijm1));
        }
        else if (isnan(Eip1jm1))
        {
            D1ijm1 = 1 / (3 * sigmaijm1 + sqrt(pow(Eij - Eim1jm1, 2) + pow(Eij - Eijm1, 2)) / (h * Eijm1));
        }
        else
        {
            D1ijm1 = 1 / (3 * sigmaijm1 + sqrt(pow((Eip1jm1 - Eim1jm1) / 2, 2) + pow(Eij - Eijm1, 2)) / (h * Eijm1));
        }
    }
    else
    {
        D1ijm1 = 1 / (3 * sigmaijm1 + sqrt(pow(Eip1jm1 - Eim1jm1, 2) + pow(Eij - Eijm2, 2)) / (2 * h * Eijm1));
    }

    if (isnan(Eijp2))
    {
        if (isnan(Eim1jp1))
        {
            D1ijp1 = 1 / (3 * sigmaijp1 + sqrt(pow(Eip1jp1 - Eim1j, 2) + pow(Eijp1 - Eij, 2)) / (h * Eijp1));
        }
        else if (isnan(Eip1jp1))
        {
            D1ijp1 = 1 / (3 * sigmaijp1 + sqrt(pow(Eip1j - Eim1jp1, 2) + pow(Eijp1 - Eij, 2)) / (h * Eijp1));
        }
        else
        {
            D1ijp1 = 1 / (3 * sigmaijp1 + sqrt(pow((Eip1jp1 - Eim1jp1) / 2, 2) + pow(Eijp1 - Eij, 2)) / (h * Eijp1));
        }
    }
    else
    {
        D1ijp1 = 1 / (3 * sigmaijp1 + sqrt(pow(Eip1jp1 - Eim1jp1, 2) + pow(Eijp2 - Eij, 2)) / (2 * h * Eijp1));
    }
    
    *fE = (
        (D1ij + D1ip1j) * (Eip1j - Eij) - (D1im1j + D1ij) * (Eij - Eim1j) +
        (D1ij + D1ijp1) * (Eijp1 - Eij) - (D1ijm1 + D1ij) * (Eij - Eijm1)
        ) / (2 * h * h) + sigmaij * (pow(Tij, 4) - Eij);
    *fT = k * (
        (pow(Tij, 2.5) + pow(Tip1j, 2.5)) * (Tip1j - Tij) - (pow(Tim1j, 2.5) + pow(Tij, 2.5)) * (Tij - Tim1j) +
        (pow(Tij, 2.5) + pow(Tijp1, 2.5)) * (Tijp1 - Tij) - (pow(Tijm1, 2.5) + pow(Tij, 2.5)) * (Tij - Tijm1)
        ) / (2 * h * h) - sigmaij * (pow(Tij, 4) - Eij);
}


/* right-hand-side function */
static void fraddiff2d(const unsigned* n, const double* x, const double* y, double* fy)
{
    const unsigned N = (unsigned)(0.5 + sqrt(*n / 2));
    const double h = 1. / N;

    unsigned i, j, jN = 0;

    i = 0;
    fij(h,
        y[jN + i], y[jN + i + 1], NAN, NAN, y[jN + i + 2], y[jN + i + 3],
        NAN, NAN, y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        NAN, y[jN + i + 4], NAN, y[jN + 4 * N + i],
        NAN, NAN, NAN, y[jN + 2 * N + i + 2],
        &fy[jN + i], &fy[jN + i + 1]);

    i = 2;
    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
        NAN, NAN, y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        NAN, y[jN + i + 4], NAN, y[jN + 4 * N + i],
        NAN, y[jN + 2 * N + i - 2], NAN, y[jN + 2 * N + i + 2],
        &fy[jN + i], &fy[jN + i + 1]);

    for (i = 4; i < 2 * N - 4; i += 2)
    {
        fij(h,
            y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
            NAN, NAN, y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
            y[jN + i - 4], y[jN + i + 4], NAN, y[jN + 4 * N + i],
            NAN, y[jN + 2 * N + i - 2], NAN, y[jN + 2 * N + i + 2],
            &fy[jN + i], &fy[jN + i + 1]);
    }

    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
        NAN, NAN, y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        y[jN + i - 4], NAN, NAN, y[jN + 4 * N + i],
        NAN, y[jN + 2 * N + i - 2], NAN, y[jN + 2 * N + i + 2],
        &fy[jN + i], &fy[jN + i + 1]);

    i += 2;
    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], NAN, NAN,
        NAN, NAN, y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        y[jN + i - 4], NAN, NAN, y[jN + 4 * N + i],
        NAN, y[jN + 2 * N + i - 2], NAN, NAN,
        &fy[jN + i], &fy[jN + i + 1]);

    jN = 2 * N;
    
    i = 0;
    fij(h,
        y[jN + i], y[jN + i + 1], NAN, NAN, y[jN + i + 2], y[jN + i + 3],
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        NAN, y[jN + i + 4], NAN, y[jN + 4 * N + i],
        NAN, NAN, y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
        &fy[jN + i], &fy[jN + i + 1]);

    i = 2;
    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        NAN, y[jN + i + 4], NAN, y[jN + 4 * N + i],
        y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
        &fy[jN + i], &fy[jN + i + 1]);

    for (i = 4; i < 2 * N - 4; i += 2)
    {
        fij(h,
            y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
            y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
            y[jN + i - 4], y[jN + i + 4], NAN, y[jN + 4 * N + i],
            y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
            &fy[jN + i], &fy[jN + i + 1]);
    }

    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        y[jN + i - 4], NAN, NAN, y[jN + 4 * N + i],
        y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
        &fy[jN + i], &fy[jN + i + 1]);

    i += 2;
    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], NAN, NAN,
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        y[jN + i - 4], NAN, NAN, y[jN + 4 * N + i],
        y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], NAN, NAN,
        &fy[jN + i], &fy[jN + i + 1]);

    for (j = 4; j < 2 * N - 4; j += 2)
    {
        jN = j * N;

        i = 0;
        fij(h,
            y[jN + i], y[jN + i + 1], NAN, NAN, y[jN + i + 2], y[jN + i + 3],
            y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
            NAN, y[jN + i + 4], y[jN - 4 * N + i], y[jN + 4 * N + i],
            NAN, NAN, y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
            &fy[jN + i], &fy[jN + i + 1]);

        i = 2;
        fij(h,
            y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
            y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
            NAN, y[jN + i + 4], y[jN - 4 * N + i], y[jN + 4 * N + i],
            y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
            &fy[jN + i], &fy[jN + i + 1]);

        for (i = 4; i < 2 * N - 4; i += 2)
        {
            fij(h, 
                y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
                y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
                y[jN + i - 4], y[jN + i + 4], y[jN - 4 * N + i], y[jN + 4 * N + i],
                y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
                &fy[jN + i], &fy[jN + i + 1]);
        }

        fij(h,
            y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
            y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
            y[jN + i - 4], NAN, y[jN - 4 * N + i], y[jN + 4 * N + i],
            y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
            &fy[jN + i], &fy[jN + i + 1]);

        i += 2;
        fij(h,
            y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], NAN, NAN,
            y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
            y[jN + i - 4], NAN, y[jN - 4 * N + i], y[jN + 4 * N + i],
            y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], NAN, NAN,
            &fy[jN + i], &fy[jN + i + 1]);
    }

    jN = j * N;

    i = 0;
    fij(h,
        y[jN + i], y[jN + i + 1], NAN, NAN, y[jN + i + 2], y[jN + i + 3],
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        NAN, y[jN + i + 4], y[jN - 4 * N + i], NAN,
        NAN, NAN, y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
        &fy[jN + i], &fy[jN + i + 1]);

    i = 2;
    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        NAN, y[jN + i + 4], y[jN - 4 * N + i], NAN,
        y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
        &fy[jN + i], &fy[jN + i + 1]);

    for (i = 4; i < 2 * N - 4; i += 2)
    {
        fij(h,
            y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
            y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
            y[jN + i - 4], y[jN + i + 4], y[jN - 4 * N + i], NAN,
            y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
            &fy[jN + i], &fy[jN + i + 1]);
    }

    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        y[jN + i - 4], NAN, y[jN - 4 * N + i], NAN,
        y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], y[jN - 2 * N + i + 2], y[jN + 2 * N + i + 2],
        &fy[jN + i], &fy[jN + i + 1]);

    i += 2;
    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], NAN, NAN,
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], y[jN + 2 * N + i], y[jN + 2 * N + i + 1],
        y[jN + i - 4], NAN, y[jN - 4 * N + i], NAN,
        y[jN - 2 * N + i - 2], y[jN + 2 * N + i - 2], NAN, NAN,
        &fy[jN + i], &fy[jN + i + 1]);

    j += 2;
    jN = j * N;

    i = 0;
    fij(h,
        y[jN + i], y[jN + i + 1], NAN, NAN, y[jN + i + 2], y[jN + i + 3],
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], NAN, NAN,
        NAN, y[jN + i + 4], y[jN - 4 * N + i], NAN,
        NAN, NAN, y[jN - 2 * N + i + 2], NAN,
        &fy[jN + i], &fy[jN + i + 1]);

    i = 2;
    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], NAN, NAN,
        NAN, y[jN + i + 4], y[jN - 4 * N + i], NAN,
        y[jN - 2 * N + i - 2], NAN, y[jN - 2 * N + i + 2], NAN,
        &fy[jN + i], &fy[jN + i + 1]);

    for (i = 4; i < 2 * N - 4; i += 2)
    {
        fij(h,
            y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
            y[jN - 2 * N + i], y[jN - 2 * N + i + 1], NAN, NAN,
            y[jN + i - 4], y[jN + i + 4], y[jN - 4 * N + i], NAN,
            y[jN - 2 * N + i - 2], NAN, y[jN - 2 * N + i + 2], NAN,
            &fy[jN + i], &fy[jN + i + 1]);
    }

    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], y[jN + i + 2], y[jN + i + 3],
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], NAN, NAN,
        y[jN + i - 4], NAN, y[jN - 4 * N + i], NAN,
        y[jN - 2 * N + i - 2], NAN, y[jN - 2 * N + i + 2], NAN,
        &fy[jN + i], &fy[jN + i + 1]);

    i += 2;
    fij(h,
        y[jN + i], y[jN + i + 1], y[jN + i - 2], y[jN + i - 1], NAN, NAN,
        y[jN - 2 * N + i], y[jN - 2 * N + i + 1], NAN, NAN,
        y[jN + i - 4], NAN, y[jN - 4 * N + i], NAN,
        y[jN - 2 * N + i - 2], NAN, NAN, NAN,
        &fy[jN + i], &fy[jN + i + 1]);
}


static void rhoraddiff2d(const unsigned* n, const double* x, const double* y, double* eigmax)
{
    const unsigned N = (unsigned)(0.5 + sqrt(*n / 2));
    const double h = 1. / N;

    *eigmax = 8 / (h * h) + 6000;
}


void get_raddiff2d(ProblemParams** params, FcnEqDiff* fcn, Rho* rho)
{
    *params = new RadDiff2dParams();

    *fcn = fraddiff2d;
    *rho = rhoraddiff2d;
}