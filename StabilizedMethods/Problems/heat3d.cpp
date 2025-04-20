#include "heat3d.h"


static void fheat3d(const unsigned* n, const double* x, const double* y, double* fy)
{
    const unsigned ns = (unsigned)(0.5 + pow(*n, 0.3333333333333333));
	const unsigned ns_sqr = ns * ns;
	const unsigned nsp1 = ns + 1;
	const double h = 1. / nsp1;

	double zz, yy, xx;
	double ykjip1, ykjim1, ykjp1i, ykjm1i, ykp1ji, ykm1ji;

	for (unsigned k = 0; k < ns; k++)
	{
		zz = (k + 1.) / nsp1;
		for (unsigned j = 0; j < ns; j++)
		{
			yy = (j + 1.) / nsp1;
			for (unsigned i = 0; i < ns; i++)
			{
				xx = (i + 1.) / nsp1;

				if (k == 0)
				{
					ykm1ji = tanh(5 * (xx + 2 * yy - *x - 0.5));
					ykp1ji = y[(k + 1) * ns_sqr + j * ns + i];
				}
				else
				{
					ykm1ji = y[(k - 1) * ns_sqr + j * ns + i];
					
					if (k == ns - 1)
					{
						ykp1ji = tanh(5 * (xx + 2 * yy + 1.5 - *x - 0.5));
					}
					else
					{
						ykp1ji = y[(k + 1) * ns_sqr + j * ns + i];
					}
				}

				if (j == 0)
				{
					ykjm1i = tanh(5 * (xx + 1.5 * zz - *x - 0.5));
					ykjp1i = y[k * ns_sqr + (j + 1) * ns + i];
				}
				else
				{
					ykjm1i = y[k * ns_sqr + (j - 1) * ns + i];

					if (j == ns - 1)
					{
						ykjp1i = tanh(5 * (xx + 2. + 1.5 * zz - *x - 0.5));
					}
					else
					{
						ykjp1i = y[k * ns_sqr + (j + 1) * ns + i];
					}
				}

				if (i == 0)
				{
					ykjim1 = tanh(5 * (2 * yy + 1.5 * zz - *x - 0.5));
					ykjip1 = y[k * ns_sqr + j * ns + i + 1];
				}
				else
				{
					ykjim1 = y[k * ns_sqr + j * ns + i - 1];

					if (i == ns - 1)
					{
						ykjip1 = tanh(5 * (1 + 2 * yy + 1.5 * zz - *x - 0.5));
					}
					else
					{
						ykjip1 = y[k * ns_sqr + j * ns + i + 1];
					}
				}

				fy[k * ns_sqr + j * ns + i] = 
					(
						ykjip1 - 2 * y[k * ns_sqr + j * ns + i] + ykjim1 + 
						ykjp1i - 2 * y[k * ns_sqr + j * ns + i] + ykjm1i +
						ykp1ji - 2 * y[k * ns_sqr + j * ns + i] + ykm1ji
					) / (h * h) 
					+
					pow(cosh(5 * (xx + 2 * yy + 1.5 * zz - *x - 0.5)), -2) * (-5 + 362.5 * tanh(5 * (xx + 2 * yy + 1.5 * zz - *x - 0.5)));
			}
		}
	}
}


static void rhoheat3d(const unsigned* n, const double* x, const double* y, double* eigmax)
{
	const unsigned ns = (unsigned)(0.5 + pow(*n, 0.3333333333333333));
	const double h = 1. / (ns + 1);

	*eigmax = 12 / (h * h);
}


void get_heat3d(ProblemParams** params, FcnEqDiff* fcn, Rho* rho)
{
    *params = new Heat3dParams();

    *fcn = fheat3d;
	*rho = rhoheat3d;
}