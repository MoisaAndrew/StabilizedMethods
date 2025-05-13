#include "../problems.h"


struct Heat3dParams : ProblemParams
{
	Heat3dParams()
	{
		nDefault = 39 * 39 * 39;

		x0 = 0;
		h0 = 1e-5;
		xend = 0.7;

		isRhoDefined = true;
		isJacConst = true;
	}

	double* y0(const unsigned n) const
	{
		const double ans = pow(n, 0.3333333333333333);
		const unsigned N = (unsigned)(ans + 0.5);
		const unsigned N_sqr = N * N;

		double* y0 = (double*)malloc(n * sizeof(double));
		for (unsigned k = 0; k < N; k++)
		{
			double zz = (k + 1) / (N + 1.);
			for (unsigned j = 0; j < N; j++)
			{
				double yy = (j + 1) / (N + 1.);
				for (unsigned i = 0; i < N; i++)
				{
					double xx = (i + 1) / (N + 1.);
					y0[k * N_sqr + j * N + i] = tanh(5 * (xx + 2 * yy + 1.5 * zz - 0.5));
				}
			}
		}
		return y0;
	}

	double* y_exact(const unsigned n) const
	{
		const double ans = pow(n, 0.3333333333333333);
		const unsigned N = (unsigned)(ans + 0.5);
		const unsigned N_sqr = N * N;

		double* y = (double*)malloc(n * sizeof(double));
		for (unsigned k = 0; k < N; k++)
		{
			double zz = (k + 1) / (N + 1.);
			for (unsigned j = 0; j < N; j++)
			{
				double yy = (j + 1) / (N + 1.);
				for (unsigned i = 0; i < N; i++)
				{
					double xx = (i + 1) / (N + 1.);
					y[k * N_sqr + j * N + i] = tanh(5 * (xx + 2 * yy + 1.5 * zz - 0.5 - 0.7));
				}
			}
		}
		return y;
	}
};