#include "../problems.h"


static const double L = 0.9;
static const double R = 5.0;
static const double alpha = 1.0;
static const double delta = 20.0;


struct Comb3dParams : ProblemParams
{
	Comb3dParams()
	{
		nDefault = 2 * 40 * 40 * 40;
		nFine = 2 * 364 * 364 * 364;

		x0 = 0;
		h0 = 1e-5;
		xend = 0.3;

		isRhoDefined = false;
		isJacConst = false;
	}

	double* y0(const unsigned n) const
	{
		double* y0 = (double*)malloc(n * sizeof(double));
		for (unsigned i = 0; i < n; i++)
		{
			y0[i] = 1.0;
		}
		return y0;
	}

	double* remap_solution(const unsigned nFrom, const unsigned nTo, const double* y) const
	{
		const unsigned N_from = (unsigned)(pow(0.5 * nFrom, 0.3333333333333333) + 0.5);
		const unsigned N_to = (unsigned)(pow(0.5 * nTo, 0.3333333333333333) + 0.5);
		const unsigned N_from_sqr = N_from * N_from;
		const unsigned N_to_sqr = N_to * N_to;

		const unsigned rate = (unsigned)((N_from - 1.) / N_to + 0.5);
		const unsigned offset = (rate - 1) / 2;

		double* res = (double*)malloc(nTo * sizeof(double));
		for (unsigned k = 0; k < N_to; k++)
		{
			for (unsigned j = 0; j < N_to; j++)
			{
				for (unsigned i = 0; i < N_to; i++)
				{
					res[k * 2 * N_to_sqr + j * 2 * N_to + 2 * i] =
						y[(k * rate + offset) * 2 * N_from_sqr + (j * rate + offset) * 2 * N_from + 2 * (i * rate + offset)];
					res[k * 2 * N_to_sqr + j * 2 * N_to + 2 * i + 1] =
						y[(k * rate + offset) * 2 * N_from_sqr + (j * rate + offset) * 2 * N_from + 2 * (i * rate + offset) + 1];
				}
			}
		}

		return res;
	}
};