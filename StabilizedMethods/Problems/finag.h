#include "../problems.h"


static const double alpha = 0.139;
static const double beta = 2.54;
static const double eta = 0.008;


struct FinagParams : ProblemParams
{
	FinagParams()
	{
		nDefault = 400;

		x0 = 0;
		h0 = 1e-4;
		xend = 400;

		isRhoDefined = false;
		isJacConst = false;
	}

	double* y0(const unsigned n) const
	{
		double* y0 = (double*)malloc(n * sizeof(double));
		for (unsigned i = 0; i < n; i++)
		{
			y0[i] = 0;
		}
		return y0;
	}
};