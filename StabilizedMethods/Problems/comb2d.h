#include "../problems.h"


static const double d = 1.0;
static const double R = 5.0;
static const double alpha = 1.0;
static const double delta = 20.0;


struct Comb2dParams : ProblemParams
{
	Comb2dParams()
	{
		nDefault = 6400;

		x0 = 0;
		h0 = 1e-7;
		xend = 0.32;

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
};