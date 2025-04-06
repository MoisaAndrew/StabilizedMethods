#include "../problems.h"


static const double L = 0.9;
static const double R = 5.0;
static const double alpha = 1.0;
static const double delta = 20.0;

static const double D = R * exp(delta) / (alpha * delta);


struct Comb3dParams : ProblemParams
{
	Comb3dParams()
	{
		nDefault = 128000;

		x0 = 0;
		h0 = 1e-7;
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
};