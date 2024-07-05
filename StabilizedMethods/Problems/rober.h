#include "../problems.h"


struct RoberParams : ProblemParams
{
	RoberParams()
	{
		nDefault = 3;

		x0 = 0;
		h0 = 1e-4;
		xend = 1.0e6;

		isRhoDefined = false;
		isJacConst = false;
	}

	double* y0(const unsigned n) const
	{
		double* y0 = (double*)malloc(n * sizeof(double));
		y0[0] = 1; y0[1] = 0; y0[2] = 0;
		return y0;
	}
};