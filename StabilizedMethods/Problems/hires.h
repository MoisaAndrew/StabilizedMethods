#include "../problems.h"


struct HiresParams : ProblemParams
{
	HiresParams()
	{
		nDefault = 8;

		x0 = 0;
		h0 = 1e-5;
		xend = 421.8122;

		isRhoDefined = false;
		isJacConst = false;
	}

	double* y0(const unsigned n) const
	{
		double* y0 = (double*)malloc(n * sizeof(double));
		y0[0] = 1.0; y0[1] = 0.0; y0[2] = 0.0; y0[3] = 0.0;
		y0[4] = 0.0; y0[5] = 0.0; y0[6] = 0.0; y0[7] = 0.0057;
		return y0;
	}
};