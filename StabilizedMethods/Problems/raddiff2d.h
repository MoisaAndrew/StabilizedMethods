#include "../problems.h"


static const double k = 5.0e-3;


struct RadDiff2dParams : ProblemParams
{
	RadDiff2dParams()
	{
		nDefault = 20000;

		x0 = 0;
		h0 = 1e-6;
		xend = 3;

		isRhoDefined = true;
		isJacConst = true;
	}

	double* y0(const unsigned n) const
	{
		double* y0 = (double*)malloc(n * sizeof(double));
		for (unsigned i = 0; i < n; i += 2)
		{
			y0[i] = 1.0e-5;
			y0[i + 1] = 0.05623413251903491;
		}
		return y0;
	}
};