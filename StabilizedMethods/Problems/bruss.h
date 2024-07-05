#include "../problems.h"


static const double alpha = 0.02;


struct BrussParams : ProblemParams
{
	BrussParams()
	{
		nDefault = 1000;

		x0 = 0;
		h0 = 1e-5;
		xend = 10;

		isRhoDefined = false;
		isJacConst = false;
	}

	double* y0(const unsigned n) const
	{
		const unsigned N = n / 2;
		double* y0 = (double*)malloc(n * sizeof(double));
		for (unsigned i = 1; i <= N; i++)
		{
			y0[2 * i - 2] = 1.0 + sin(2 * M_PI * i / (N + 1));
			y0[2 * i - 1] = 3.0;
		}
		return y0;
	}
};