#include "../problems.h"


static const double mu = 0.0003;


struct BurgersParams : ProblemParams
{
	BurgersParams()
	{
		nDefault = 500;

		x0 = 0;
		h0 = 1e-4;
		xend = 2.5;

		isRhoDefined = false;
		isJacConst = false;
	}

	double* y0(const unsigned n) const
	{
		double* y0 = (double*)malloc(n * sizeof(double));
		for (unsigned i = 0; i < n; i++)
		{
			double temp = (i + 1) / (n + 1.);
			y0[i] = 1.5 * temp * (1 - temp) * (1 - temp);
		}
		return y0;
	}
};