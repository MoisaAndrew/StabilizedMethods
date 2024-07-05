#include "../problems.h"


static const double mu = 0.0003;


struct HeatParams : ProblemParams
{
	HeatParams()
	{
		nDefault = 159;

		x0 = 0;
		h0 = 1e-5; //1 / 20.;
		xend = 1;

		isRhoDefined = false;
		isJacConst = false;
	}

	double* y0(const unsigned n) const
	{
		double* y0 = (double*)malloc(n * sizeof(double));
		for (unsigned i = 0; i < n; i++)
		{
			double x = (i + 1) / (n + 1.);
			y0[i] = x * (x + 1);
		}
		return y0;
	}
};