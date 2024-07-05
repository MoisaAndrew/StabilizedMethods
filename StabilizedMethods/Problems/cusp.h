#include "../problems.h"


static const double sigma = 1. / 144.;


struct CuspParams : ProblemParams
{
	CuspParams()
	{
		nDefault = 96;

		x0 = 0;
		h0 = 1e-8;
		xend = 1.1;

		isRhoDefined = false;
		isJacConst = false;
	}

	double* y0(const unsigned n) const
	{
		const unsigned N = n / 3;
		
		double* y0 = (double*)malloc(n * sizeof(double));
		for (unsigned j = 0; j < N; j++) {
			y0[3 * j] = 0;
			y0[3 * j + 1] = -2 * cos(2 * j * M_PI / N);
			y0[3 * j + 2] = 2 * sin(2 * j * M_PI / N);
		}
		return y0;
	}
};