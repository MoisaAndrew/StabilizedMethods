#include "../problems.h"


static const double alpha = 0.1;


struct Bruss2dParams : ProblemParams
{
	Bruss2dParams()
	{
		nDefault = 32768;

		x0 = 0;
		h0 = 1e-7;
		xend = 11.5;

		isRhoDefined = false;
		isJacConst = false;
	}

	double* y0(const unsigned n) const
	{
		const unsigned N_sqr = n / 2;
		const double ans = sqrt(N_sqr);
		const unsigned N = (unsigned)(ans + 0.5);

		double* y0 = (double*)malloc(n * sizeof(double));
		for (unsigned j = 0; j < N; j++) {
			const double yy = j / ans;
			for (unsigned i = 0; i < N; i++) {
				y0[j * N + i] = 22.0 * yy * pow(1.0 - yy, 1.5);
			}
		}
		for (unsigned i = 0; i < N; i++) {
			const double xx = i / ans;
			for (unsigned j = 0; j < N; j++) {
				y0[N_sqr + j * N + i] = 27.0 * xx * pow(1.0 - xx, 1.5);
			}
		}
		return y0;
	}
};