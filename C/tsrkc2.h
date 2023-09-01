#include "methods_common.h"


int tsrkc2(const unsigned n, 
	const double x0, const double x1, const double xend, 
	double* h, 
	double* y0, double* y1,
	const FcnEqDiff f, const Rho rho, const SolTrait solout,
	const double* atol, const double* rtol,
	unsigned iwork[12]);