#include "methods_common.h"


int rkcc(const unsigned n, const double x, const double xend, double* y,
	const FcnEqDiff f, const Rho rho, const SolTrait solout,
	const double* atol, const double rtol,
	unsigned iwork[10]);