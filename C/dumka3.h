#include "methods_common.h"


int dumka3(const unsigned neqn, double* time, const double tend, const double h0,
	const double atol, const double rtol,
	const FcnEqDiff f, const Rho cour,
	double* y, double* z0, double* z1, double* z2,
	double* z3, double* oldEigenVector,
	unsigned iwork[12]);
