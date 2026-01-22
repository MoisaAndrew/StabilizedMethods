#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef void (*Fcn)(const unsigned* n, const double* x, const double* y, double* f);

typedef void (*Rho)(const unsigned* n, const double* x, const double* y, double* eigmax);

typedef void (*Jac)
(
	Fcn fcn, const unsigned* n, const double* x, const double* y,
	const double* yy, const double* rewt, const double* f, const double* ff,
	const double* hrl1, double* wp, int* iwp, int* idid
);

typedef void (*PSol)
(
	const unsigned* n, const double* x, const double* y,
	const double* f, double* wk,
	const double* hrl1, double* p, int* iwp,
	double* b, const unsigned* lr, int* idid
);

typedef void (*SolTrait)(const unsigned n, const double xold, const double x, const double* y);