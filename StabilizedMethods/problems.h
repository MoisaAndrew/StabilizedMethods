#include "methods_common.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>


struct ProblemParams
{
	unsigned nDefault;
	unsigned nFine; // for discretization errors computation

	double x0;
	double h0;
	double xend;

	bool isRhoDefined;
	bool isSpradConst;

	virtual double* y0(const unsigned n) const { return NULL; };

	virtual double* y_exact(const unsigned n) const { return NULL; };

	virtual double* remap_solution(const unsigned nFrom, const unsigned nTo, const double* y) const { return NULL; };
};


void get_bruss(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_bruss2d(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_burgers(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_comb2d(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_comb3d(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_cusp(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_finag(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_heat(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_heat3d(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_hires(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_raddiff(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_raddiff2d(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);

void get_rober(ProblemParams** params, Fcn* fcn, Rho* rho, Jac* jac, PSol* psol);


static void psoldiag(const unsigned* n, const double* x, const double* y,
	const double* f, double* wk,
	const double* hrl1, double* p, int* iwp,
	double* b, const unsigned* lr, int* idid)
{
	for (unsigned i = 0; i < *n; i++)
	{
		b[i] /= p[i];
	}
}