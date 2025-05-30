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
	bool isJacConst;

	virtual double* y0(const unsigned n) const { return NULL; };

	virtual double* y_exact(const unsigned n) const { return NULL; };

	virtual double* remap_solution(const unsigned nFrom, const unsigned nTo, const double* y) const { return NULL; };
};


void get_bruss(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_bruss2d(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_burgers(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_comb2d(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_comb3d(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_cusp(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_finag(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_heat(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_heat3d(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_hires(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_raddiff(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_raddiff2d(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);

void get_rober(ProblemParams** params, FcnEqDiff* fcn, Rho* rho);