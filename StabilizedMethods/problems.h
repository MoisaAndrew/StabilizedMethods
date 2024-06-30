#include "methods_common.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>


void get_bruss(unsigned* n, FcnEqDiff* fcn, Rho* rho, double* x, double* h, double* xend, double** y);

void get_bruss2d(unsigned* n, FcnEqDiff* fcn, Rho* rho, double* x, double* h, double* xend, double** y);

void get_burgers(unsigned* n, FcnEqDiff* fcn, Rho* rho, double* x, double* h, double* xend, double** y);

void get_comb(unsigned* n, FcnEqDiff* fcn, Rho* rho, double* x, double* h, double* xend, double** y);

void get_cusp(unsigned* n, FcnEqDiff* fcn, Rho* rho, double* x, double* h, double* xend, double** y);

void get_finag(unsigned* n, FcnEqDiff* fcn, Rho* rho, double* x, double* h, double* xend, double** y);

void get_hires(unsigned* n, FcnEqDiff* fcn, Rho* rho, double* x, double* h, double* xend, double** y);

void get_rober(unsigned* n, FcnEqDiff* fcn, Rho* rho, double* x, double* h, double* xend, double** y);

void get_heat(unsigned* n, FcnEqDiff* fcn, Rho* rho, double* x, double* h, double* xend, double** y);

void get_raddiff(unsigned* n, FcnEqDiff* fcn, Rho* rho, double* x, double* h, double* xend, double** y);

void get_raddiff2d(unsigned* n, FcnEqDiff* fcn, Rho* rho, double* x, double* h, double* xend, double** y);