#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef void (*FcnEqDiff)(const unsigned* n, const double* x, const double* y, double* f);
typedef void (*Rho)(const unsigned* n, const double* x, const double* y, double* eigmax);
typedef void (*SolTrait)(const unsigned n, const double xold, const double x, const double* y);