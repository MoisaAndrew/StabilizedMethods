//	Author : A.Moisa
//		e - mail : andrey.moysa@gmail.com
//
//	Version of September 2024
//
//	A.V. Moisa, A family of two-step second order Runge�Kutta�Chebyshev methods
//	Journal of Computational and Applied Mathematics 446 (2024)
//	https://doi.org/10.1016/j.cam.2024.115868


#include "methods_common.h"


int tsrkc2trho(const unsigned n, const double x,
	const double* yn, const double* fn,
	const double uround,
	const FcnEqDiff f,
	double* eigvec, double* z, double* fz,
	double* eigmax,
	unsigned iwork[12])
{
	const unsigned maxiter = 50;
	const double safe = 1.2;

	const double sqrtu = sqrt(uround);
	double ynor = 0, znor = 0, dzyn, quot, dfzfn, eigmaxo;
	unsigned ind, nind, ntest, i, iter;

	if (iwork[5] == 0)
	{
		memcpy(fz, fn, n * sizeof(double));
		f(&n, &x, fz, z);
		iwork[8]++;
	}
	else
	{
		memcpy(z, eigvec, n * sizeof(double));
	}

	for (i = 0; i < n; i++)
	{
		ynor += yn[i] * yn[i];
		znor += z[i] * z[i];
	}
	ynor = sqrt(ynor);
	znor = sqrt(znor);

	if (ynor != 0)
	{
		if (znor != 0)
		{
			dzyn = ynor * sqrtu;
			quot = dzyn / znor;
			for (i = 0; i < n; i++)
			{
				z[i] = yn[i] + z[i] * quot;
			}
		}
		else
		{
			dzyn = ynor * sqrtu;
			for (i = 0; i < n; i++)
			{
				z[i] = yn[i] + yn[i] * sqrtu;
			}
		}
	}
	else
	{
		if (znor != 0)
		{
			quot = uround / znor;
			for (i = 0; i < n; i++)
			{
				z[i] *= quot;
			}
		}
		else
		{
			for (i = 0; i < n; i++)
			{
				z[i] = uround;
			}
		}
	}

	*eigmax = 0;
	for (iter = 1; iter <= maxiter; iter++)
	{
		f(&n, &x, z, fz);
		iwork[8]++;
		dfzfn = 0;
		for (i = 0; i < n; i++)
		{
			dfzfn += pow(fz[i] - fn[i], 2);
		}
		dfzfn = sqrt(dfzfn);
		eigmaxo = *eigmax;
		*eigmax = dfzfn / dzyn;
		*eigmax *= safe;

		if (iter >= 2 && fabs(*eigmax - eigmaxo) <= *eigmax * 0.05)
		{
			for (i = 0; i < n; i++)
			{
				eigvec[i] = z[i] - yn[i];
			}
			return 1;
		}

		if (dfzfn != 0)
		{
			quot = dzyn / dfzfn;
			for (i = 0; i < n; i++)
			{
				z[i] = yn[i] + (fz[i] - fn[i]) * quot;
			}
		}
		else
		{
			nind = n;
			ntest = 0;
			ind = iter % nind;
			if (z[ind] != yn[ind] || ntest == 10)
			{
				z[ind] = yn[ind] - (z[ind] - yn[ind]);
			}
			else
			{
				nind = n + ind;
				ntest++;
			}
		}
	}

	printf("convergence failure in the spectral radius computation");
	return -3;
}


void tsrkcstep(const unsigned n,
	const double x0, const double x1, const double h,
	const double* y0, const double* y1, const double* f1,
	const unsigned mdeg,
	const double w, const double beta, const double gamma, const double acosht,
	const FcnEqDiff f,
	double* y2, double* yjm1, double* yjm2)
{
	double* yswap;
	double* res = y2;
	const double hbeta = h * beta;
	const double acoshtds = acosht / mdeg;
	const double onemgamma = 1 - gamma;
	double coshims = 1, coshi = w, tim1wdtiw;
	double temp1 = hbeta / w, temp2, temp3;
	double ci1 = x1 + temp1, ci2 = x1 + temp1, ci3 = x1;
	double err = 0;
	unsigned i, j;

	memcpy(yjm2, y1, n * sizeof(double));
	for (i = 0; i < n; i++)
	{
		yjm1[i] = yjm2[i] + temp1 * f1[i];
	}

	for (i = 2; i <= mdeg; i++)
	{
		coshims = coshi;
		coshi = cosh(i * acoshtds);
		tim1wdtiw = coshims / coshi;
		temp1 = 2 * hbeta * tim1wdtiw;
		temp2 = 2 * w * tim1wdtiw;
		temp3 = 1 - temp2;
		f(&n, &ci1, yjm1, y2);
		ci1 = temp1 + temp2 * ci2 + temp3 * ci3;
		for (j = 0; j < n; j++)
		{
			y2[j] = temp1 * y2[j] + temp2 * yjm1[j] + temp3 * yjm2[j];
		}

		yswap = yjm2;
		yjm2 = yjm1;
		yjm1 = y2;
		y2 = yswap;

		ci3 = ci2;
		ci2 = ci1;
	}

	for (j = 0; j < n; j++)
	{
		res[j] = gamma * y0[j] + onemgamma * yjm1[j];
	}
}


int tsrkc2core(const unsigned n,
	double x0, double x1, const double xend,
	double* h,
	double* y0, double* y1,
	const FcnEqDiff f, const Rho rho, const SolTrait solout,
	const double* atol, const double* rtol,
	const double uround, unsigned iwork[12])
{
	double facmax = 5, xold = x1, hp = x1 - x0, err = 0;
	double errp = 0, hnew, fac, facp, eigmax;
	unsigned nrej = 0, mdeg, mdego = 0, nrho = 0;
	const double tsw = 1.1, acosht = acosh(1.1);
	double s, s2dls;
	double w, w2m1, dtsw, d2tsw;
	double q, onemq;
	double beta, alpha, gamma;
	bool last = false, reject = false;
	double* y2 = (double*)malloc(n * sizeof(double));
	double* swap;
	double *res[3] = { y0, y1, y2 };
	double* f1 = (double*)malloc(n * sizeof(double));
	double* f2 = (double*)malloc(n * sizeof(double));
	double* yjm1 = (double*)malloc(n * sizeof(double));
	double* yjm2 = (double*)malloc(n * sizeof(double));
	double* eigvec = (double*)malloc(n * sizeof(double));
	unsigned j;
	double temp1, temp2;
	int idid = 1;

	f(&n, &x1, y1, f1);
	iwork[4]++;

	while (true)
	{
		if (1.1 * (*h) >= fabs(xend - x1))
		{
			*h = fabs(xend - x1);
			last = true;
		}
		if (x1 + 10 * *h <= x1)
		{
			printf("stepsize becomes too small\r\n");
			idid = -2;
			break;
		}

		if (nrho == 0)
		{
			if (iwork[5] == 0 || iwork[1] == 0)
			{
				if (iwork[0] == 1)
				{
					rho(&n, &x1, y1, &eigmax);
				}
				else
				{
					idid = tsrkc2trho(n, x1, y1, f1, uround, f, eigvec, yjm1, yjm2, &eigmax, iwork);
					if (idid == -3)
					{
						return -3;
					}
				}

				if (eigmax > iwork[10])
				{
					iwork[10] = (unsigned)(eigmax + 0.5);
				}
				if (iwork[5] == 0)
				{
					iwork[11] = iwork[10];
				}
				if (eigmax < iwork[11])
				{
					iwork[11] = (unsigned)(eigmax + 0.5);
				}
			}
		}

		q = hp / *h;
		onemq = 1 - q;

		s2dls = 0.759782816506459 * (onemq + sqrt(1 + q * (q - 0.598626091572911)));
		mdeg = 1 + sqrt(*h * eigmax * s2dls);
		if (mdeg < 2)
		{
			mdeg = 2;
		}
		if (mdeg != mdego)
		{
			s = (double)mdeg;
			w = cosh(acosht / s);
			w2m1 = w * w - 1;
			dtsw = s * sinh(acosht) / sqrt(w2m1);
			d2tsw = (s * s * tsw - w * dtsw) / w2m1;
		}

		beta = (onemq * dtsw + sqrt(pow(onemq * dtsw, 2) + 4 * q * tsw * d2tsw)) / (2 * d2tsw);
		alpha = (1 + q) / (q * tsw + beta * dtsw);
		gamma = 1 - alpha * tsw;

		if (mdeg > iwork[9])
		{
			iwork[9] = mdeg;
		}

		tsrkcstep(n, x0, x1, *h, y0, y1, f1, mdeg, w, beta, gamma, acosht, f, y2, yjm1, yjm2);

		temp2 = x1 + *h;
		f(&n, &temp2, y2, f2);

		err = 0;
		if (iwork[3] == 0)
		{
			for (j = 0; j < n; j++)
			{
				temp1 = y1[j] - y2[j] + *h * f2[j];
				temp2 = fmax(fabs(y1[j]), fabs(y2[j])) * (*rtol);
				err += pow(temp1 / (*atol + temp2), 2);
			}
		}
		else
		{
			for (j = 0; j < n; j++)
			{
				temp1 = y1[j] - y2[j] + *h * f2[j];
				temp2 = fmax(fabs(y1[j]), fabs(y2[j])) * rtol[j];
				err += pow(temp1 / (atol[j] + temp2), 2);
			}
		}

		err = sqrt(err / n);
		err *= 0.3333333333333333;

		mdego = mdeg;
		iwork[5]++;
		iwork[4] += mdeg;
		
		if (isfinited(err))
		{
			fac = sqrt(1 / err);
			if (errp != 0 && !reject)
			{
				facp = sqrt(errp) / (q * err);
				fac = fmin(fac, facp);
			}
			if (reject)
			{
				facmax = 1;
			}
			fac = fmin(facmax, fmax(0.1, 0.8 * fac));
			hnew = *h * fac;
		}
		else
		{
			hnew = *h * 0.1;
			err = 2;
		}

		if (err < 1)
		{
			if (iwork[2] == 1)
			{
				solout(n, x1, x1 + (*h), y2);
			}

			iwork[6]++;
			facmax = 2;
			x0 = x1;
			x1 += *h;
			if (reject)
			{
				hnew = fmin(hnew, *h);
				if (xend < x1)
				{
					hnew = fmax(hnew, *h);
				}
				reject = false;
				nrej = 0;
			}
			hp = *h;
			*h = hnew;
			nrho = (nrho + 1) % 25;
			if (last)
			{
				idid = 1;
				break;
			}
			else
			{
				swap = y0;
				y0 = y1;
				y1 = y2;
				y2 = swap;

				swap = f1;
				f1 = f2;
				f2 = swap;

				errp = err;
			}
		}
		else
		{
			iwork[7]++;
			reject = true;
			last = false;
			*h = 0.8 * hnew;
			if (xold == x1)
			{
				nrej++;
				if (nrej == 10)
				{
					*h = 1e-5;
				}
			}
			xold = x1;

			if (nrho != 0)
			{
				nrho = 0;
			}
			else
			{
				nrho = 13;
			}
		}
	}

	if (res[1] != y2)
	{
		if (res[0] != y2)
		{
			memcpy(res[0], y1, n * sizeof(double));
			memcpy(res[1], y2, n * sizeof(double));
		}
		else
		{
			memcpy(res[1], y2, n * sizeof(double));
			memcpy(res[0], y1, n * sizeof(double));
		}
	}

	free(res[2]);
	free(f1);
	free(f2);
	free(yjm1);
	free(yjm2);
	free(eigvec);

	return idid;
}


/// <summary>
/// Two-Step RKC method of order 2
/// </summary>
/// <param name="n">Number of differential equations of the system</param>
/// <param name="x0">Initial point of integration</param>
/// <param name="x1">Additional point of integration. Set x1 = x0 if solution is unknown elsewhere</param>
/// <param name="xend">End of the interval of integration, may be less than x0</param>
/// <param name="h">Initial step size guess (usually between 1e-4 and 1e-6)</param>
/// <param name="y0">Initial value of the solution (array of length n)</param>
/// <param name="y1">Additional value of the solution (array of length n)</param>
/// <param name="f">Name (external) of function computing the value of f(x, y)</param>
/// <param name="rho">Name (external) of a function giving the spectral radius of the Jacobian matrix.
/// Supply a NULL pointer if iwork[0] == 1, and TSRKC2 will compute spectral radius internally</param>
/// <param name="solout">Name (external) of a function providing the numerical solution during integration. 
/// Supply a NULL pointer if iwork[2] = 0.</param>
/// <param name="atol">Absolute error tolerance. Can be scalar or vector of length n</param>
/// <param name="rtol">Relative error tolerance. Can be scalar or vector of length n</param>
/// <param name="iwork">Integer array of length 12 that 
/// gives information on how the problem is to be solved and 
/// communicates statistics about the integration process.
/// <para> iwork[0] </para>
///	<para> = 0 TSRKC2 attempts to compute the spectral radius internally; </para>
///	<para> = 1 rho returns an upper bound of the spectral radius of the Jacobian matrix </para>
/// <para> iwork[1] </para>
///	<para> = 0 The Jacobian is not constant; </para>
///	<para> = 1 The Jacobian is constant </para>
/// <para> iwork[2] </para>
///	<para> = 0 function solout is called after every successful step </para>
///	<para> = 1 function solout is never called </para>
/// <para> iwork[3] </para>
///	<para> = 0 Atol and rtol are scalar </para>
///	<para> = 1 Atol and rtol are arrays of length n </para>
/// </param>
/// <returns>
/// <para> 1 Successful computation x = xend; </para>
/// <para> -1 Invalid input parameters; </para>
/// <para> -2 Stepsize becomes too small; </para>
/// <para> -3 The method used in TSRKC2 to estimate the spectral radius did not converge. </para>
/// <para> iwork[4] - Number of function evaluations </para>
/// <para> iwork[5] - Number of steps </para>
/// <para> iwork[6] - Number of accepted steps </para>
/// <para> iwork[7] - Number of rejected steps </para>
/// <para> iwork[8] - Number of evaluations of f used to estimate the spectral radius (equal to zero if iwork[0] = 1) </para>
/// <para> iwork[9] - Maximum number of stages used </para>
/// <para> iwork[10] - Maximum value of the estimated bound for the spectral radius </para>
/// <para> iwork[11] - Minimum value of the estimated bound for the spectral radius </para>
/// </returns>
int tsrkc2(const unsigned n,
	const double x0, const double x1, const double xend,
	double* h,
	double* y0, double* y1,
	const FcnEqDiff f, const Rho rho, const SolTrait solout,
	const double* atol, const double* rtol,
	unsigned iwork[12])
{
	// may depends on the machines
	const double uround = 1e-16;

	if (*h > fabs(xend - x1))
	{
		printf("initial step is longer than the integration interval\r\n");
		return -1;
	}
	if (*h < 10 * uround)
	{
		printf("initial step-size is too small\r\n");
		return -1;
	}
	if (iwork[3] == 0)
	{
		if (atol[0] <= 0 || rtol[0] <= 10 * uround)
		{
			printf("tolerances are too small\r\n");
			return -1;
		}
	}
	else
	{
		for (unsigned i = 0; i < n; i++)
		{
			if (atol[i] <= 0 || rtol[i] <= 10 * uround)
			{
				printf("tolerances are too small\r\n");
				return -1;
			}
		}
	}

	iwork[4] = 0, iwork[5] = 0, iwork[6] = 0, iwork[7] = 0,
	iwork[8] = 0, iwork[9] = 0, iwork[10] = 0, iwork[11] = 0;

	return tsrkc2core(n, x0, x1, xend, h, y0, y1, f, rho, solout, atol, rtol, uround, iwork);
}