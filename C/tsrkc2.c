#include "tsrkc2.h"


int tsrkctrho(const unsigned n, const double x,
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


double tsrkcstep(const unsigned n, 
	const double x0, const double x1, const double h, 
	const double* y0, const double* y1, const double* f1, 
	const double* atol, const double* rtol, const int ntol, 
	const double mdeg,
	const double w, const double beta, const double gamma, const double gammaemb, const double acosht,
	const FcnEqDiff f, 
	double* y2, double* yjm1, double* yjm2)
{
	double* yswap;
	double *res = y2;
	const double s = (double)mdeg, hbetads2 = h * beta / (s * s);
	const double acoshtds = acosht / s;
	const double onemgamma = 1 - gamma, onemgammaemb = 1 - gammaemb;
	double coshims, coshi = cosh(acoshtds), tim1wdtiw;
	double temp1 = hbetads2 / w, temp2, temp3;
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
		temp1 = 2 * hbetads2 * tim1wdtiw;
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
		y2[j] = gamma * y0[j] + onemgamma * yjm1[j];
		yjm1[j] = gammaemb * y0[j] + onemgammaemb * yjm2[j];
	}

	if (ntol == 0)
	{
		for (j = 0; j < n; j++)
		{
			temp3 = y2[j] - yjm1[j];
			ci1 = fmax(fabs(y1[j]), fabs(y2[j])) * (*rtol);
			err += pow(temp3 / (*atol + ci1), 2);
		}
	}
	else
	{
		for (j = 0; j < n; j++)
		{
			temp3 = y2[j] - yjm1[j];
			ci1 = fmax(fabs(y1[j]), fabs(y2[j])) * rtol[j];
			err += pow(temp3 / (atol[j] + ci1), 2);
		}
	}

	err = sqrt(err / n);

	if (res != y2)
	{
		memcpy(res, y2, n * sizeof(double));
	}

	return (0.11399435157357929 * s + 0.02397234947266258) * err;
}


int tsrkccore(const unsigned n,
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
	double s, s2, sm1, sm1ds, s2dls;
	double w, w2m1, dtsw, d2tsw, tsm1w, dtsm1w;
	double q, onemq;
	double beta, alpha, alphaemb, gamma, gammaemb;
	bool last = false, reject = false;
	double* y2 = (double*)malloc(n * sizeof(double));
	double* yswap;
	double *res[3] = { y0, y1, y2 };
	double* f1 = (double*)malloc(n * sizeof(double));
	double* yjm1 = (double*)malloc(n * sizeof(double));
	double* yjm2 = (double*)malloc(n * sizeof(double));
	double* eigvec = (double*)malloc(n * sizeof(double));
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
		if (*h < 10 * uround)
		{
			printf("tolerances are too small\n");
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
					idid = tsrkctrho(n, x1, y1, f1, uround, f, eigvec, yjm1, yjm2, &eigmax, iwork);
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
			s2 = s * s;
			w = cosh(acosht / s);
			w2m1 = w * w - 1;
			dtsw = s * sinh(acosht) / sqrt(w2m1);
			d2tsw = (s2 * tsw - w * dtsw) / w2m1;

			sm1 = s - 1;
			sm1ds = sm1 / s;
			tsm1w = cosh(sm1ds * acosht);
			dtsm1w = sm1 * sinh(sm1ds * acosht) / sqrt(w2m1);
		}

		beta = s2 * (onemq * dtsw + sqrt(pow(onemq * dtsw, 2) + 4 * q * tsw * d2tsw)) / (2 * d2tsw);
		alpha = s2 * (1 + q) / (q * s2 * tsw + beta * dtsw);
		gamma = 1 - alpha * tsw;

		alphaemb = s2 * (1 + q) / (q * s2 * tsm1w + beta * dtsm1w);
		gammaemb = 1 - alphaemb * tsm1w;

		if (mdeg > iwork[9])
		{
			iwork[9] = mdeg;
		}

		err = tsrkcstep(n, x0, x1, *h, y0, y1, f1, atol, rtol, iwork[3], 
			mdeg, w, beta, gamma, gammaemb, acosht, f, y2, yjm1, yjm2);
		mdego = mdeg;
		iwork[5]++;
		iwork[4] += mdeg - 1;

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
				yswap = y0;
				y0 = y1;
				y1 = y2;
				y2 = yswap;
				f(&n, &x1, y1, f1);
				iwork[4]++;
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
	free(yjm1);
	free(yjm2);
	free(eigvec);

	return idid;
}


int tsrkc2(const unsigned n,
	const double x0, const double x1, const double xend,
	double* h,
	double* y0, double* y1,
	const FcnEqDiff f, const Rho rho, const SolTrait solout,
	const double* atol, const double* rtol,
	unsigned iwork[12])
{
	const double uround = 1e-16;

	if (*h > fabs(xend - x1))
	{
		printf("initial step is longer than the integration interval\n");
		return -1;
	}
	if (*h < 10 * uround)
	{
		printf("initial step-size is too small\n");
		return -1;
	}
	if (iwork[3] == 0)
	{
		if (atol[0] <= 0 || rtol[0] <= 10 * uround)
		{
			printf("tolerances are too small\n");
			return -1;
		}
	}
	else
	{
		for (unsigned i = 0; i < n; i++)
		{
			if (atol[i] <= 0 || rtol[i] <= 10 * uround)
			{
				printf("tolerances are too small\n");
				return -1;
			}
		}
	}

	iwork[4] = 0, iwork[5] = 0, iwork[6] = 0, iwork[7] = 0,
	iwork[8] = 0, iwork[9] = 0, iwork[10] = 0, iwork[11] = 0;

	return tsrkccore(n, x0, x1, xend, h, y0, y1, f, rho, solout, atol, rtol, uround, iwork);
}