#include "rkc.h"


int rkcrho(const unsigned n, const double x, const FcnEqDiff f,
	double* yn, double* fn, double* v, double* fv, double* work,
	const double hmax, const double uround, double* sprad,
	unsigned iwork[12])
{
	const double safe = 1.2;
	const double sqrtu = sqrt(uround), small = 1. / hmax;
	double ynrm = 0., vnrm = 0., sigma, sigmal, dynrm, dfnrm;
	unsigned i, iter, itmax = 50, index;

	if (iwork[5] == 0)
	{
		memcpy(v, fn, n * sizeof(double));
	}
	else
	{
		memcpy(v, work, n * sizeof(double));
	}

	for (i = 0; i < n; i++)
	{
		ynrm += yn[i] * yn[i];
		vnrm += v[i] * v[i];
	}
	ynrm = sqrt(ynrm);
	vnrm = sqrt(vnrm);
	
	if (ynrm != 0.)
	{
		dynrm = ynrm * sqrtu;
		if (vnrm != 0.)
		{
			for (i = 0; i < n; i++)
			{
				v[i] = yn[i] + v[i] * (dynrm / vnrm);
			}
		}
		else
		{
			for (i = 0; i < n; i++)
			{
				v[i] = yn[i] + yn[i] * sqrtu;
			}
		}
	}
	else
	{
		dynrm = uround;
		if (vnrm != 0.)
		{
			dynrm = uround;
			for (i = 0; i < n; i++)
			{
				v[i] *= (dynrm / vnrm);
			}
		}
		else
		{
			for (i = 0; i < n; i++)
			{
				v[i] = dynrm;
			}
		}
	}

	sigma = 0.;
	for (iter = 1; iter <= itmax; iter++)
	{
		f(&n, &x, v, fv);
		iwork[8]++;			
		dfnrm = 0.;
		for (i = 0; i < n; i++)
		{
			dfnrm += pow(fv[i] - fn[i], 2);
		}
		dfnrm = sqrt(dfnrm);
		sigmal = sigma;
		sigma = dfnrm / dynrm;
		*sprad = safe * sigma;

		if (iter >= 2 && fabs(sigma - sigmal) <= 0.01 * fmax(sigma, small))
		{
			for (i = 0; i < n; i++)
			{
				work[i] = v[i] - yn[i];
			}
			return 1;
		}

		if (dfnrm != 0.)
		{
			for (i = 0; i < n; i++)
			{
				v[i] = yn[i] + (fv[i] - fn[i]) * (dynrm / dfnrm);
			}
		}
		else
		{
			index = 1 + iter % n;
			v[index] = yn[index] - (v[index] - yn[index]);
		}
	}

	return 6;
}


void step(const unsigned n, const double x, const FcnEqDiff f,
	double* yn, double* fn, const double h, const unsigned m, 
	double* y, double* yjm1, double* yjm2)
{
	unsigned i;

	double w0 = 1. + 2. / (13. * m * m);
	double temp1 = w0 * w0 - 1., temp2 = sqrt(temp1);
	double arg = m * log(w0 + temp2);
	double w1 = sinh(arg) * temp1 / (cosh(arg) * m * temp2 - w0 * sinh(arg));
	double bjm1 = 1. / (4. * w0 * w0), bjm2 = bjm1;

	memcpy(yjm2, yn, n * sizeof(double));
	double mus = w1 * bjm1;
	for (i = 0; i < n; i++)
	{
		yjm1[i] = yn[i] + h * mus * fn[i];
	}
	double thjm2 = 0., thjm1 = mus;
	double zjm1 = w0, zjm2 = 1.;
	double dzjm1 = 1., dzjm2 = 0.;
	double d2zjm1 = 0.,	d2zjm2 = 0.;

	double zj, dzj, d2zj, bj, ajm1, mu, nu, thj, cj;
	for (unsigned j = 2; j <= m; j++)
	{
		zj = 2. * w0 * zjm1 - zjm2;
		dzj = 2. * w0 * dzjm1 - dzjm2 + 2. * zjm1;
		d2zj = 2. * w0 * d2zjm1 - d2zjm2 + 4. * dzjm1;
		bj = d2zj / (dzj * dzj);
		ajm1 = 1. - zjm1 * bjm1;
		mu = 2. * w0 * bj / bjm1;
		nu = -bj / bjm2;
		mus = mu * w1 / w0;

		cj = x + h * thjm1;
		f(&n, &cj, yjm1, y);
		for (i = 0; i < n; i++)
		{
			y[i] = mu * yjm1[i] + nu * yjm2[i] + (1. - mu - nu) * yn[i] + h * mus * (y[i] - ajm1 * fn[i]);
		}
		thj = mu * thjm1 + nu * thjm2 + mus * (1. - ajm1);

		if (j < m)
		{
			memcpy(yjm2, yjm1, n * sizeof(double));
			memcpy(yjm1, y, n * sizeof(double));
			thjm2 = thjm1;
			thjm1 = thj;
			bjm2 = bjm1;
			bjm1 = bj;
			zjm2 = zjm1;
			zjm1 = zj;
			dzjm2 = dzjm1;
			dzjm1 = dzj;
			d2zjm2 = d2zjm1;
			d2zjm1 = d2zj;
		}
	}
}


int rkclow(const unsigned n, double x, const double xend, double* y,
	const FcnEqDiff f, const Rho rho, const SolTrait solout,
	unsigned iwork[10], const double rtol, const double* atol, const double uround)
{
	bool newspc = true, jacatt = false, last;
	unsigned mmax = fmax(sqrt(rtol / (10. * uround)), 2);
	unsigned nstsig = 0;
	const double tdir = xend - x > 0 ? 1. : -1.;
	const double hmax = tdir * xend - x;
	double hmin = 10. * uround * fmax(fabs(x), hmax);
	double absh, h, hold, est, sprad, wt, at, err, errold, fac, temp1, temp2, ci;
	unsigned i, m;

	double* yn = (double*)malloc(n * sizeof(double));
	double* fn = (double*)malloc(n * sizeof(double));
	double* vtemp1 = (double*)malloc(n * sizeof(double));
	double* vtemp2 = (double*)malloc(n * sizeof(double));
	double* work = (double*)malloc(n * sizeof(double));
	int idid = 1;

	memcpy(yn, y, n * sizeof(double));
	f(&n, &x, yn, fn);
	iwork[4]++;

	while (true)
	{
		if (newspc)
		{
			if (iwork[1] == 1)
			{
				rho(&n, &x, yn, &sprad);
			}
			else
			{
				idid = rkcrho(n, x, f, yn, fn, vtemp1, vtemp2, work, hmax, uround, &sprad, iwork);
				if (idid == 6)
				{
					return idid;
				}
			}
			jacatt = true;
		}

		if (iwork[5] == 0)
		{
			absh = hmax;
			if (sprad * absh > 1.)
			{
				absh = 1. / sprad;
			}
			absh = fmax(absh, hmin);

			for (i = 0; i < n; i++)
			{
				vtemp1[i] = yn[i] + absh * fn[i];
			}
			ci = x + absh;
			f(&n, &ci, vtemp1, vtemp2);
			iwork[4]++;
			est = 0.;
			at = atol[0];
			for (i = 0; i < n; i++)
			{
				if (iwork[3])
				{
					at = atol[i];
				}
				wt = at + rtol * fabs(yn[i]);
				if (wt == 0.)
				{
					return 3;
				}
				est += pow((vtemp2[i] - fn[i]) / wt, 2);
			}
			est = absh * sqrt(est / n);
			if (0.1 * absh < hmax * sqrt(est))
			{
				absh = fmax(0.1 * absh / sqrt(est), hmin);
			}
			else
			{
				absh = hmax;
			}
		}

		last = false;
		if (1.1 * absh >= fabs(xend - x))
		{
			absh = fabs(xend - x);
			last = true;
		}
		m = 1 + (unsigned)sqrt(1.54 * absh * sprad + 1.);

		if (m > mmax)
		{
			m = mmax;
			absh = (m * m - 1.) / (1.54 * sprad);
			last = false;
		}
		iwork[9] = m > iwork[9] ? m : iwork[9];

		h = tdir * absh;
		hmin = 10. * uround * fmax(fabs(x), fabs(x + h));
		step(n, x, f, yn, fn, h, m, y, vtemp1, vtemp2);
		ci = x + h;
		f(&n, &ci, y, vtemp1);
		iwork[4] += m;
		iwork[5]++;

		err = 0.;
		at = atol[0];
		for (i = 0; i < n; i++)
		{
			if (iwork[3])
			{
				at = atol[i];
			}
			wt = at + rtol * fmax(fabs(y[i]), fabs(yn[i]));
			if (wt == 0)
			{
				return 3;
			}
			est = 0.8 * (yn[i] - y[i]) + 0.4 * h * (fn[i] + vtemp1[i]);
			err = err + pow(est / wt, 2);
		}
		err = sqrt(err / n);

		if (!isfinited(err))
		{
			return 4;
		}
		if (err > 1.)
		{
			iwork[7]++;
			absh = 0.8 * absh / pow(err, 0.3333333333333333);
			if (absh < hmin)
			{
				return 4;
			}
			else
			{
				newspc = !jacatt;
				continue;
			}
		}

		if (iwork[0] == 0)
		{
			solout(n, x, x + h, y);
		}

		iwork[6]++;
		x += h;
		jacatt = iwork[2] == 1;
		nstsig = (nstsig + 1) % 25;
		if (iwork[1] == 1 || nstsig == 0)
		{
			newspc = !jacatt;
		}
		else
		{
			newspc = false;
		}

		memcpy(vtemp2, fn, n * sizeof(double));
		memcpy(fn, vtemp1, n * sizeof(double));
		memcpy(vtemp1, yn, n * sizeof(double));
		memcpy(yn, y, n * sizeof(double));
		fac = 10.;
		if (iwork[6] == 1)
		{
			temp2 = pow(err, 0.3333333333333333);
			if (0.8 < fac * temp2)
			{
				fac = 0.8 / temp2;
			}
		}
		else
		{
			temp1 = 0.8 * absh * pow(errold, 0.3333333333333333);
			temp2 = fabs(hold) * pow(err, 0.6666666666666667);
			if (temp1 < fac * temp2)
			{
				fac = temp1 / temp2;
			}
		}
		absh = fmax(0.1, fac) * absh;
		absh = fmax(hmin, fmin(hmax, absh));
		errold = err;
		hold = h;
		h = tdir * absh;
		if (last)
		{
			return 1;
		}
	}

	free(yn);
	free(fn);
	free(vtemp1);
	free(vtemp2);
	free(work);
}


int rkcc(const unsigned n, const double x, const double xend, double* y,
	const FcnEqDiff f, const Rho rho, const SolTrait solout,
	const double* atol, const double rtol,
	unsigned iwork[10])
{
	const double uround = 2.22e-16;

	if (n < 0)
	{
		return 5;
	}
	if (rtol > 0.1 || rtol < 10. * uround)
	{
		return 5;
	}
	if (atol[0] < 0)
	{
		return 5;
	}
	if (iwork[3])
	{
		for (unsigned i = 1; i < n; i++)
		{
			if (atol[i] < 0)
			{
				return 5;
			}
		}
	}

	iwork[4] = 0, iwork[5] = 0, iwork[6] = 0,
	iwork[7] = 0, iwork[8] = 0, iwork[9] = 0;

	return rkclow(n, x, xend, y, f, rho, solout, iwork, rtol, atol, uround);
}