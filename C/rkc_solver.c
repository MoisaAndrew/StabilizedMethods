//	Version of February 2025


#include "methods_common.h"


static double rkc_rho(const unsigned n, const double x, const FcnEqDiff f,
	double* yn, double* fn, double* v, double* fv, double* eigvec,
	const double hmax, const double uround, unsigned iwork[10])
{
	const double safe = 1.2;
	const unsigned itmax = 50;
	const double sqrtu = sqrt(uround), small = 1. / hmax;
	double ynrm = 0., vnrm = 0., sprad, sigma, sigmal, dynrm, dfnrm;
	unsigned i, iter, index;

	if (iwork[5] == 0)
	{
		memcpy(v, fn, n * sizeof(double));
	}
	else
	{
		memcpy(v, eigvec, n * sizeof(double));
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
		sprad = safe * sigma;

		if (iter >= 2 && fabs(sigma - sigmal) <= 0.01 * fmax(sigma, small))
		{
			for (i = 0; i < n; i++)
			{
				eigvec[i] = v[i] - yn[i];
			}
			return sprad;
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

	return -1;
}


static void step_rkc(const unsigned n, const double x, const FcnEqDiff f,
	const double* yn, const double* fn, const double h, const unsigned m,
	double* y, double* yjm1, double* yjm2)
{
	unsigned i;
	double* swap, * res = y;

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
	double d2zjm1 = 0., d2zjm2 = 0.;

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
			swap = yjm2;
			yjm2 = yjm1;
			yjm1 = y;
			y = swap;

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

	if (res != y)
	{
		memcpy(res, y, n * sizeof(double));
	}
}


static void step_tsrkc2(const unsigned n,
	const double x1, const double h,
	const double* y0, const double* y1, const double* f1,
	const unsigned m,
	const double w, const double beta, const double gamma, const double acosht,
	const FcnEqDiff f,
	double* y2, double* yjm1, double* yjm2)
{
	double* yswap;
	double* res = y2;
	const double hbeta = h * beta;
	const double acoshtds = acosht / m;
	double coshims = 1, coshi = w, tim1wdtiw;
	double temp1 = hbeta / w, temp2, temp3;
	double ci1 = x1 + temp1, ci2 = x1 + temp1, ci3 = x1;
	unsigned i, j;

	memcpy(yjm2, y1, n * sizeof(double));
	for (j = 0; j < n; j++)
	{
		yjm1[j] = yjm2[j] + temp1 * f1[j];
	}

	for (i = 2; i <= m; i++)
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

	temp3 = 1 - gamma;
	for (j = 0; j < n; j++)
	{
		res[j] = temp3 * y0[j] + gamma * yjm1[j];
	}
}


static void step_tsrkc3(const unsigned n,
	const double x1, const double h,
	const double* y0, const double* y1, const double* f1,
	const unsigned m,
	const double w, const double acosht,
	const double beta, const double deltadchi,
	const double gamma, const double chi,
	const FcnEqDiff f,
	double* y2, double* yjm1, double* yjm2)
{
	double* yswap;
	double* res = y2;
	const double hbeta = h * beta, acoshtds = acosht / m;
	const double chim1 = chi - 1;
	double coshims = 1, coshi = w, tim1wdtiw;
	double temp1 = chi * hbeta / w, temp2, temp3;
	double ci1 = x1 + temp1, ci2 = x1 + temp1, ci3 = x1;
	unsigned i, j;

	memcpy(yjm2, y1, n * sizeof(double));
	for (j = 0; j < n; j++)
	{
		yjm1[j] = yjm2[j] + temp1 * f1[j];
	}

	for (i = 2; i <= m; i++)
	{
		coshims = coshi;
		coshi = cosh(i * acoshtds);
		tim1wdtiw = coshims / coshi;
		temp1 = 2 * hbeta * tim1wdtiw;
		temp2 = 2 * w * tim1wdtiw;
		temp3 = 1 - temp2;
		f(&n, &ci1, yjm1, y2);
		ci1 = chi * temp1 + temp2 * ci2 + temp3 * ci3;
		for (j = 0; j < n; j++)
		{
			y2[j] = temp1 * (y2[j] + chim1 * f1[j]) + temp2 * yjm1[j] + temp3 * yjm2[j];
		}

		yswap = yjm2;
		yjm2 = yjm1;
		yjm1 = y2;
		y2 = yswap;

		ci3 = ci2;
		ci2 = ci1;
	}

	temp2 = 1 - gamma - deltadchi;
	for (j = 0; j < n; j++)
	{
		res[j] = temp2 * y1[j] + gamma * y0[j] + deltadchi * yjm1[j];
	}
}


static double calc_chi(const double tsw, const unsigned mdeg, const double q,
	const double w, const double acosht,
	const double beta, const double delta, const double gamma)
{
	const double acoshtds = acosht / mdeg, sinhacoshtds = sinh(acoshtds);
	double coshims = 1, coshi = w, tim1wdtiw;

	double sumbc2 = 0, ci3 = 0, ci2 = beta / w, ci1, m;

	double b = sinh((mdeg - 1) * acoshtds) * w / sinhacoshtds;

	sumbc2 += b * ci2 * ci2;

	for (unsigned i = 2; i < mdeg; i++)
	{
		coshims = coshi;
		coshi = cosh(i * acoshtds);
		tim1wdtiw = coshims / coshi;
		m = 2 * w * tim1wdtiw;

		ci1 = m * ci2 + (1 - m) * ci3 + 2 * beta * tim1wdtiw;
		b = sinh((mdeg - i) * acoshtds) * coshi / sinhacoshtds;

		sumbc2 += b * ci1 * ci1;

		ci3 = ci2;
		ci2 = ci1;
	}

	sumbc2 *= 6 * beta * delta / tsw;

	return (1 + gamma * pow(q, 3)) / sumbc2;
}


static int rkc_core(const unsigned n, double x, const double xend, double* y,
	const FcnEqDiff f, const Rho rho, const SolTrait solout,
	const double* atol, const double rtol, const double uround,
	double* work, unsigned iwork[10], const int method)
{
	bool newspc = true, jacatt = false, last;
	unsigned mmax = fmax(sqrt(rtol / (10. * uround)), 2);
	unsigned nstsig = 0;
	const double tdir = xend - x > 0 ? 1. : -1.;
	const double hmax = tdir * (xend - x);
	const double facmax = method == 0 ? 2 : 10;
	const double errpow = method == 1 ? 0.5 : 0.3333333333333333;
	double hmin = 10. * uround * fmax(fabs(x), hmax);
	double absh, h, abshold, est, sprad, wt, at, err, errold, fac, temp, ci;
	double q, onemq, onepq, onepq2, s, s2, w, w2, w2m1, dtsw, d2tsw, d3tsw;
	double alpha, beta, gamma, delta, chi;
	unsigned i, m, mold = 0;

	double* yn = &work[0];
	double* fn = &work[n];
	double* vtemp1 = &work[2 * n];
	double* vtemp2 = &work[3 * n];
	
	double* yold = NULL;
	double* eigvec = NULL;
	if (method != 2) // TSRKC3 or TSRKC2
	{
		yold = &work[4 * n];
		if (iwork[0] == 0)
		{
			eigvec = &work[5 * n];
		}
	} 
	else if (iwork[0] == 0) // RKC
	{
		eigvec = &work[4 * n];
	}

	double* swap, * res = y;

	const double tsw = method == 0 ? 1.25 : 1.1;
	const double acosht = acosh(tsw);

	memcpy(yn, y, n * sizeof(double));
	f(&n, &x, yn, fn);
	iwork[4]++;

	if (iwork[0] == 1)
	{
		rho(&n, &x, yn, &sprad);
	}
	else
	{
		sprad = rkc_rho(n, x, f, yn, fn, vtemp1, vtemp2, eigvec, hmax, uround, iwork);
		if (sprad < 0)
		{
			return 6;
		}
	}
	jacatt = true;

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
	step_rkc(n, x, f, yn, fn, h, m, y, vtemp1, vtemp2);
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

		if (method == 0)
		{
			est = 1.2 * (yn[i] - y[i]) + 0.6 * h * (fn[i] + vtemp1[i]);
		}
		else if (method == 1)
		{
			est = 0.3333333333333333 * (yn[i] - y[i] + h * vtemp1[i]);
		}
		else
		{
			est = 0.8 * (yn[i] - y[i]) + 0.4 * h * (fn[i] + vtemp1[i]);
		}

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
		absh *= 0.8 * pow(err, -0.3333333333333333);
		if (absh < hmin)
		{
			return 4;
		}
		else
		{
			newspc = !jacatt;
		}
	}

	if (newspc)
	{
		if (iwork[2] == 1)
		{
			solout(n, x, x + h, y);
		}

		iwork[6]++;
		x += h;

		if (last)
		{
			if (res != y)
			{
				memcpy(res, y, n * sizeof(double));
			}
			return 1;
		}

		jacatt = iwork[1] == 1;
		nstsig = (nstsig + 1) % 25;
		if (iwork[0] == 1 || nstsig == 0)
		{
			newspc = !jacatt;
		}
		else
		{
			newspc = false;
		}

		swap = vtemp2;
		vtemp2 = fn;
		fn = vtemp1;
		if (method == 2) // RKC
		{
			vtemp1 = yn;
		}
		else // TSRKC3 or TSRKC2
		{
			vtemp1 = yold;
			yold = yn;
		}
		yn = y;
		y = swap;

		abshold = absh;
		fac = pow(err, -errpow);
		absh *= fmax(0.1, fmin(0.8 * fac, facmax));
		absh = fmax(hmin, fmin(hmax, absh));
		errold = err;
	}

	while (true)
	{
		if (newspc)
		{
			if (iwork[0] == 1)
			{
				rho(&n, &x, yn, &sprad);
			}
			else
			{
				sprad = rkc_rho(n, x, f, yn, fn, vtemp1, vtemp2, eigvec, hmax, uround, iwork);
				if (sprad < 0)
				{
					return 6;
				}
			}
			jacatt = true;
		}

		last = false;
		if (1.1 * absh >= fabs(xend - x))
		{
			absh = fabs(xend - x);
			last = true;
		}

		q = abshold / absh;
		onemq = 1 - q;

		if (method == 2 || iwork[6] == 0) // RKC
		{
			m = 1 + (unsigned)sqrt(1.54 * absh * sprad + 1.);
			if (m > mmax)
			{
				m = mmax;
				absh = (m * m - 1.) / (1.54 * sprad);
				q = abshold / absh;
				last = false;
			}

			h = tdir * absh;
			step_rkc(n, x, f, yn, fn, h, m, y, vtemp1, vtemp2);
		}
		else if (method == 0) // TSRKC3
		{
			onepq = 1 + q;
			onepq2 = onepq * onepq;

			m = 1 + (unsigned)sqrt(absh * sprad * 1.267029788142009 * (onemq + sqrt(1 + q * (0.44256220745562963 + q))));
			if (m < 3)
			{
				m = 3;
			}
			else if (m > mmax)
			{
				m = mmax;
				absh = 0.789247426823654 * m * (m + 2.534059576284018 * abshold * sprad) 
					/ (sprad * (2 * m + 3.094799076236184 * abshold * sprad));
				q = abshold / absh;
				onemq = 1 - q;
				onepq = 1 + q;
				onepq2 = onepq * onepq;
				last = false;
			}

			if (m != mold)
			{
				s = (double)m;
				s2 = s * s;
				w = cosh(acosht / s);
				w2 = w * w;
				w2m1 = w2 - 1;
				dtsw = s * sinh(acosht) / sqrt(w2m1);
				d2tsw = (s2 * tsw - w * dtsw) / w2m1;
				d3tsw = ((1 + 2 * w2 + s2 * w2m1) * dtsw - 3 * s2 * w * tsw) / (w2m1 * w2m1);
			}

			beta = (onemq * d2tsw + sqrt(pow(onemq * d2tsw, 2) + 4 * q * dtsw * d3tsw)) / (2 * d3tsw);
			alpha = onepq / (q * dtsw + beta * d2tsw);
			gamma = (alpha * dtsw - 1) / q;
			delta = alpha * tsw / beta;
			chi = calc_chi(tsw, m, q, w, acosht, beta, delta, gamma);

			h = tdir * absh;
			step_tsrkc3(n, x, h, yold, yn, fn, m, w, acosht, beta, delta / chi, gamma, chi, f, y, vtemp1, vtemp2);
		}
		else if (method == 1) // TSRKC2
		{
			m = 1 + (unsigned)sqrt(absh * sprad * 0.759782816506459 * (onemq + sqrt(1 + q * (q - 0.598626091572911))));

			if (m < 2)
			{
				m = 2;
			}
			else if (m > mmax)
			{
				m = mmax;
				absh = 1.316165591370016 * m * (m + 1.519565633012918 * abshold * sprad)
					/ (sprad * (2 * m + 1.9743914509024376 * abshold * sprad));
				q = abshold / absh;
				onemq = 1 - q;
				last = false;
			}

			if (m != mold)
			{
				s = (double)m;
				w = cosh(acosht / s);
				w2m1 = w * w - 1;
				dtsw = s * sinh(acosht) / sqrt(w2m1);
				d2tsw = (s * s * tsw - w * dtsw) / w2m1;
			}

			beta = (onemq * dtsw + sqrt(pow(onemq * dtsw, 2) + 4 * q * tsw * d2tsw)) / (2 * d2tsw);
			gamma = (1 + q) * tsw / (q * tsw + beta * dtsw);

			h = tdir * absh;
			step_tsrkc2(n, x, h, yold, yn, fn, m, w, beta, gamma, acosht, f, y, vtemp1, vtemp2);
		}

		hmin = 10. * uround * fmax(fabs(x), fabs(x + h));
		mold = m;
		iwork[9] = m > iwork[9] ? m : iwork[9];

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

			if (method == 0)
			{
				if (iwork[6] == 0)
				{
					est = 1.2 * (yn[i] - y[i]) + 0.6 * h * (fn[i] + vtemp1[i]);
				}
				else
				{
					est = 0.6 * (yn[i] / q - yold[i] / (q * onepq2) - y[i] * (2 + q) / onepq2 + h * vtemp1[i] / onepq);
				}
			}
			else if (method == 1)
			{
				est = 0.3333333333333333 * (yn[i] - y[i] + h * vtemp1[i]);
			}
			else if (method == 2)
			{
				est = 0.8 * (yn[i] - y[i]) + 0.4 * h * (fn[i] + vtemp1[i]);
			}

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
				
			absh *= 0.8 * pow(err, -errpow);
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

		if (iwork[2] == 1)
		{
			solout(n, x, x + h, y);
		}

		iwork[6]++;
		x += h;

		if (last)
		{
			if (res != y)
			{
				memcpy(res, y, n * sizeof(double));
			}
			return 1;
		}

		jacatt = iwork[1] == 1;
		nstsig = (nstsig + 1) % 25;
		if (iwork[0] == 1 || nstsig == 0)
		{
			newspc = !jacatt;
		}
		else
		{
			newspc = false;
		}

		swap = vtemp2;
		vtemp2 = fn;
		fn = vtemp1;
		if (method == 2)
		{
			vtemp1 = yn;
		}
		else
		{
			vtemp1 = yold;
			yold = yn;
		}
		yn = y;
		y = swap;

		abshold = absh;
		fac = pow(err, -errpow);
		if (iwork[6] > 1)
		{
			temp = pow(errold, errpow) * fac * fac / q;
			fac = fmin(fac, temp);
		}
		absh *= fmax(0.1, fmin(0.8 * fac, facmax));
		absh = fmax(hmin, fmin(hmax, absh));
		errold = err;
	}
}


int rkc_solver(const unsigned n, const double x, const double xend, double* y,
	const FcnEqDiff f, const Rho rho, const SolTrait solout,
	const double* atol, const double rtol,
	unsigned iwork[10], const int method)
{
	const double uround = 1e-16;

	if (n <= 0 || method < 0 || method > 2)
	{
		return 5;
	}

	if (atol[0] <= 0 || rtol < 10. * uround)
	{
		return 5;
	}

	if (iwork[3])
	{
		for (unsigned i = 1; i < n; i++)
		{
			if (atol[i] <= 0)
			{
				return 5;
			}
		}
	}

	iwork[4] = 0, iwork[5] = 0, iwork[6] = 0,
	iwork[7] = 0, iwork[8] = 0, iwork[9] = 0;

	unsigned long long work_length = 4;
	if (iwork[0] == 0)
	{
		work_length++;
	}
	if (method != 2)
	{
		work_length++;
	}
	double* work = (double*)malloc(work_length * n * sizeof(double));

	int idid = rkc_core(n, x, xend, y, f, rho, solout, atol, rtol, uround, work, iwork, method);

	free(work);

	return idid;
}