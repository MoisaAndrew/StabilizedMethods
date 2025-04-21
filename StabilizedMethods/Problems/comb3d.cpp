#include "comb3d.h"


/* right-hand-side function */
static void fcomb3d(const unsigned* n, const double* x, const double* y, double* fy)
{
	const unsigned N = (unsigned)(0.5 + pow(*n / 2, 0.3333333333333333));
	const unsigned NN = N * N;
	const double h = 1 / (0.5 + N);

	double ckji, ckjip1, ckjim1, ckjp1i, ckjm1i, ckp1ji, ckm1ji;
	double tkji, tkjip1, tkjim1, tkjp1i, tkjm1i, tkp1ji, tkm1ji;

	unsigned k2NN, j2N, i;

	for (k2NN = 0; k2NN < *n; k2NN += 2 * NN)
	{
		for (j2N = 0; j2N < 2 * NN; j2N += 2 * N)
		{
			for (i = 0; i < 2 * N; i += 2)
			{
				ckji = y[k2NN + j2N + i];
				tkji = y[k2NN + j2N + i + 1];

				if (k2NN == 0)
				{
					ckm1ji = ckji;
					ckp1ji = y[(k2NN + 2 * NN) + j2N + i];
					tkm1ji = tkji;
					tkp1ji = y[(k2NN + 2 * NN) + j2N + i + 1];
				}
				else
				{
					ckm1ji = y[(k2NN - 2 * NN) + j2N + i];
					tkm1ji = y[(k2NN - 2 * NN) + j2N + i + 1];

					if (k2NN == 2 * NN * (N - 1))
					{
						ckp1ji = 1;
						tkp1ji = 1;
					}
					else
					{
						ckp1ji = y[(k2NN + 2 * NN) + j2N + i];
						tkp1ji = y[(k2NN + 2 * NN) + j2N + i + 1];
					}
				}

				if (j2N == 0)
				{
					ckjm1i = ckji;
					ckjp1i = y[k2NN + (j2N + 2 * N) + i];
					tkjm1i = tkji;
					tkjp1i = y[k2NN + (j2N + 2 * N) + i + 1];
				}
				else
				{
					ckjm1i = y[k2NN + (j2N - 2 * N) + i];
					tkjm1i = y[k2NN + (j2N - 2 * N) + i + 1];

					if (j2N == 2 * N * (N - 1))
					{
						ckjp1i = 1;
						tkjp1i = 1;
					}
					else
					{
						ckjp1i = y[k2NN + (j2N + 2 * N) + i];
						tkjp1i = y[k2NN + (j2N + 2 * N) + i + 1];
					}
				}

				if (i == 0)
				{
					ckjim1 = ckji;
					ckjip1 = y[k2NN + j2N + i + 2];
					tkjim1 = tkji;
					tkjip1 = y[k2NN + j2N + i + 3];
				}
				else
				{
					ckjim1 = y[k2NN + j2N + i - 2];
					tkjim1 = y[k2NN + j2N + i - 1];

					if (i == 2 * (N - 1))
					{
						ckjip1 = 1;
						tkjip1 = 1;
					}
					else
					{
						ckjip1 = y[k2NN + j2N + i + 2];
						tkjip1 = y[k2NN + j2N + i + 3];
					}
				}

				fy[k2NN + j2N + i] =
					(
						ckjip1 - 2 * ckji + ckjim1 +
						ckjp1i - 2 * ckji + ckjm1i +
						ckp1ji - 2 * ckji + ckm1ji
						) / (h * h)
					-
					D * ckji * exp(-delta / tkji);
				fy[k2NN + j2N + i + 1] =
					1.1111111111111111 * (
						(
							tkjip1 - 2 * tkji + tkjim1 +
							tkjp1i - 2 * tkji + tkjm1i +
							tkp1ji - 2 * tkji + tkm1ji
							) / (h * h)
						+
						alpha * D * ckji * exp(-delta / tkji)
						);
			}
		}
	}
}


void get_comb3d(ProblemParams** params, FcnEqDiff* fcn, Rho* rho)
{
    *params = new Comb3dParams();

    *fcn = fcomb3d;
}