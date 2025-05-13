#include "utils.h"

#include <math.h>
#include <stdlib.h>


double euclidean_norm(const unsigned n, const double* v)
{
	double norm = 0.;
	for (unsigned j = 0; j < n; j++)
	{
		norm += v[j] * v[j];
	}
	return sqrt(norm);
}


double maximum_norm(const unsigned n, const double* v)
{
	double norm = 0., absv;
	for (unsigned j = 0; j < n; j++)
	{
		absv = fabs(v[j]);
		if (absv > norm)
		{
			norm = absv;
		}
	}
	return norm;
}


double get_absolute_error(const unsigned n, const double* model_solution, const double* y)
{
	double* errvec = (double*)malloc(n * sizeof(double));
	unsigned i;

	for (i = 0; i < n; i++)
	{
		errvec[i] = model_solution[i] - y[i];
	}

	const double error = maximum_norm(n, errvec);

	free(errvec);

	return error;
}


double get_relative_error(const unsigned n, const double* model_solution, const double* y)
{
	double* errvec = (double*)malloc(n * sizeof(double));
	unsigned i;

	for (i = 0; i < n; i++)
	{
		errvec[i] = (model_solution[i] - y[i]) / model_solution[i];
	}

	const double error = maximum_norm(n, errvec);

	free(errvec);

	return error;
}


void print_error(const unsigned n, const double* model_solution, const double* y)
{
	double* error = (double*)malloc(n * sizeof(double));
	unsigned i;

	for (i = 0; i < n; i++)
	{
		error[i] = model_solution[i] - y[i];
	}

	printf("aerr = %.2e", maximum_norm(n, error));

	for (i = 0; i < n; i++)
	{
		error[i] /= model_solution[i];
	}
	printf(", rerr = %.2e", maximum_norm(n, error));

	free(error);
}


void to_wolfram_style(const unsigned buffer_length, char* buffer)
{
	int i, j;
	for (i = buffer_length - 1; i >= 0; i--)
	{
		if (buffer[i] == 'e')
		{
			for (j = buffer_length - 2; j > i; j--)
			{
				buffer[j + 1] = buffer[j];
			}
			buffer[i + 1] = '^';
			buffer[i] = '*';
		}
	}
}


void print_report(const unsigned report_length, double* report[2])
{
	const unsigned buffer_length = 25;
	char* buffer = (char*)malloc(buffer_length * sizeof(char));
	
	if (isfinite(report[0][0]) && report[0][0] > 0)
	{
		snprintf(buffer, buffer_length, "{{%1.2e,%.3f}", report[0][0], report[1][0]);
	}
	else
	{
		snprintf(buffer, buffer_length, "{{NaN,%.3f}", report[1][0]);
	}
	to_wolfram_style(buffer_length, buffer);
	printf("%s", buffer);
	for (unsigned i = 1; i < report_length; i++)
	{
		if (isfinite(report[0][i]) && report[0][i] > 0)
		{
			snprintf(buffer, buffer_length, ",{%1.2e,%.3f}", report[0][i], report[1][i]);
		}
		else
		{
			snprintf(buffer, buffer_length, ",{NaN,%.3f}", report[1][i]);
		}
		to_wolfram_style(buffer_length, buffer);
		printf("%s", buffer);
	}
	printf("}\r\n");

	free(buffer);
}