#include "methods_common.h"

#include "utils.h"


void solout_h(const unsigned n, const double xold, const double x, const double* y)
{
	const unsigned buffer_length = 25;
	char* buffer = (char*)malloc(buffer_length * sizeof(char));

	snprintf(buffer, buffer_length, "{%1.3f,%1.2e},", xold, x - xold);
	to_wolfram_style(buffer_length, buffer);
	printf("%s", buffer);

	free(buffer);
}