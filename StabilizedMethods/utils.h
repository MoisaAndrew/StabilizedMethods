#include <stdio.h>


double get_absolute_error(const unsigned n, const double* model_solution, const double* y);

double get_relative_error(const unsigned n, const double* model_solution, const double* y);

void print_error(const unsigned n, const double* model_solution, const double* y);

void to_wolfram_style(const unsigned buffer_length, char* buffer);

void print_report(const unsigned report_length, double* report[2]);