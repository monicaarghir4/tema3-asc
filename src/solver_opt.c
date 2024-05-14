/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
#include "func_opt.c"

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	double* C;

	double* A_T = transpose(N, A);
	double* B_T = transpose(N, B);

	double* first_calc = multiply_lower_with_normal(N, A_T, B);
	double* first_calc_T = transpose(N, first_calc);

	double* second_calc = sum(N, first_calc, first_calc_T);

	C = multiply_normal(N, second_calc, B_T);
	return C;
}
