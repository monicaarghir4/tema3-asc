/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
#include <cblas.h>

/* 
 * Add your BLAS implementation here
 */

double* my_solver(int N, double *A, double *B) {
	double* first_calc = (double *)malloc(N * N * sizeof(double));

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			first_calc[i * N + j] = B[i * N + j];
		}
	}

	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
                N, N, 1.0, A, N, first_calc, N);

	double* second_calc = (double *)malloc(N * N * sizeof(double));

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			second_calc[i * N + j] = B[i * N + j];
		}
	}

	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
				N, N, 1.0, A, N, second_calc, N);

	cblas_daxpy(N * N, 1.0, first_calc, 1, second_calc, 1);

	double* C = (double *)malloc(N * N * sizeof(double));


	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0,
					second_calc, N, B, N, 0.0, C, N);

	free(first_calc);
	free(second_calc);
	return C;
}
