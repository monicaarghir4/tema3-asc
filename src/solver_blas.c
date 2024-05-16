/*
 * Tema 3 ASC
 * 2024 Spring
 */
#include "utils.h"
#include <cblas.h>

/* 
	Implementare folosind biblioteca blas pentru C = (AT * B + B * A) * BT
 */

double* my_solver(int N, double *A, double *B) {
	// alocarea de memorie pentru matricea rezultata in primul calcul
	double* first_calc = (double *)malloc(N * N * sizeof(double));

	// daca alocarea esueaza, intoarcem NULL
    if (first_calc == NULL) {
        return NULL;
    }

	// copierea matricei B in first_calc
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			first_calc[i * N + j] = B[i * N + j];
		}
	}

	// calcul AT * B
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
                N, N, 1.0, A, N, first_calc, N);

	// alocarea de memorie pentru matricea rezultata in al doilea calcul
	double* second_calc = (double *)malloc(N * N * sizeof(double));

	// daca alocarea esueaza, intoarcem NULL
    if (second_calc == NULL) {
        return NULL;
    }

	// copierea matricei B in second_calc
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			second_calc[i * N + j] = B[i * N + j];
		}
	}

	// calcul B * A
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
				N, N, 1.0, A, N, second_calc, N);

	// calcul C = AT * B + B * A
	cblas_daxpy(N * N, 1.0, first_calc, 1, second_calc, 1);

	// alocarea de memorie pentru matricea rezultata
	double* C = (double *)malloc(N * N * sizeof(double));

	// daca alocarea esueaza, intoarcem NULL
    if (C == NULL) {
        return NULL;
    }

	// calcul C = (AT * B + B * A) * BT
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0,
					second_calc, N, B, N, 0.0, C, N);

	// eliberam memoria
	free(first_calc);
	free(second_calc);
	
	return C;
}
