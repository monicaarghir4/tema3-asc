/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
// #include "func_neopt.h"

double* transpose(int N, double *M) {
    double *T = (double*)malloc(N * N * sizeof(double));
    if (T == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            T[i * N + j] = M[j * N + i];
        }
    }

    return T;
}

double* sum(int N, double *M1, double *M2) {
    double *R = (double*)malloc(N * N * sizeof(double));
    if (R == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * N + j] = M1[i * N + j] + M2[i * N + j];
        }
    }

    return R;
}

double* multiply_normal(int N, double *M1, double *M2) {
    double *R = (double*)malloc(N * N * sizeof(double));
    if (R == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * N + j] = 0;
            for (int k = 0; k < N; k++) {
                R[i * N + j] += M1[i * N + k] * M2[k * N + j];
            }
        }
    }

    return R;
}

double* multiply_lower_with_normal(int N, double *LM, double *M) {
    double *R = (double*)malloc(N * N * sizeof(double));
    if (R == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * N + j] = 0;
            for (int k = i; k <= N; k++) {
                R[i * N + j] += LM[i * N + k] * M[k * N + j];
            }
        }
    }

    return R;
}

/*
 * Add your unoptimized implementation here
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
