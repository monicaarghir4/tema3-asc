/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"

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
    double *R = (double*)calloc(N * N, sizeof(double));

    if (R == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
		register double *orig_p1 = &M1[i * N];

        for (int j = 0; j < N; j++) {
			register double *p1 = orig_p1;
			register double *p2 = &M2[j];
            register double r = 0.0;

            for (int k = 0; k < N; k++) {
                r += *p1 * *p2;
				p1++;
				p2 += N;
            }

			R[i * N + j] = r;
        }
    }

    return R;
}

double* multiply_lower_with_normal(int N, double *LM, double *M) {
    double *R = (double*)calloc(N * N, sizeof(double));

    if (R == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
		register double *orig_p1 = &LM[i * N];
        for (int j = 0; j < N; j++) {
			register double *p1 = orig_p1;
			register double *p2 = &M[j];
            register double r = 0.0;
            for (int k = 0; k <= i; k++) {
				r += *p1 * *p2;
				p1++;
				p2 += N;
            }
			R[i * N + j] = r;
        }
    }

    return R;
}

double* multiply_normal_with_upper(int N, double *M, double *UM) {
    double *R = (double*)calloc(N * N, sizeof(double));

	if (R == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
		register double *orig_p1 = &M[i * N];
        for (int j = 0; j < N; j++) {
			register double *p1 = orig_p1;
			register double *p2 = &UM[j];
            register double r = 0.0;
            for (int k = 0; k <= j; k++) {
				r += *p1 * *p2;
				p1++;
				p2 += N;
            }
			R[i * N + j] = r;
        }
    }

  	return R;
}

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	double *C;

	double *A_T = transpose(N, A);
	double *B_T = transpose(N, B);

	double *first_calc = multiply_lower_with_normal(N, A_T, B);

	double *second_calc = multiply_normal_with_upper(N, B, A);

	double *third_calc = sum(N, first_calc, second_calc);

	C = multiply_normal(N, third_calc, B_T);

	free(A_T);
	free(B_T);
	free(first_calc);
	free(second_calc);
	free(third_calc);

	return C;
}
