/*
 * Tema 3 ASC
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
			register int idx1 = i * N + j;
			register int idx2 = j * N + i;
            T[idx1] = M[idx2];
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
			register int idx = i * N + j;
            R[idx] = M1[idx] + M2[idx];
        }
    }

    return R;
}

double* multiply_normal(int N, double *M1, double *M2) {
    double *R = (double*)calloc(N * N, sizeof(double));

    if (R == NULL) {
        return NULL;
    }

	// loop unrolling
	for (int i = 0; i < N; i++) {
		double *orig_p1 = &M1[i * N];

		for (int j = 0; j < N; j++) {
			double *p1 = orig_p1;
			double *p2 = &M2[j];

			register double r0 = 0.0;
			register double r1 = 0.0;
			register double r2 = 0.0;
			register double r3 = 0.0;

			int k = 0;

			while (k < N - 3) {
				r0 += p1[0] * p2[0 * N];
				r1 += p1[1] * p2[1 * N];
				r2 += p1[2] * p2[2 * N];
				r3 += p1[3] * p2[3 * N];
	
				p1 += 4;
				p2 += 4 * N;

				k += 4;
			}

			while (k < N) {
				r0 += *p1 * *p2;

				p1++;
				p2 += N;

				k++;
			}

			R[i * N + j] = r0 + r1 + r2 + r3;
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
		double *orig_p1 = &LM[i * N];

		for (int j = 0; j < N; j++) {
			double *p1 = orig_p1;
			double *p2 = &M[j];

			register double r0 = 0.0;
			register double r1 = 0.0;
			register double r2 = 0.0;
			register double r3 = 0.0;

			int k = 0;

			while (k <= i - 3) {
				r0 += p1[0] * p2[0 * N];
				r1 += p1[1] * p2[1 * N];
				r2 += p1[2] * p2[2 * N];
				r3 += p1[3] * p2[3 * N];

				p1 += 4;
				p2 += 4 * N;

				k += 4;
			}
			
			while (k <= i) {
				r0 += *p1 * *p2;

				p1++;
				p2 += N;

				k++;
			}

			R[i * N + j] = r0 + r1 + r2 + r3;
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
		double *orig_p1 = &M[i * N];

		for (int j = 0; j < N; j++) {
			double *p1 = orig_p1;
			double *p2 = &UM[j];

			register double r0 = 0.0;
			register double r1 = 0.0;
			register double r2 = 0.0;
			register double r3 = 0.0;

			int k = 0;

			while (k <= j - 3) {
				r0 += p1[0] * p2[0 * N];
				r1 += p1[1] * p2[1 * N];
				r2 += p1[2] * p2[2 * N];
				r3 += p1[3] * p2[3 * N];

				p1 += 4;
				p2 += 4 * N;

				k += 4;
			}

			while (k <= j) {
				r0 += *p1 * *p2;

				p1++;
				p2 += N;

				k++;
			}

			R[i * N + j] = r0 + r1 + r2 + r3;
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
