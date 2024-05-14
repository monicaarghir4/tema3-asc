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

    for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
			printf("%lf ", R[i * N + j]);
		printf("\n");
	}

	printf("\n");

    return R;
}

double* multiply_normal_with_upper(int N, double *M, double *UM) {
    double *R = (double*)malloc(N * N * sizeof(double));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * N + j] = 0;
            for (int k = i; k <= N; k++) {
                R[i * N + j] += M[i * N + k] * UM[k * N + j];
            }
        }
    }

  	return R;
}

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	double* C = (double*)malloc(N * N * sizeof(double));

	double *A_T = (double*)malloc(N * N * sizeof(double));
    if (A_T == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A_T[i * N + j] = A[j * N + i];
        }
    }

	double *B_T = (double*)malloc(N * N * sizeof(double));
    if (B_T == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B_T[i * N + j] = B[j * N + i];
        }
    }

	// for (int i = 0; i < N; i++) {
	// 	for (int j = 0; j < N; j++)
	// 		printf("%lf ", A_T[i * N + j]);
	// 	printf("\n");
	// }

	// printf("\n");

	// for (int i = 0; i < N; i++) {
	// 	for (int j = 0; j < N; j++)
	// 		printf("%lf ", B_T[i * N + j]);
	// 	printf("\n");
	// }

	// printf("\n");

	double* first_calc = (double*)malloc(N * N * sizeof(double));

	for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            first_calc[i * N + j] = 0;
            for (int k = 0; k <= i; k++) {
                first_calc[i * N + j] += A_T[i * N + k] * B[k * N + j];
            }
        }
    }

	// for (int i = 0; i < N; i++) {
	// 	for (int j = 0; j < N; j++)
	// 		printf("%lf ", first_calc[i * N + j]);
	// 	printf("\n");
	// }

	// printf("\n");

	double* second_calc = (double*)malloc(N * N * sizeof(double));

	for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            second_calc[i * N + j] = 0;
            for (int k = 0; k <= j; k++) {
                second_calc[i * N + j] += B[i * N + k] * A[k * N + j];
            }
        }
    }

	// for (int i = 0; i < N; i++) {
	// 	for (int j = 0; j < N; j++)
	// 		printf("%lf ", second_calc[i * N + j]);
	// 	printf("\n");
	// }

	// printf("\n");

	double* third_calc = (double*)malloc(N * N * sizeof(double));

	for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            third_calc[i * N + j] = first_calc[i * N + j] + second_calc[i * N + j];
        }
    }

	// for (int i = 0; i < N; i++) {
	// 	for (int j = 0; j < N; j++)
	// 		printf("%lf ", third_calc[i * N + j]);
	// 	printf("\n");
	// }

	// printf("\n");

	for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i * N + j] = 0;
            for (int k = 0; k < N; k++) {
                C[i * N + j] += third_calc[i * N + k] * B_T[k * N + j];
            }
        }
    }

	// for (int i = 0; i < N; i++) {
	// 	for (int j = 0; j < N; j++)
	// 		printf("%lf ", C[i * N + j]);
	// 	printf("\n");
	// }

	// printf("\n");

	// double* second_calc = multiply_normal_with_upper(N, B, A);

	// double* third_calc = sum(N, first_calc, second_calc);

	// C = multiply_normal(N, third_calc, B_T);
	return C;
}


// int main() {
// 	int N = 3;
// 	double A[] = {1, 2, 3, 0, 4, 5, 0, 0, 6};
// 	double B[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

// 	double *C = my_solver(N, A, B);

// 	for (int i = 0; i < N; i++) {
// 		for (int j = 0; j < N; j++)
// 			printf("%lf ", A[i * N + j]);
// 		printf("\n");
// 	}

// 	printf("\n");

// 	for (int i = 0; i < N; i++) {
// 		for (int j = 0; j < N; j++)
// 			printf("%lf ", B[i * N + j]);
// 		printf("\n");
// 	}

// 	printf("\n");

// 	for (int i = 0; i < N; i++) {
// 		for (int j = 0; j < N; j++)
// 			printf("%lf ", C[i * N + j]);
// 		printf("\n");
// 	}

// 	printf("\n");

// 	return 0;
// }