/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "func_opt.h"

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
