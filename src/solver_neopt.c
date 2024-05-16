/*
 * Tema 3 ASC
 * 2024 Spring
 */
#include "utils.h"

/*
    Functie care calculeaza transpusa unei matrici
*/
double* transpose(int N, double *M) {
    // alocarea memoriei pentru matrice
    double *T = (double*)malloc(N * N * sizeof(double));

    // daca alocarea esueaza, intoarcem NULL
    if (T == NULL) {
        return NULL;
    }

    // setam elementele transpusei
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            T[i * N + j] = M[j * N + i];
        }
    }

    return T;
}

/*
    Functie care calculeaza suma a doua matrici
*/
double* sum(int N, double *M1, double *M2) {
    // alocarea memoriei pentru matricea rezultat
    double *R = (double*)malloc(N * N * sizeof(double));

    // daca alocarea esueaza, intoarcem NULL
    if (R == NULL) {
        return NULL;
    }

    // punem in matricea rezultat, suma elementelor de pe aceleasi
    // pozitii din cele doua matrici
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * N + j] = M1[i * N + j] + M2[i * N + j];
        }
    }

    return R;
}

/*
    Functie care calculeaza produsul a doua matrici
*/
double* multiply_normal(int N, double *M1, double *M2) {
    // alocarea memoriei pentru matricea rezultat
    double *R = (double*)malloc(N * N * sizeof(double));

    // daca alocarea esueaza, intoarcem NULL
    if (R == NULL) {
        return NULL;
    }

    // punem in matricea rezultat, elementele rezultate in urma
    // inmultirii elementelor celor doua matrici
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

/*
    Functie care calculeaza produsul dintre o matrice inferior
    triunghiulara si o matrice normala
*/
double* multiply_lower_with_normal(int N, double *LM, double *M) {
    // alocarea memoriei pentru matricea rezultat
    double *R = (double*)malloc(N * N * sizeof(double));
    
    // daca alocarea esueaza, intoarcem NULL
    if (R == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * N + j] = 0;
            // realizam inmultirile doar pentru elementele de sub diagonala
            // principala
            for (int k = 0; k <= i; k++) {
                R[i * N + j] += LM[i * N + k] * M[k * N + j];
            }
        }
    }

    return R;
}

/*
    Functie care calculeaza produsul dintre o matrice normala si 
    una superior triunghiulara
*/
double* multiply_normal_with_upper(int N, double *M, double *UM) {
    // alocarea memoriei pentru matricea rezultat
    double *R = (double*)malloc(N * N * sizeof(double));

    // daca alocarea esueaza, intoarcem NULL
    if (R == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i * N + j] = 0;
            // realizam inmultirile doar pentru elementele de deasupra
            // diagonalei principale
            for (int k = 0; k <= j; k++) {
                R[i * N + j] += M[i * N + k] * UM[k * N + j];
            }
        }
    }

  	return R;
}

/*
    Implementare neoptimizata pentru C = (AT * B + B * A) * BT
 */
double* my_solver(int N, double *A, double* B) {
	double *C;

	double *A_T = transpose(N, A);
	double *B_T = transpose(N, B);

    // AT * B
	double *first_calc = multiply_lower_with_normal(N, A_T, B);

    // B * A
	double *second_calc = multiply_normal_with_upper(N, B, A);

    // (AT * B + B * A)
	double *third_calc = sum(N, first_calc, second_calc);

    // (AT * B + B * A) * BT
	C = multiply_normal(N, third_calc, B_T);

    // eliberam memoria utilizata
	free(A_T);
	free(B_T);
	free(first_calc);
	free(second_calc);
	free(third_calc);

	return C;
}
