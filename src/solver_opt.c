/*
 * Tema 3 ASC
 * 2024 Spring
 */
#include "utils.h"

/*
	Functie care calculeaza transpusa unei matrici
*/
double* transpose(int N, double *M) {
	// alocarea memoriei pentru matricea rezultat
    double *T = (double*)malloc(N * N * sizeof(double));

	// daca alocarea esueaza, intoarcem NULL
    if (T == NULL) {
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
			// punem indecsii folositi in registre
			register int idx1 = i * N + j;
			register int idx2 = j * N + i;
            T[idx1] = M[idx2];
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

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
			// stocam indexul intr-un registru
			register int idx = i * N + j;
            R[idx] = M1[idx] + M2[idx];
        }
    }

    return R;
}

/*
	Functie care calculeaza produsul a doua matrici
*/
double* multiply_normal(int N, double *M1, double *M2) {
	// alocarea memoriei pentru matricea rezultat si initializarea cu 0
    double *R = (double*)calloc(N * N, sizeof(double));

	// daca alocarea esueaza, intoarcem NULL
    if (R == NULL) {
        return NULL;
    }

	// loop unrolling
	for (int i = 0; i < N; i++) {
		// pointer la startul liniei din prima matrice
		double *orig_p1 = &M1[i * N];

		for (int j = 0; j < N; j++) {
			// pointeri pentru traversarea coloanelor celor 2 matrice
			double *p1 = orig_p1;
			double *p2 = &M2[j];

			// stocarea sumelor partiale in 4 registrii
			register double r0 = 0.0;
			register double r1 = 0.0;
			register double r2 = 0.0;
			register double r3 = 0.0;

			int k = 0;

			while (k < N - 3) {
				// calcularea sumelor pentru cate 4 elemente
				r0 += p1[0] * p2[0 * N];
				r1 += p1[1] * p2[1 * N];
				r2 += p1[2] * p2[2 * N];
				r3 += p1[3] * p2[3 * N];
	
				// mutam pointerii pentru urmatorul set de elemente
				p1 += 4;
				p2 += 4 * N;

				// incrementam counter-ul pentru urmatoarele 4 elemente
				k += 4;
			}

			// calcularea sumelor pentru restul elementelor
			while (k < N) {
				r0 += *p1 * *p2;

				p1++;
				p2 += N;

				k++;
			}

			// stocarea sumei finale in matricea rezultat
			R[i * N + j] = r0 + r1 + r2 + r3;
		}
	}

    return R;
}

/*
	Functie care calculeaza produsul dintre o matrice inferior
    triunghiulara si o matrice normala
*/
double* multiply_lower_with_normal(int N, double *LM, double *M) {
	// alocarea memoriei pentru matricea rezultat si initializarea cu 0
    double *R = (double*)calloc(N * N, sizeof(double));

	// daca alocarea esueaza, intoarcem NULL
    if (R == NULL) {
        return NULL;
    }

	// loop unrolling
	for (int i = 0; i < N; i++) {
		// pointer la startul liniei din prima matrice
		double *orig_p1 = &LM[i * N];

		for (int j = 0; j < N; j++) {
			// pointeri pentru traversarea coloanelor celor 2 matrice
			double *p1 = orig_p1;
			double *p2 = &M[j];

			// stocarea sumelor partiale in 4 registrii
			register double r0 = 0.0;
			register double r1 = 0.0;
			register double r2 = 0.0;
			register double r3 = 0.0;

			int k = 0;

			// parcurgerea elementelor
			while (k <= i - 3) {
				// calcularea sumelor pentru cate 4 elemente
				r0 += p1[0] * p2[0 * N];
				r1 += p1[1] * p2[1 * N];
				r2 += p1[2] * p2[2 * N];
				r3 += p1[3] * p2[3 * N];

				// mutam pointerii pentru urmatorul set de elemente
				p1 += 4;
				p2 += 4 * N;

				// incrementam counter-ul pentru urmatoarele 4 elemente
				k += 4;
			}
			
			// calcularea sumelor pentru restul elementelor
			// pana la valoarea i pentru a nu lua in calcul
			// elementele egale cu 0
			while (k <= i) {
				r0 += *p1 * *p2;

				p1++;
				p2 += N;

				k++;
			}

			// stocarea sumei finale in matricea rezultat
			R[i * N + j] = r0 + r1 + r2 + r3;
		}
	}

    return R;
}

/*
	Functie care calculeaza produsul dintre o matrice normala si 
    una superior triunghiulara
*/
double* multiply_normal_with_upper(int N, double *M, double *UM) {
	// alocarea memoriei pentru matricea rezultat si initializarea cu 0
    double *R = (double*)calloc(N * N, sizeof(double));

	// daca alocarea esueaza, intoarcem NULL
	if (R == NULL) {
        return NULL;
    }

	// loop unrolling
	for (int i = 0; i < N; i++) {\
		// pointer la startul liniei din prima matrice
		double *orig_p1 = &M[i * N];

		for (int j = 0; j < N; j++) {
			// pointeri pentru traversarea coloanelor celor 2 matrice
			double *p1 = orig_p1;
			double *p2 = &UM[j];

			// stocarea sumelor partiale in 4 registrii
			register double r0 = 0.0;
			register double r1 = 0.0;
			register double r2 = 0.0;
			register double r3 = 0.0;

			int k = 0;

			// parcurgerea elementelor
			while (k <= j - 3) {
				// calcularea sumelor pentru cate 4 elemente
				r0 += p1[0] * p2[0 * N];
				r1 += p1[1] * p2[1 * N];
				r2 += p1[2] * p2[2 * N];
				r3 += p1[3] * p2[3 * N];

				// mutam pointerii pentru urmatorul set de elemente
				p1 += 4;
				p2 += 4 * N;

				// incrementam counter-ul pentru urmatoarele 4 elemente
				k += 4;
			}

			// calcularea sumelor pentru restul elementelor
			// pana la valoarea j pentru a nu lua in calcul
			// elementele egale cu 0
			while (k <= j) {
				r0 += *p1 * *p2;

				p1++;
				p2 += N;

				k++;
			}

			// stocarea sumei finale in matricea rezultat
			R[i * N + j] = r0 + r1 + r2 + r3;
		}
	}

  	return R;
}

/*
	Implementare optimizata pentru C = (AT * B + B * A) * BT
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
