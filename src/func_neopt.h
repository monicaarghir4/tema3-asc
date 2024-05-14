#ifndef FUNC_NEOPT_
#define FUNC_NEOPT_

#include <stdlib.h>
#include <stdio.h>

double* transpose(int N, double *M);
double* sum(int N, double *M1, double *M2);
double* multiply_normal(int N, double *M1, double *M2);
double* multiply_lower_with_normal(int N, double *LM, double *M);

#endif  // FUNC_NEOPT_
