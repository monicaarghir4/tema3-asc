#ifndef FUNC_OPT_
#define FUNC_OPT_

#include <stdlib.h>
#include <stdio.h>

double* transpose(int N, double *M);
double* sum(int N, double *M1, double *M2);
double* multiply_normal(int N, double *M1, double *M2);
double* multiply_lower_with_normal(int N, double *LM, double *M);

#endif  // FUNC_OPT_
