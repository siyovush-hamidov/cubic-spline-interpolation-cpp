/**
* @file Cholesky.h
* @author Melih Altun @2015
**/


#ifndef __CHOLESKY_H__
#define __CHOLESKY_H__

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "matrixOperations.h"

/*Cholesky Decomposition on symmetric Matrices: A = L L'
Parameters: (output) transpose of Lower triangular matrix L' (inputs) symmetric matrix square matrix A, number of rows
return value -1 indicates matrix is not positive definite */
int cholesky_decompose(float *L, float *A, int n);

/*Matrix inversion with Cholesky Decomposition
Parameters: (output) inverted matrix (inputs) square matrix A, number of rows */
void cholesky_inverse(float *A_inv, float *A, int n);

#endif
