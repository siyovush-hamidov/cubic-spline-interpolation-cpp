/**
* @file Cholesky.cpp
* @author Melih Altun @2015
**/

#include "../include/Cholesky.h"

/*Cholesky Decomposition on symmetric Matrices: A = L L'
Parameters: (output) transpose of Lower triangular matrix L' (inputs) symmetric matrix square matrix A, number of rows 
return value -1 indicates matrix is not positive definite */
int cholesky_decompose(float *L, float *A, int n)   // cholesky_decompose assumes input matrix is symmetric!
{
	int i, j, k;
	float sum;
	memset(L, 0, n*n*sizeof(float));
	// L_ki = (a_ki - sum(j=1 -> i-1) l_ij l_kj) / l_ii
	// L_kk = sqrt(a_kk - sum(j=1 -> k-1) l_kj^2)
	for (i = 0; i<n; i++) {
		for (j = i; j<n; j++) {
			for (sum = A[i + j*n], k = i - 1; k >= 0; k--)
				sum -= L[i + k*n] * L[j + k*n];
			if (i == j) {
				if (sum <= 0.0)
					return -1;
				L[i + i*n] = sqrt(sum);
			}
			else
				L[j + i*n] = sum / L[i + i*n];
		}
	}
	return 0;
}


/*Matrix inversion with Cholesky Decomposition
Parameters: (output) inverted matrix (inputs) square matrix A, number of rows */
void cholesky_inverse(float *A_inv, float *A, int n)
{
	int i, j, k;
	float sum;
	float *L, *AA, *A_tr, *A_inv_cpy;

	L = new float[n * n];  //c++
	AA = new float[n * n];
	A_tr = new float[n * n];
	A_inv_cpy = new float[n * n];

	bool symmetric = is_symmetric(A, n);
	// if input is not symetric -> A^-1 = A' (A A')^-1
	if (symmetric) {
		if (cholesky_decompose(L, A, n) != 0) {
			memset(A_inv, 0, n*n*sizeof(float));
			//printf("Matrix is not invertible!\n");
			goto CLEAN_UP;
		}
	} else {
		transpose(A_tr, A, n, n);
		multiply_square_matrices(AA, A, A_tr, n);
		if (cholesky_decompose(L, AA, n) != 0) {
			memset(A_inv, 0, n*n*sizeof(float));
			//printf("Matrix is not invertible!\n");
			goto CLEAN_UP;
		}
	}
	// A = L L' -> L L' A^-1 = I -> A^-1 = I / L L' = I (L L')^-1
	for (i = 0; i < n; i++) {
		for (j = 0; j <= i; j++) {
			sum = (i == j ? 1.0 : 0.0);
			for (k = i - 1; k >= j; k--)
				sum -= L[i + k*n] * A_inv[j + k*n];
			A_inv[j + i*n] = sum / L[i + i*n];
		}
	}
	for (i = n - 1; i >= 0; i--) {
		for (j = 0; j <= i; j++) {
			sum = (i < j ? 0.0 : A_inv[j + i*n]);
			for (k = i + 1; k < n; k++)
				sum -= L[k + i*n] * A_inv[j + k*n];
			A_inv[i + j*n] = A_inv[j + i*n] = sum / L[i + i*n];
		}
	}
	if (!symmetric) {
		copy_matrix(A_inv_cpy, A_inv, n, n);
		multiply_square_matrices(A_inv, A_tr, A_inv_cpy,n);
	}

CLEAN_UP:
	delete[] L;
	delete[] AA;
	delete[] A_tr;
	delete[] A_inv_cpy;
}
