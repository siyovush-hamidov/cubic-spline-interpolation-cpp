/**
* @file matrixOperations.cpp
* @author Melih Altun @2015
**/

#include "../include/matrixOperations.h"

/*multiply NxK by KxM matrices
parameters: (output) NxM matrix  (input) NxK matrix, KxM matrix, row count for input 1,
col count for input 1 which is also row count for input 2, col count for input 2 */
void multiply_matrices(float output_matrix[], float matrix_1[], float matrix_2[], int N, int K, int M)
{
	int i, j, k;

	memset(output_matrix, 0, N*M*sizeof(float));

	for (i = 0; i<N; i++) {
		for (j = 0; j<M; j++) {
			for (k = 0; k<K; k++) {
				output_matrix[(i*M) + j] += matrix_1[(i*K) + k] * matrix_2[(k*M) + j];
			}
		}
	}
}

/*multiply two NxN matrices
parameters: (output) NxN matrix  (input) NxN input 1, NxN input 2, N */
void multiply_square_matrices(float output_matrix[], float matrix_1[], float matrix_2[], int  N)
{
	int i, j, k;

	memset(output_matrix, 0, N*N*sizeof(float));

	for (i = 0; i<N; i++) {
		for (j = 0; j<N; j++) {
			for (k = 0; k<N; k++) {
				output_matrix[(i*N) + j] += matrix_1[(i*N) + k] * matrix_2[(k*N) + j];
			}
		}
	}
}

/*multiply an NxK matrix with a Kx1 vector
parameters: (output) Nx1 output vector  (input) NxK input matrix, Kx1 input vector, K, N  */
void multiply_matrix_with_vector(float output_vector[], float input_matrix[], float input_vector[], int N, int K)
{
	int i, j;
	memset(output_vector, 0, N*sizeof(float));

	for (i = 0; i<N; i++) {
		for (j = 0; j<K; j++) {
			output_vector[i] += input_matrix[lin_index(i, j, K)] * input_vector[j];
		}
	}
}

/*add two vectors of the same length
parameters: (output) Nx1 output vector (input) Nx1 input vector 1, Nx1 input vector 2, N */
void add_vectors(float output_vector[], float vector_1[], float vector_2[], int length)
{
	int i;
	for (i = 0; i<length; i++) {
		output_vector[i] = vector_1[i] + vector_2[i];
	}
}

/*subtract a vector from another with the same length
parameters: (output) Nx1 output vector (input) Nx1 input vector 1, Nx1 input vector 2, N */
void subtract_vector_from_another(float output_vector[], float vector_1[], float vector_2[], int length)
{
	int i;
	for (i = 0; i<length; i++) {
		output_vector[i] = vector_1[i] - vector_2[i];
	}
}

/*add two matrices of the same size
parameters: (output) NxM output vector (input) NxM input matrix,  NxM input matrix, N, M */
void add_matrices(float output_matrix[], float matrix_1[], float matrix_2[], int numrows, int numcols)
{
	int i, j;
	for (i = 0; i<numrows; i++) {
		for (j = 0; j<numcols; j++) {
			output_matrix[lin_index(i, j, numcols)] = matrix_1[lin_index(i, j, numcols)] + matrix_2[lin_index(i, j, numcols)];
		}
	}
}

/*subtract a matrix from another with the same size
parameters: (output) NxM output vector (input) NxM input matrix 1,  NxM input matrix 2, N, M */
void subtract_matrix_from_another(float output_matrix[], float matrix_1[], float matrix_2[], int numrows, int numcols)
{
	int i, j;
	for (i = 0; i<numrows; i++) {
		for (j = 0; j<numcols; j++) {
			output_matrix[lin_index(i, j, numcols)] = matrix_1[lin_index(i, j, numcols)] - matrix_2[lin_index(i, j, numcols)];
		}
	}
}

/* finds matrix transpose of a non square matrix
parameters: (output) MxN matrix transpose (input) NxM matrix, row count, column count */
void transpose(float transposed[], float *input, int N, int M)
{
	int   i, j;
	for (i = 0; i < N; i++)	{
		for (j = 0; j < M; j++) {
			transposed[lin_index(j, i, N)] = input[lin_index(i, j, M)];
		}
	}
}

/* outputs an identity matrix of given size
parameters: (output) identity matrix, (input) matrix size */
void set_to_identity(float matrix[], int size)
{
	int i, j;
	for (i = 0; i < size; i++){
		for (j = 0; j < size; j++) {
			if (i == j)
				matrix[lin_index(i, j, size)] = 1.0;
			else
				matrix[lin_index(i, j, size)] = 0.0;
		}
	}
}

/* copies a matrix into another
parameters: (output) matrix, (input) matrix, row count, column count */
void copy_matrix(float out[], float in[], int N, int M)
{
	int i;
	for (i = 0; i < M*N; i++)
		out[i] = in[i];
}

/*returns norm of vector
parameters: (input) vector, number of elements */
float vector_norm(float in[], int size) {
	int i;
	float sum = 0;
	for (i = 0; i < size; i++)
		sum += in[i] * in[i];
	return (float)sqrt(sum);   // <--- sqrt
}

/*diagonalizes a vector
parameters: (output) matrix, (input) vector, number of elements */
void vector_to_diagonal(float matrix[], float vector[], int size)
{
	int i, j;
	for (i = 0; i < size; i++){
		for (j = 0; j < size; j++) {
			if (i == j)
				matrix[lin_index(i, j, size)] = vector[i];
			else
				matrix[lin_index(i, j, size)] = 0;
		}
	}
}

/*extracts a vector from square matrix diagonal
parameters: (output)vector, (input)matrix, number of elements */
void diagonal_to_vector(float vector[], float matrix[], int size)
{
	int i, j;
	for (i = 0; i < size; i++){
		for (j = 0; j < size; j++){
			if (i == j)
				vector[i] = matrix[lin_index(i, j, size)];
		}
	}
}

/*scales a vector or matrix by a given coefficient
parameters: (output)vector or matrix, (input) vector or matrix, scaling coefficient, number of elements */
void scale_elements(float arrayOut[], float arrayIn[], float scaler, int size)
{
	int i;
	for (i = 0; i < size; i++)
		arrayOut[i] = arrayIn[i] * scaler;
}

/*inner product of two vectors
parameters: (inputs) vector1, vector2, number of elements
returns dot product*/
float dotProduct(float vector1[], float vector2[], int size)
{
	int i;
	float sum = 0;
	for (i = 0; i < size; i++)
		sum += vector1[i] * vector2[i];
	return sum;
}

/*element-wise multiplication of two vectors or matrices of equal size
parameters: (output) output matrix or vector, (input) 1st matrix or vector, 2nd matrix or vector, number of elements*/
void multiply_elementwise(float arrayOut[], float arrayIn1[], float arrayIn2[], int size)
{
	int i;
	for (i = 0; i < size; i++)
		arrayOut[i] = arrayIn1[i] * arrayIn2[i];
}

/*Checks symmetry of square matrix
parameters: (input) square matrix A of size n
returns symmetry result */ 
bool is_symmetric(float *A, int n)
{
	bool symmetric = true;
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			if (A[lin_index(i, j, n)] != A[lin_index(j, i, n)]) {
				symmetric = false;
				break;
			}
		}
		if (!symmetric)
			break;
	}
	return symmetric;
}

/*returns trace of a matrix
parameters: (input) input matrix, matrix size */
float trace_of_matrix(float A[], int n)
{
	int i;
	float sum = 0;
	for (i = 0; i < n; i++)
		sum += A[lin_index(i, i, n)];
	return sum;
}
