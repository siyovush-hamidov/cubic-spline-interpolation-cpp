/**
* @file matrixOperations.h
* @author Melih Altun @2015
**/

#ifndef __MATRIX_OPERATIONS_H__
#define __MATRIX_OPERATIONS_H__

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define lin_index(i, j, numCol)  ( ((i)*(numCol))+(j) )   //2D to 1D array

/*multiply NxK by KxM matrices
parameters: (output) NxM matrix  (input) NxK matrix, KxM matrix, row count for input 1,
col count for input 1 which is also row count for input 2, col count for input 2 */
void multiply_matrices(float output_matrix[], float matrix_1[], float matrix_2[], int N, int K, int M);

/*multiply two NxN matrices
parameters: (output) NxN matrix  (input) NxN input 1, NxN input 2, N */
void multiply_square_matrices(float output_matrix[], float matrix_1[], float matrix_2[], int N);

/*multiply an NxK matrix with a Kx1 vector
parameters: (output) Nx1 output vector  (input) NxK input matrix, Kx1 input vector, K, N  */
void multiply_matrix_with_vector(float output_vector[], float input_matrix[], float input_vector[], int N, int K);

/*add two vectors of the same length
parameters: (output) Nx1 output vector (input) Nx1 input vector 1, Nx1 input vector 2, N */
void add_vectors(float output_vector[], float vector_1[], float vector_2[], int length);

/*subtract a vector from another with the same length
parameters: (output) Nx1 output vector (input) Nx1 input vector 1, Nx1 input vector 2, N */
void subtract_vector_from_another(float output_vector[], float vector_1[], float vector_2[], int length);

/*add two matrices of the same size
parameters: (output) NxM output vector (input) NxM input matrix,  NxM input matrix, N, M */
void add_matrices(float output_matrix[], float matrix_1[], float matrix_2[], int numrows, int numcols);

/*subtract a matrix from another with the same size
parameters: (output) NxM output vector (input) NxM input matrix,  NxM input matrix, N, M */
void subtract_matrix_from_another(float output_matrix[], float matrix_1[], float matrix_2[], int numrows, int numcols);

/* finds matrix transpose of a non square matrix
parameters: (output) MxN matrix transpose (input) NxM matrix, row count, column count */
void transpose(float transposed[], float input[], int N, int M);

/* outputs an identity matrix of given size */
void set_to_identity(float matrix[], int size);

/* copies a matrix into another */
void copy_matrix(float out[], float in[], int N, int M);

/*returns norm of vector*/
float vector_norm(float in[], int size);

/*diagonalizes a vector*/
void vector_to_diagonal(float matrix[], float vector[], int size);

/*extract a vector from square matrix diagonal*/
void diagonal_to_vector(float vector[], float matrix[], int size);

/*scales a vector or matrix by a given coefficient*/
void scale_elements(float arrayOut[], float arrayIn[], float scaler, int size);

/*inner product of two vectors*/
float dotProduct(float vector1[], float vector2[], int size);

/*element-wise multiplication of two vectors or matrices of equal size*/
void multiply_elementwise(float arrayOut[], float arrayIn1[], float arrayIn2[], int size);

/*Checks symmetry of square matrix A of size n*/
bool is_symmetric(float A[], int n);

/*returns trace of a matrix*/
float trace_of_matrix(float A[], int n);

#endif
