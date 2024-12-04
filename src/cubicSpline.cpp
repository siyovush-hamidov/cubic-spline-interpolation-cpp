#include "../include/cubicSpline.h"

void cubicSpline(float c0[], float c1[], float c2[], float c3[], float x[], float y[], int n)
{
	int i;
	float *diff_x, *diff_y, *A, *A_inv, *m, *b;
	diff_x = new float[n - 1];
	diff_y = new float[n - 1];
	A = new float[(n - 2)*(n - 2)];
	A_inv = new float[(n - 2)*(n - 2)];
	b = new float[n - 2];
	m = new float[n];

	memset(A, 0, (n - 2)*(n - 2)*sizeof(float));

	for (i = 0; i < n - 1; i++) {
		diff_x[i] = x[i + 1] - x[i];
		diff_y[i] = (y[i + 1] - y[i]) / diff_x[i];
	}


	for (i = 0; i < n - 2; i++) {
		if (i>0)
			A[lin_index(i, i - 1, n-2)] = diff_x[i];
		A[lin_index(i, i, n-2)] = 2 * (diff_x[i] + diff_x[i + 1]);
		if (i < n - 3)
			A[lin_index(i, i + 1, n-2)] = diff_x[i + 1];
		b[i] = 6 * (diff_y[i + 1] - diff_y[i]);
	}

	cholesky_inverse(A_inv, A, n - 2);
	multiply_matrix_with_vector(m, A_inv, b, n - 2, n - 2);

	for (i = n - 2; i>0; i--)
		m[i] = m[i - 1];
	m[0] = m[n - 1] = 0;

	for (i = 0; i < n-1; i++) {
		c0[i] = y[i];
		c1[i] = diff_y[i] - diff_x[i] * (2 * m[i] + m[i + 1]) / 6;
		c2[i] = m[i] / 2;
		c3[i] = (m[i + 1] - m[i]) / (6 * diff_x[i]);
	}

	delete[] diff_x;
	delete[] diff_y;
	delete[] A;
	delete[] A_inv;
	delete[] m;
	delete[] b;
}

float evaluate_segment(float x, float coeffs_segment[]) {
	float y = coeffs_segment[3] * pow(x, 3) + coeffs_segment[2] * pow(x, 2) + coeffs_segment[1] * x + coeffs_segment[0];
	return y;
}

int getSegmentIndex(float x1[], int n1, float x)
{
	for (int i = n1 - 1; i >= 0; i--) {
		if (x >= x1[i])
			return i;
	}
	return 0;
}

void evaluate_polynomial(float y2[], float x2[], int n2, float x1[], int n1, float c0[], float c1[], float c2[], float c3[])
{
	float x;
	int n_segments, segment_index, start_index, end_index;
	float segment_coeffs[4];
	n_segments = n1 - 1;
	for (int i = 0; i < n2; i++) {
		x = x2[i];
		segment_index = getSegmentIndex(x1, n1, x);

		if (segment_index == n_segments)
			segment_index = n_segments - 1;

		segment_coeffs[0] = c0[segment_index];
		segment_coeffs[1] = c1[segment_index];
		segment_coeffs[2] = c2[segment_index];
		segment_coeffs[3] = c3[segment_index];
	
		y2[i] = evaluate_segment(x-x1[segment_index], segment_coeffs);
	}
}