
#include "stdafx.h"
#include "gslExMatrix.h"


int gslExMultiplyMatrix(double* matrix_A, int row_A, int col_A, double* matrix_B, int row_B, int col_B, double* result_matrix)
{
	gsl_matrix_view A_v = gsl_matrix_view_array(matrix_A, row_A, col_A);
	gsl_matrix_view B_v = gsl_matrix_view_array(matrix_B, row_B, col_B);

	gsl_matrix_view result_v = gsl_matrix_view_array(result_matrix, row_A, row_B);

	int res = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0f, &A_v.matrix, &B_v.matrix, 0.0, &result_v.matrix);

	return res;
}

int gslExSVDDecompMatrix(double* matrix_A, int row_A, int col_A,
						 double* matrix_U, double* matrix_D, double* matrix_V)// row_A >= col_A
{
	memcpy(matrix_U, matrix_A, sizeof(double) * row_A * col_A);

	gsl_matrix_view A_v = gsl_matrix_view_array(matrix_U, row_A, col_A);
	gsl_matrix_view V_v = gsl_matrix_view_array(matrix_V, col_A, col_A);

	double* temp = new double[col_A];
	gsl_vector_view S_v = gsl_vector_view_array(temp, col_A);

	double* work = new double[col_A];
	gsl_vector_view Work_v = gsl_vector_view_array(work, col_A);

	gsl_linalg_SV_decomp(&A_v.matrix, &V_v.matrix, &S_v.vector, &Work_v.vector);

	memset(matrix_D, 0x00, sizeof(double) * col_A);
	for (int i = 0; i < col_A; ++ i)
		matrix_D[i * col_A + i] = temp[i];

	delete[] work;
	delete[] temp;

	return 1;
}

int gslExInverseMatrix(double* matrix, int nDim)
{
	gsl_matrix_view A = gsl_matrix_view_array(matrix, nDim, nDim);

	double* b = new double[nDim * nDim];
	gsl_matrix_view B = gsl_matrix_view_array(b, nDim, nDim);

	int sig = 0;
	gsl_permutation* perm = gsl_permutation_alloc(nDim);

	int nRes = 0;
	nRes += gsl_linalg_LU_decomp(&A.matrix, perm, &sig);
	nRes += gsl_linalg_LU_invert(&A.matrix, perm, &B.matrix);

	gsl_permutation_free(perm);

	memcpy(matrix, b, sizeof(double) * nDim * nDim); 			// M = M(-) 

	delete[] b;

	return nRes;
}

void gslExLeastSquareFit_Linear(double* X, double* y, double* c, int row, int col)
{
	gsl_matrix_view X_v = gsl_matrix_view_array(X, row, col);
	gsl_vector_view y_v = gsl_vector_view_array(y, row);

	gsl_vector_view c_v = gsl_vector_view_array(c, col);

	double* cov = new double[col * col];
	gsl_matrix_view cov_v = gsl_matrix_view_array(cov, col, col);

	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(row, col);

	double chisq;
	gsl_multifit_linear(&X_v.matrix, &y_v.vector, &c_v.vector, &cov_v.matrix, &chisq, work);

	gsl_multifit_linear_free(work);

	delete[] cov;
}