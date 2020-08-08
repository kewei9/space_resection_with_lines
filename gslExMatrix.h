
#ifndef _GSL_EX_MATRIX_OPERATION_FILE_HEADER_
#define _GSL_EX_MATRIX_OPERATION_FILE_HEADER_


// matrix operations
int gslExMultiplyMatrix(double* matrix_A, int row_A, int col_A, double* matrix_B, int row_B, int col_B, double* result_matrix);

int gslExInverseMatrix(double* matrix, int nDim);

int gslExSVDDecompMatrix(double* matrix_A, int row_A, int col_A,
						 double* matrix_U, double* matrix_D, double* matrix_V);// row_A >= col_A

void gslExLeastSquareFit_Linear(double* X, double* y, double* c, int row, int col);
//

#endif