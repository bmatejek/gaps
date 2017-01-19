// Include file for Lapack wrapper functions


// Solve an system of equations Ax=b with n variables and n equations
// A is a n by n matrix 
// b has n right hand sides
// x has n variables -- this vector will be filled in with the answer
// nrhs has the number of right hand sides -- i.e., the nuumber of columns in b and x
int RNSolveLinearSystem(double *A, double *x, double *b, int n, int rhs = 1);


// Solve an over-constrained system of equations Ax=b with m variables and n equations
// A is a m by n matrix 
// b has m right hand sides
// x has n variables -- this vector will be filled in with the answer
// nrhs has the number of right hand sides -- i.e., the nuumber of columns in b and x
int RNSolveLeastSquares(double *A, double *x, double *b, int m, int n, int nrhs = 1);


// Perform SVD on matrix a
int RNDecomposeSVD(const double *a, double *u, double *w, double *vt, int m, int n);


// Compute eigenvalues and eigenvectors of a symmetric matrix a
int RNDecomposeEigen(const double *a, double *eigenvalues, double *eigenvectors, int n);
