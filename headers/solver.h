#ifndef SOLVER_H
#define SOLVER_H

void get_Diagonal_Inv_Square(double *Dinv, double *A, int m);
void get_Negative_Diagonal_Inf_Square(double *L, double *A, int m);
void get_Negative_Diagonal_Sup_Square(double *U, double *A, int m);
void Diagonal_Mult(double *Dinv, double *L, double *U, int m, double *T);
void SolveEquationforx(double *T, double *C, double *xk, double *x, int m, double *realdelta);
void solver_Jacobi_tridiagonal(int m, double *realdelta, double *A, double *x, double *xk, double *b);

#endif