#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "../headers/solver.h"

#define AC_MAT(A, n, i, j) A[(i)*(n)+(j)]


void get_Diagonal_Inv_Square(double *Dinv, double *A, int m) {
    for (int i = 0; i < m; i++)
        Dinv[i] = 1.0 / A[m - 1 + i]; // diagonal central
        
    
}

void get_Negative_Diagonal_Inf_Square(double *L, double *A, int m) {
    for (int i = 0; i < m - 1; i++)
        L[i] = -A[i]; // diagonal inferior
}

void get_Negative_Diagonal_Sup_Square(double *U, double *A, int m) {
    for (int i = 0; i < m - 1; i++)
        U[i] = -A[2 * m - 1 + i]; // diagonal superior
     
}

void Diagonal_Mult(double *Dinv, double *L, double *U, int m, double *T) {
    for (int i = 0; i < m - 1; i++) {
        AC_MAT(T, m - 1, 0, i) = Dinv[i] * U[i];       // upper
        AC_MAT(T, m - 1, 1, i) = Dinv[i + 1] * L[i];   // lower
    }
  
}

void SolveEquationforx(double *T, double *C, double *xk, double *x, int m, double *realdelta) {
    double error = 0.0;
    
    x[0] = AC_MAT(T, m - 1, 0, 0) * xk[1] + C[0];
    error = fabs(x[0] - xk[0]);

    x[m - 1] = AC_MAT(T, m - 1, 1, m - 2) * xk[m - 2] + C[1];
    if (fabs(x[m - 1] - xk[m - 1]) > error)
        error = fabs(x[m - 1] - xk[m - 1]);

    for (int i = 1; i < m - 1; i++) {
        x[i] = AC_MAT(T, m - 1, 0, i) * xk[i + 1] + AC_MAT(T, m - 1, 1, i-1) * xk[i - 1];
        //printf("\n %lf    =   %lf * %lf +%lf * %lf",x[i],AC_MAT(T, m - 1, 0, i),xk[i + 1],AC_MAT(T, m - 1, 1, i-1),xk[i - 1]);
        //fflush(stdout);
        double local_err = fabs(x[i] - xk[i]);
        if (local_err > error)
            error = local_err;
    }
    
    *realdelta = error;
}

void solver_Jacobi_tridiagonal(int m, double *realdelta, double *A, double *x, double *xk, double *b) {

    double *Dinv = (double *)malloc(m * sizeof(double));
    double *L = (double *)malloc((m - 1) * sizeof(double));
    double *U = (double *)malloc((m - 1) * sizeof(double));
    double *T = (double *)malloc(2 * (m - 1) * sizeof(double));
    double C[2];
    int iter = 0;

    get_Diagonal_Inv_Square(Dinv, A, m);
    get_Negative_Diagonal_Inf_Square(L, A, m);
    get_Negative_Diagonal_Sup_Square(U, A, m);
    Diagonal_Mult(Dinv, L, U, m, T);
    
    C[0] = Dinv[0] * b[0];
    C[1] = Dinv[m - 1] * b[1];
   
        SolveEquationforx(T, C, xk, x, m, realdelta);
    

    free(Dinv); free(L); free(U); free(T); 
}
