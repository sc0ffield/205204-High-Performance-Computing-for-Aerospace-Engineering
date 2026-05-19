#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "mpi.h"

#define AC_MAT(A, n, i, j) A[(i)*(n)+(j)]

void solver_Jacobi_tridiagonal(int m, double delta, double *A, double *x, double *b);
void Initialize_xk(double *xk, int m, double initial_value);
void get_Diagonal_Inv_Square(double *Dinv, double *A, int m);
void get_Negative_Diagonal_Inf_Square(double *L, double *A, int m);
void get_Negative_Diagonal_Sup_Square(double *U, double *A, int m);
void Diagonal_Mult(double *Dinv, double *L, double *U, int m, double *T);
void SolveEquationforx(double* T, double *C, double *xk, double *x, int m, double *realdelta);
void printmatrix(double *A, int m, int n);  
void equal_vector(double *A, double *B, int m);

int main(){
    int m=2; 

    double *T,*Coeff, *BC, delta=10e-9;

    T=(double*)malloc(m*sizeof(double));
    Coeff=(double*)malloc((3*m-2)*sizeof(double));
    BC=(double*)malloc(2*sizeof(double));

    Coeff[0]=5;
    Coeff[1]=2;
    Coeff[2]=7;
    Coeff[3]=1;

    BC[0]=11;
    BC[1]=13;

    solver_Jacobi_tridiagonal(m,delta,Coeff,T,BC);

    printmatrix(T,m,1); //the solution has to be: {7.111, -3.222}
}

void solver_Jacobi_tridiagonal(int m, double delta, double *A, double *x, double *b){   //A represents the 3 diagonals from  
                                                                                        //tridiagonal coefficient matrix(from down to up)
    double *Dinv, *L, *U, *T, *C, *xk, realdelta=10, initial_value=1;                   //Adapt diagonals!!
                                                                                        //b is for Halo. Matrix 2xNproc.

    Dinv=(double*)malloc(m*sizeof(double));
    L=(double*)malloc((m-1)*sizeof(double));
    U=(double*)malloc((m-1)*sizeof(double));
    T=(double*)malloc(2*(m-1)*sizeof(double));
    C=(double*)malloc(2*sizeof(double));
    xk=(double*)malloc(m*sizeof(double));

    Initialize_xk(xk,m,initial_value);
    get_Diagonal_Inv_Square(Dinv,A,m);
    get_Negative_Diagonal_Inf_Square(L,A,m);
    get_Negative_Diagonal_Sup_Square(U,A,m);

    free(A);

    Diagonal_Mult(Dinv,L,U,m,T);

    free(U);
    free(L);

    C[0]=Dinv[0]*b[0];
    C[1]=Dinv[m-1]*b[1]; 

    free(Dinv);

    int iter=0;
    while(realdelta>delta){
        iter++;

        SolveEquationforx(T,C,xk,x,m,&realdelta);
        equal_vector(xk,x,m);
    }
}

void Initialize_xk(double *xk, int m, double initial_value){
    for(int i=0; i<m; i++)
        xk[i]=initial_value;
}

void get_Diagonal_Inv_Square(double *Dinv, double *A, int m){
    for(int i=0; i<m; i++)
        Dinv[i]=1/A[m-1+i];
}

void get_Negative_Diagonal_Inf_Square(double *L, double *A, int m){
    for(int i=0; i<m-1; i++)
        L[i]=-A[i];
}

void get_Negative_Diagonal_Sup_Square(double *U, double *A, int m){
    for(int i=0; i<m-1; i++)
        U[i]=-A[2*m-1+i];
}

void Diagonal_Mult(double *Dinv, double *L, double *U, int m, double *T){
    for(int i=0; i<m-1; i++){
        AC_MAT(T,m-1,0,i)=Dinv[i]*U[i];   //Upper diagonal
        AC_MAT(T,m-1,1,i)=Dinv[i+1]*L[i];     //Lower diagonal
    }
}

void SolveEquationforx(double* T, double *C, double *xk, double *x, int m, double *realdelta){
    double error;
    
    x[0]=AC_MAT(T,m-1,0,0)*xk[1]+C[0];
    error=fabs(x[0]-xk[0]); 

    x[m-1]=AC_MAT(T,m-1,1,m-2)*xk[m-2]+C[1];
    if(fabs(x[m-1]-xk[m-1])>error)
        error=fabs(x[m-1]-xk[m-1]);
    for(int i=1; i<m-1; i++){
        x[i]=AC_MAT(T,m-1,0,i)*xk[i+1]+AC_MAT(T,m-1,1,i)*xk[i-1];
        if(error<(fabs(x[i]-xk[i])))
            error=fabs(x[i]-xk[i]);
    }
    *realdelta=error;
}

void printmatrix(double *A, int m, int n){        
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            printf("%.2lf    ", AC_MAT(A,n,i,j));   
        }
        printf("\n");
    }
    fflush(stdout);
 
}

void equal_vector(double *A, double *B, int m){
    for(int i=0; i<m; i++)
        A[i]=B[i];
}