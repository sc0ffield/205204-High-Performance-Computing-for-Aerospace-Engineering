#include <stdlib.h>
#include <stdio.h>
#include "../headers/tests.h"

void test_meshing(int N, double L, double dx, double lambda1, double lambda2, double xchange){ 
    double *x = (double*)malloc((N+2)*sizeof(double));
    double *x_cv = (double*)malloc((N+2)*sizeof(double));
    double *lambda = (double*)malloc((N+2)*sizeof(double));
    
    printf("X CONTROL VOLUMES (TEST): ");
    for(int i=0; i<N+2; i++){
        x_cv[i] = dx*i;
        printf("%f ", x_cv[i]);
    }
    printf("\n");
    
    printf("X NODES (TEST): ");
    for(int i=0; i<(N+2); i++){
        x[i] = -dx/2+i*dx;
        printf("%f ", x[i]);
    }
    printf("\n");

    printf("LAMBDA X (TEST): ");
    for(int i=0; i<(N+2); i++){
        if(x[i]<xchange)
            lambda[i]=lambda1;
        else
            lambda[i]=lambda2;
        printf("%f ", lambda[i]);
    }
    printf("\n");

    free(x); free(x_cv); free(lambda);
}