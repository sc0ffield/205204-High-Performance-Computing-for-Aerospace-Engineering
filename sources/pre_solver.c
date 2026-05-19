#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include "../headers/pre_solver.h"

void halo_update(int myrank, int Ni, int Nproc, double *Ti_s) {
    
                
            
                
                
                if (myrank > 0) {
                    // Send first node to left
                    MPI_Send(&Ti_s[1], 1, MPI_DOUBLE,myrank - 1, 0, MPI_COMM_WORLD);
                }
                if (myrank < Nproc - 1) {
                    // Left halo recieves from right
                    MPI_Recv(&Ti_s[Ni+1], 1, MPI_DOUBLE,myrank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                if (myrank < Nproc - 1) {
                    // Send last node to right
                    MPI_Send(&Ti_s[Ni], 1, MPI_DOUBLE,myrank + 1, 1, MPI_COMM_WORLD);
                }
                if (myrank > 0) {
                    // Right halo recieves from left
                    MPI_Recv(&Ti_s[0], 1, MPI_DOUBLE,myrank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

             }

double lambdaharm(double lambdai, double lambdaj){
    return ((2*lambdai*lambdaj)/(lambdai+lambdaj));
}

void calc_coefficients(double *lambdas, double lambda1, double lambda2, double *a, int Ni, int P, 
                        int mystart, int myrank){      

    double* lambdaE = (double*)malloc(Ni*sizeof(double));
    double* lambdaW = (double*)malloc(Ni*sizeof(double));
    double* lambdaP = (double*)malloc(Ni*sizeof(double));

    for(int i=0;i<Ni;i++){
        lambdaW[i] = (myrank==0&&(i==0))?lambda1:lambdaharm(lambdas[i+mystart], lambdas[i+mystart-1]);
        lambdaE[i] = (myrank==P-1&&(i==Ni-1))?lambda2:lambdaharm(lambdas[i+mystart], lambdas[i+mystart+1]);
        lambdaP[i] = (lambdaW[i] + lambdaE[i]);
    }

    int idx = 0;

    for(int i = 0; i < Ni; i++) 
        a[idx++] = -lambdaW[i];

    a[idx++]=0;
    a[idx++]=1;

    for(int i = 0; i < Ni; i++)
        a[idx++] = lambdaP[i]; 

    a[idx++]=1;
    a[idx++]=0;

    for(int i = 0; i < Ni; i++)
        a[idx++] = -lambdaE[i];

    free(lambdaE); free(lambdaW); free(lambdaP);
}