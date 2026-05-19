#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include "../headers/mesh_generator.h"

void worksplit(int *mystart, int *myend, int myrank, int P, int N){
    int op=N/P;
    int shift=N-op*P;

    if(P<N){
        *mystart=myrank*op;
        *myend=*mystart+op;
        if(shift>myrank){
            *mystart+=myrank;
            *myend+=myrank+1;
        }   
        else{
            *mystart+=shift;
            *myend+=shift;
        }
    }
    else{
        if(myrank<N){
            *mystart =myrank;
            *myend=myrank+1;
        }
        else{
            *mystart=N;
            *myend=N;
        }
    }
}

void generate_mesh(int P, int Ni, double *x, double *xi, double dx, int mystart, double xchange,
                    double *lambdasi,double *lambdas, double lambda1, double lambda2){
    //local mesh
    for(int i=0;i<Ni;i++){
        xi[i] = dx*(mystart+1+i)-dx/2;
        if(xi[i]<xchange)
            lambdasi[i]=lambda1;
        else
            lambdasi[i]=lambda2;
    }

    //global mesh
    allgathering(Ni, P, x, xi);
    allgathering(Ni, P, lambdas, lambdasi);
} 

void allgathering(int send_count, int P, double *X, double *sendbuf){
    //index adjustment
    int *recvcounts=(int*)malloc(P*sizeof(int));
    int *displs=(int*)malloc(P*sizeof(int));
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(&send_count, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
    displs[0]=0;
    for (int i=1;i<P;++i){
        displs[i]=displs[i-1]+recvcounts[i-1];
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgatherv(sendbuf, send_count, MPI_DOUBLE, X, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    free(recvcounts); free(displs);
}