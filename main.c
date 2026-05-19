#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#include "headers/mesh_generator.h"
#include "headers/pre_solver.h"
#include "headers/solver.h"
#include "headers/print.h"
#include "headers/tests.h"

//CHECK MPI STATUS
void checkr(int r,char *txt) {
    if (r!=MPI_SUCCESS) {
        fprintf(stderr,"Error: %s\n",txt);
        exit(-1);
    }   
}

int main(int argc, char **argv) {
    int myrank, mystart, myend, P, r, l; //variables for parallelizing

    //PROBLEM KNOWN CONSTANTS
    const int N = 1000; const double L = 5.0, dx = L/N, xchange = 2.8;
    const double T0 = 100, TL = 200, lambda1 = 120, lambda2 = 180;

    const double halo_left = -dx/2, halo_right = L + dx/2;
    const double delta = 1e-5, T_ini = (T0+TL)*0.5; int Ni;
    
    r=MPI_Init(&argc,&argv); checkr(r,"init");
    r=MPI_Comm_rank(MPI_COMM_WORLD,&myrank); checkr(r,"rank");
    r=MPI_Comm_size(MPI_COMM_WORLD,&P); checkr(r,"size");

    //1.MESHING
    worksplit(&mystart, &myend, myrank, P, N);
    Ni = myend-mystart; //number of nodes for each process
    double *x = (double*)malloc(N*sizeof(double)); //node positions 
    double *xi = (double*)malloc(Ni*sizeof(double)); //node positions for each process
    double *lambdas = (double*)malloc(N*sizeof(double)); //lambdas for each node
    double *lambdasi = (double*)malloc(Ni*sizeof(double)); //lambdas for each node of each process
    generate_mesh(P, Ni, x, xi, dx, mystart, xchange, lambdasi, lambdas, lambda1, lambda2);
    

    free (xi); free(lambdasi);

    //2.INITIAL FIELD 
    double *Ti_s = (double*)malloc((Ni+2)*sizeof(double)); //supposed real node temperatures for each process
    double *Ti_c = (double*)malloc((Ni+2)*sizeof(double)); //calculated real nodes temperatures for each process

    for(int i=0;i<Ni+2;i++)
        Ti_s[i]=T_ini; //Supposed temperatures

    //3.BOUNDARY CONDITIONS
    double b[2];

    if(myrank==0)
    Ti_s[0]=T0;
    if(myrank==P-1)
    Ti_s[Ni+1]=TL;

    //pre iteration
    double *a = (double*)malloc((3*Ni+4)*sizeof(double)); //a represents the 3 diagonals from tridiagonal coefficient matrix(from down to up)
   
    double realdelta, global_delta, dub_time, time1, time2, elapsedg=0;
    int iter = 1;

    //4.ITERATION
    time1=MPI_Wtime();
    while(1){
        //4.1.HALO UPDATE
       
         

        halo_update(myrank,Ni,P,Ti_s);
        
       

        b[0]=Ti_s[0];   b[1]=Ti_s[Ni+1];

        //4.2.COEFFICIENTS
        calc_coefficients(lambdas, lambda1, lambda2, a, Ni, P, mystart, myrank);
    
        //5.2.JACOBI SOLVER

       
        solver_Jacobi_tridiagonal(Ni+2,&realdelta, a, Ti_c, Ti_s,b);
      
        for(int i=0; i<Ni+2; i++)
            Ti_s[i]=Ti_c[i]; //update real nodes
        
        //5.3.CONVERGENCE CRITERION
         

        MPI_Allreduce(&realdelta, &global_delta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if(global_delta<delta)
            break;

        iter++;

    }
    time2=MPI_Wtime();

    dub_time=time2-time1;
    MPI_Reduce(&dub_time, &elapsedg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
   



    free(Ti_c);free(a);
    
    //6.PRINTING
    
    MPI_Barrier(MPI_COMM_WORLD);

    double *T_global = NULL, *Ti_f;
    Ti_f = malloc(Ni* sizeof(double));

    if (myrank == 0) 
        T_global = malloc(N * sizeof(double));

    for(int i=1; i<Ni+1; i++)
        Ti_f[i-1]=Ti_s[i];

    reunir_temperaturas(Ni, Ti_f, T_global, N, MPI_COMM_WORLD);
    

    free(Ti_s); free(Ti_f);
    
    if(myrank==0){
        print_statement(L,xchange,lambda1,lambda2,T0,TL);
        print_discretization(N,dx,x,lambdas, halo_left, halo_right);
        print_initialf(T_ini);
        print_bc(TL, T0);
        print_iteration();
        print_results(T_global, N, iter,delta);
        free(T_global);
    }
    
    free(x); free(lambdas);

    if(myrank==0)
    printf("Time-----------------> %lf\n",dub_time);

    MPI_Finalize();
    /*
    SERIAL TEST AREA 
    test_meshing(N,L,dx,lambda1,lambda2,xchange);
    */
}