#include <stdlib.h>
#include <stdio.h>
#include "../headers/print.h"

void print_statement(double L, double xchange, double lambda1, 
                        double lambda2, double T0, double TL){
    printf("1D HEAT CONDUCTION PROBLEM IN A COMPOSITE WALL\n\n"
            "PROBLEM STATEMENT\n\n");
    printf("A wall of total length %fm is composed of "
            "two materials joined at x_c = %fm. The first segment "
            "(0<=x<x_c) is made of a material A with thermal conductivity "
            "lambda_1 = %fW/mK, and the second segment (x_c<=x<L) is made "
            "of material B with thermal conductivity lambda_1 = %fW/mK.\n\n"
            "The system is in steady-state, and heat conduction is "
            "one-dimensional along the x-axis. The temperature at "
            "the left end of the wall is maintained at T_0 = %fK, and the "
            "temperature at the right end is maintained at T_L = %fK.\n\n",
            L,xchange,lambda1,lambda2,T0,TL);
    printf("Assuming:\n-Steady-state heat conduction\n-No internal heat " 
            "generation\n-Perfect thermal contact between the two materials "
            "\n-Constant thermal conductivity in each material\n-One-dimensional "
            "conduction only.\n\nDetermine the temperature distribution T(x) along "
            "the lenght of the wall.\n");
}

void print_discretization(int N, double dx, double* x, double *lambdas, double halo_left, double halo_right){
    printf("\n\n1. DISCRETIZATION OF THE MESH\n\n");
    printf("Contemplating %d internal nodes:\ndx = %f \nNode positions [m] = {",N,dx);
    for(int i=0;i<N;i++){
        printf("%f",x[i]);
        if(i<(N-1))
            printf(", ");
    }
    printf("}\nThermal conductivities per node [W/mK] = {");
    for(int i=0;i<N;i++){
        printf("%f",lambdas[i]);
        if(i<(N-1))
            printf(", ");
    }
    printf("}\n\nTwo external nodes (HALO) are also considered at x[0] = "
            "%fm and x[N+1] = %fm.\n", halo_left, halo_right);
}

void print_initialf(double Tini){
    printf("\n\n2. INITIAL MAP DEFINITION\n\n");
    printf("T_supposed[i]=T_initial=%fK\n",Tini);
}

void print_bc(double TL, double T0){
    printf("\n\n3. BOUNDARY CONDITIONS\n\n");
    printf("T[0]=T_0=%fK\n",T0);
    printf("T[N+1]=T_L=%fK\n",TL);
}

void print_iteration(){
    printf("\n\n4. ITERATION\n\n");
    printf("   4.1. Halo update.\n   4.2. Calculation of coefficients (aP,aW,aE,b)."
            "\n   4.3. Map resolution by the Jacobi method.\n"
            "   4.4. Is max|T_calculated[i] - T_supposed[i]|<delta?\n       NO: "
            "T_supposed[i] = T_calculated[i] and back to point 4.1.\n       YES: end of the iteration.");
}

void print_results(double *T, int N, int iter, double delta){
    printf("\n\n6. RESULTS\n\n");
    printf("Contemplating delta = %lf, and with %d iterations:\n",delta,iter);
    printf("T[i] [K] (internal nodes)= {");
    for(int i=0;i<N;i++){
        printf("%f",T[i]);
        if(i<(N-1))
            printf(", ");
    }
    printf("}\n");
}