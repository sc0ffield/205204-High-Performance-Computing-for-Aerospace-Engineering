#ifndef PRE_SOLVER_H
#define PRE_SOLVER_H

void halo_update(int myrank, int Ni, int Nproc, double *Ti_c);
double lambdaharm(double lambdai, double lambdaj);
void calc_coefficients(double *lambdas, double lambda1, double lambda2, double *a, int Ni, int P, 
    int mystart, int myrank);

#endif