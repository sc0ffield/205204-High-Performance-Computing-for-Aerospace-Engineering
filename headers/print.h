#ifndef PRINT_H
#define PRINT_H

void print_statement(double L, double xchange, double lambda1, 
    double lambda2, double T0, double TL);
void print_discretization(int N, double dx, double* x, double *lambdas, double halo_left, double halo_right);
void print_bc(double TL, double T0);
void print_initialf(double Tini);
void print_iteration();
void print_results(double *T, int N, int iter, double delta);

#endif