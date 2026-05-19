#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

void worksplit(int *mystart, int *myend, int myrank, int P, int N);
void generate_mesh(int P, int Ni, double *x, double *xi, double dx, int mystart, double xchange,
    double *lambdasi, double *lambdas, double lambda1, double lambda2);
void allgathering(int send_count, int P, double *X, double *sendbuf);

#endif