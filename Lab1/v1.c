#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#define N 2500
#define t -0.01
#define e 10e-6
void fill_vector(double *vector, double value, int n)
{
    for (int i=0; i<n; ++i) {
        vector[i] = value;
    }
}
int main(int argc,char *argv[]) {
    int size, rank;
    double starttime = 0, endtime = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int *array = malloc(sizeof(int)*size);
    int res=N%size;
    int Nx = sqrt(N);
    int Ny = NNx;
    int coef = Nsize;
    int end = 0;
    array[0] = 0;
    for (int i = 1; i < size; i++){
        array[i] = array[i - 1] + coef + (i - 1 < res);
    }
    if (rank == size-1)
        end = N;
    else end = array[rank+1];
    if (rank < res) coef++;
    int *array_coef = malloc (sizeof(int)*size);
    MPI_Allgather(&coef, 1, MPI_INT, array_coef, 1, MPI_INT, MPI_COMM_WORLD);
    double *A = malloc(sizeof(double) * N * coef);
    for (int i = 0; i < coef; ++i){
        for (int j = 0; j < N; ++j){
            A[i*N+j] = 0;
        }
    }
    int l = 1, p = 1;
    for (int r = 0; r < size; r++) {
        if(rank==r) {
            for (int i = 0; i < coef; i++) {
                for (int j = 0; j < N; ++j) {
                    if (i + array[rank] == j) A[i * N + j] = -4;
                    if (i + array[rank] + 1 == j && j != Nx * l) {
                        A[i * N + j] = 1;
                    }
                    if (i + array[rank] + 1 == j && j == Nx * l)
                        l++;
                    if (i + array[rank] - 1 == j && j != Nx * p - 1) {
                        A[i * N + j] = 1;
                    }
                    if (i + array[rank] - 1 == j && j == Nx * p - 1)
                        p++;
                }
            }
            MPI_Send(&l, 1, MPI_INT, (rank+1)%size, 0, MPI_COMM_WORLD);
            MPI_Send(&p, 1, MPI_INT, (rank+1)%size, 1, MPI_COMM_WORLD);
        }
        if(rank == (r+1)%size){
            MPI_Recv(&l, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&p, 1, MPI_INT, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    for (int k = 1+(Nysize)*rank; k < Ny; ++k){
        for(int i = 0; i < array_coef[rank]; i++) {
            for (int j = Nx*k; j < Nx*(k+1); ++j) {
                if(i + array[rank] + Nx == j) A[i*N+j] = 1;
            }
        }
    }
    for (int k = (Nysize)*rank - 1 + (rank == 0); k < Ny - 1; ++k){
        for(int i = 0; i < array_coef[rank]; i++) {
            for (int j = Nx*k; j < Nx * (k + 1); ++j) {
                if(i + array[rank] - Nx == j) {
                    A[i * N + j] = 1;
                }
            }
        }
    }
    double *b = malloc(sizeof(double) * N);
    for (int i = 0; i < N; ++i){
        switch (i % 5) {
            case 0:
                b[i] = -50;
                break;
            case 1:
                b[i] = -25;
                break;
            case 2:
                b[i] = 0;
                break;
            case 3:
                b[i] = 25;
                break;
            case 4:
                b[i] = 50;
        }
    }
    double *x = malloc(sizeof(double) * N);
    fill_vector(x,0, N);
    double *result = malloc(sizeof(double)*N);
    double module = 0.0, sum_b;
    for (int i = array[rank]; i < end; ++i) {
        module += b[i] * b[i];
    }
    MPI_Allreduce(&module, &sum_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sum_b = sqrt(sum_b);
    int count = 0;
    starttime = MPI_Wtime();
    while (1)
    {
        double module = 0, sum = 0;
        for (int i = array[rank]; i < end; ++i) {
            result[i] = 0;
            for (int j = 0; j < N; ++j) {
                result[i] += A[(i - array[rank]) * N + j] * x[j];
            }
            result[i] = result[i] - b[i];
            module += result[i] * result[i];
            result[i] = x[i] - t * result[i];
        }
        MPI_Allreduce(&module, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sum = sqrt(sum);
        double div = sum  sum_b;
        printf("%f\n", div);
        if (div < e)
            break;
        count++;
        MPI_Allgatherv(result+array[rank], coef, MPI_DOUBLE, x, array_coef, array, MPI_DOUBLE, MPI_COMM_WORLD);
    }
    endtime = MPI_Wtime();
    if (rank == 0) {
        for (int i = 0; i < N; ++i) {
            printf("%f ", x[i]);
        }
        printf("\nCount of operation %d\n", count);
        printf("Time of running in sec %f\n", endtime - starttime);
    }
    free(A);
    free(b);
    free(x);
    free(result);
    MPI_Finalize();
    return 0;
}

