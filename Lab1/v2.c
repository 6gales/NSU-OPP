//#include<stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <mpi.h>
//#include <time.h>
//#include <memory.h>
//#define N 2500
//#define t -0.01
//#define e 1.0e-6
//void fill_vector(double *vector, double value, int n) {
//    for (int i = 0; i < n; ++i)
//        vector[i] = value;
//}
//int main(int argc,char *argv[]) {
//    int size, rank;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
//	int Nx = sqrt(N);
//	int Ny = N / Nx;
//	double starttime = 0, endtime = 0;
//    int *array = malloc(sizeof(int)*size);
////    if (rank == 0)
////        freopen("output0.txt", "w", stdout);
////    if (rank == 1)
////        freopen("output1.txt", "w", stdout);
////    if (rank == 2)
////        freopen("output2.txt", "w", stdout);
////    if (rank == 3)
////        freopen("output3.txt", "w", stdout);
////    if (rank == 4)
////        freopen("output4.txt", "w", stdout);
////    if (rank == 5)
////        freopen("output5.txt", "w", stdout);
////    if (rank == 6)
////        freopen("output6.txt", "w", stdout);
////    if (rank == 7)
////        freopen("output7.txt", "w", stdout);
//    int res=N%size;
//    int coef = N/size;
//    array[0] = 0;
//    for (int i = 1; i < size; i++){
//        array[i] = array[i - 1] + coef + (i - 1 < res);
//    }
//    if (rank == 0) {
//        for (int i = 0; i < size; i++) {
//            printf("%d ", array[i]);
//        }
//    }
//    if (rank < res) coef++;
//    //printf ("#%d %d\n", rank, coef);
//    int *array_coef = malloc (sizeof(int)*size);
//    MPI_Allgather(&coef, 1, MPI_INT, array_coef, 1, MPI_INT, MPI_COMM_WORLD);
//        if (rank == 0) {
////        for (int i = 0; i < size; i++) {
////            printf("%d ", array_coef[i]);
////        }
//    }
//    double *A = malloc(sizeof(double) * N * coef);
//    for (int i = 0; i < coef; ++i){
//        for (int j = 0; j < N; ++j){
//            A[i*N+j] = 0;
//        }
//    }
//   // printf("#%d aaaa\n", rank);
//    int l = 1, p = 1;
//    for (int r = 0; r < size; r++) {
//        if(rank==r) {
//            for (int i = 0; i < coef; i++) {
//                for (int j = 0; j < N; ++j) {
//                    if (i + array[rank] == j) A[i * N + j] = -4;
//                    if (i + array[rank] + 1 == j && j != Nx * l) {
//                        A[i * N + j] = 1;
//                    }
//                    if (i + array[rank] + 1 == j && j == Nx * l)
//                        l++;
//                    if (i + array[rank] - 1 == j && j != Nx * p - 1) {
//                        A[i * N + j] = 1;
//                    }
//                    if (i + array[rank] - 1 == j && j == Nx * p - 1)
//                        p++;
//                }
//            }
//            if (rank != size - 1) {
//                printf("#%d here sending\n", rank);
//                MPI_Send(&l, 1, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD);
//                MPI_Send(&p, 1, MPI_INT, (rank + 1) % size, 1, MPI_COMM_WORLD);
//            }
//        }
//        if(rank == (r+1)%size && rank != 0){
//            printf("#%d here recv\n", rank);
//            MPI_Recv(&l, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            MPI_Recv(&p, 1, MPI_INT, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        printf("#%d here waiting on iter %d\n", rank, r);
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
////    printf("#%d here passed\n", rank);
//    for (int k = 1+(Ny/size)*rank; k < Ny; ++k){
//        for(int i = 0; i < array_coef[rank]; i++) {
//            for (int j = Nx*k; j < Nx*(k+1); ++j) {
//                if(i + array[rank] + Nx == j) A[i*N+j] = 1;
//            }
//        }
//    }
//////printf("#%d hhhh\n", rank);
//    for (int k = (Ny/size)*rank - 1 + (rank == 0); k < Ny - 1; ++k){
//        for(int i = 0; i < array_coef[rank]; i++) {
//            for (int j = Nx*k; j < Nx * (k + 1); ++j) {
//                if(i + array[rank] - Nx == j) {
//                    A[i * N + j] = 1;
//                }
//            }
//        }
//    }
//    //for (int r = 0; r < size; ++r) {
//    //    if (r == rank) {
//    //        for (int i = 0; i < coef; ++i) {
//    //            for (int j = 0; j < N; ++j) {
//    //                printf("#%d %f ", rank, A[i * N + j]);
//    //            }
//    //            printf("\n");
//    //        }
//    //    }
//    //    MPI_Barrier(MPI_COMM_WORLD);
//    //}
//    double *b = malloc(sizeof(double) * coef);
//
//            for (int i = 0; i < coef; ++i) {
//                //printf("%d %d\n", array[rank], i);
//                switch ((i + array[rank]) % 5) {
//                    case 0:
//                        b[i] = -50;
//                        break;
//                    case 1:
//                        b[i] = -25;
//                        break;
//                    case 2:
//                        b[i] = 0;
//                        break;
//                    case 3:
//                        b[i] = 25;
//                        break;
//                    case 4:
//                        b[i] = 50;
//                }
//            }
////    for (int r = 0; r < size; ++r){
////        if (r == rank) {
////            for (int i = 0; i < coef; ++i) {
////                printf("%f ", b[i]);
////            }
////        }
////        MPI_Barrier(MPI_COMM_WORLD);
////    }
//    double *x = malloc(sizeof(double) * coef);
//    fill_vector(x,0, coef);
////        //if (rank == 0) {
////        for (int i = 0; i < coef; i++) {
////            printf("%f ", x[i]);
////        //}
////    }
//    double *result = malloc(sizeof(double)*coef);
//    double module = 0.0, sum_b;
//    for (int i = 0; i < coef; ++i) {
//        module += b[i] * b[i];
//    }
//    //printf("#%d %f ", rank, module);
//    MPI_Allreduce(&module, &sum_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    sum_b = sqrt(sum_b);
////    printf("#%d %f\n", rank, sum_b);
//    double *tmp = malloc(sizeof(double)*array_coef[0]);
//    int count = 0;
//    starttime = MPI_Wtime();
//	
//	while (1)
//    {
// //       printf("#%d before\n", rank);
//        for (int i = 0; i < coef; ++i)
//        {
//            result[i] = 0;
//        }
//
//        for (int k = 0; k < size; ++k)
//        {
// //           printf("#%d is on memory on iteration %d\n", rank, k);
//            memcpy(tmp, x, coef*sizeof(double));
////            printf("#%d before on iteration %d\n", rank, k);
//            MPI_Bcast(tmp, array_coef[0], MPI_DOUBLE, k, MPI_COMM_WORLD);
//  //          printf("#%d on iteration %d\n", rank, k);
//            for (int i = 0; i < coef; ++i) {
//                for (int j = 0; j < array_coef[k]; ++j) {
//                    result[i] += A[i * N + array[k] + j] * tmp[j];
//    //                if (rank == 2) printf("______________i: %d, j: %d\n", i, j);
//                }
//            }
//      //      printf("#%d after for on iteration %d\n", rank, k);
//        }
//      //  printf("#%d middle\n", rank);
//        double sum = 0, module = 0;
//        for (int i = 0; i < coef; ++i){
//            result[i] = result[i] - b[i];
//            module += result[i] * result[i];
//            result[i] = x[i] - t * result[i];
//        }
//        MPI_Allreduce(&module, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        sum = sqrt(sum);
//        double div = sum / sum_b;
//   //     printf("%f\n", div);
//        if (div < e)
//            break;
//        //printf("%f\n", div);
//        count ++;
//        double *swap = result;
//        result = x;
//        x = swap;
//    }
//    double *ptr = NULL;
//    if (rank == 0)
//        ptr = malloc(sizeof(double)*N);
//   // MPI_Gather(result, coef, MPI_DOUBLE, ptr, coef, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    MPI_Gatherv(result, coef, MPI_DOUBLE, ptr, array_coef, array, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    endtime = MPI_Wtime();
//    if (rank == 0) {
//        for (int i = 0; i < N; ++i) {
//           printf("%f ", ptr[i]);
//        }
//        printf("Count of operation %d\n", count);
//        printf("Time of running %f\n", endtime - starttime);
//        free(ptr);
//    }
//    free(A);
//    free(b);
//    free(x);
//    free(result);
//    free(tmp);
//    free(array);
//    free(array_coef);
//    MPI_Finalize();
//    return 0;
//}