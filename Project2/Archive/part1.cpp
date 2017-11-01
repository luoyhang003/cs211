#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "lapacke.h"
#include "blas.h"

// using namespace std;

int mydgetrf(double* m, int N, int* pvt, double* tempv);
//int mydtrsm(int N, double* A, int* pvt, double* b, double* x, double* y);
int transpose(double *a, int N);
double err(double *a, double *b, int n);
int mydtrsm_forward(int N, double* A, int* pvt, double* b, double* x, double *y);
int mydtrsm_back(int N, double* A, int* pvt, double* b, double* x, double* y);


int main(int argc, char* argv[]) {
        if(argc > 1) {
                int N=0;
                N = atoi(argv[1]);
                // printf("N=%d\n", N);

                char TRANS = 'N';
                int INFO = N;
                int LDA = N;
                int LDB = N;
                // N = 3;
                int NRHS = 1;
                int *IPIV = (int *)calloc(sizeof(int), N);

                double* B = (double*)calloc(sizeof(double), N);
                double* B1 = (double*)calloc(sizeof(double), N);



                char SIDE = 'L';
                char UPLO = 'L';
                char DIAG = 'U';
                int M    = 1;
                double a    = 1.0;

                double *A = NULL;
                A = (double*)calloc(sizeof(double), N*N);

                double *A1 = NULL;
                A1 = (double*)calloc(sizeof(double), N*N);

                double *b = NULL;
                b = (double*)calloc(sizeof(double), N);

                double *x = NULL;
                x = (double*)calloc(sizeof(double), N);

                double* tempv = NULL;
                tempv = (double*)calloc(sizeof(double), N);

                double max, min;
                max = N;
                min = 0.0;
                double range = 0.0+(max - min);
                double div = RAND_MAX / range;
                srand(time(NULL));

                for(int i = 0; i < N; i++) {
                        B[i] = min + (rand()/div);
                        B1[i] = B[i];
                }

                for(int index = 0; index < N*N; index++) {
                        A[index] = min + (rand()/div);
                        A1[index] = A[index];
                }

                for(int index = 0; index < N; index++) {
                        b[index] = min + (rand()/div);
                }



                int *pvt = (int*)calloc(sizeof(int), N);

                for(int i = 0; i < N; i++) {
                        pvt[i] = i;
                }

                // for(int i = 0; i < N; i++) {
                //         for(int j = 0; j < N; j++) {
                //                 printf("%f ",A[i*N+j]);
                //         }
                //         printf("\n");
                // }

                struct timespec start, end;
                double running;
                double gflops;

                clock_gettime(CLOCK_MONOTONIC, &start);

                int n = mydgetrf(A, N, pvt, tempv);

                clock_gettime(CLOCK_MONOTONIC, &end);

                running = ((double)end.tv_sec + 1.0e-9*end.tv_nsec) -
                          ((double)start.tv_sec + 1.0e-9*start.tv_nsec);
                gflops = (2.0*powf(N,3.0))/(3.0*running*powf(10.0,9.0));
                printf("MY Running: %.6f seconds.\n", running);
                printf("MY GFLOPS: %f\n", gflops);
                // printf("Error: %d\n", n);

                // for(int i = 0; i < N; i++) {
                //         for(int j = 0; j < N; j++) {
                //                 printf("%f ",A[i*N+j]);
                //         }
                //         printf("\n");
                // }
                //
                // for(int j = 0; j < N; j++) {
                //         printf("%f ",pvt[j]);
                // }

                int trans = transpose(A1, N);

                clock_gettime(CLOCK_MONOTONIC, &start);

                // LU factorization
                LAPACK_dgetrf(&N,&N,A1,&LDA,IPIV,&INFO);

                clock_gettime(CLOCK_MONOTONIC, &end);

                running = ((double)end.tv_sec + 1.0e-9*end.tv_nsec) -
                          ((double)start.tv_sec + 1.0e-9*start.tv_nsec);
                gflops = (2.0*powf(N,3.0))/(3.0*running*powf(10.0,9.0));
                printf("LAPACK Running: %.6f seconds.\n", running);
                printf("LAPACK GFLOPS: %f\n", gflops);

                trans = transpose(A1, N);

                double e = err(A, A1, N);

                printf("err: %f\n", e);

                for(int i = 0; i < N; i++)
                {
                        double tmp = B[IPIV[i]-1];
                        B[IPIV[i]-1] = B[i];
                        B[i] = tmp;
                }

                // forward  L(Ux) = B => y = Ux
                dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A1, &N, B, &N);
                UPLO = 'U';
                DIAG = 'N';
                // backward Ux = y
                dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A1, &N, B, &N);

        }else {
                printf("Please input the size of the MATRIX.\n");
        }
}

int mydgetrf(double* m, int N, int *pvt, double* tempv) {
        double max;
        int maxind, temp;
        for(int i = 0; i < N-1; i++) {
                maxind = i;
                max = fabs(m[i*N+i]);
                for(int t = i+1; t < N; t++) {
                        if(fabs(m[t*N+i]) > max) {
                                maxind = t;
                                max = fabs(m[t*N+i]);
                        }
                }
                if(max == 0) {
                        printf("LUfactoration failed: coefficient matrix is singular.\n");
                        return -1;
                }else {
                        if(maxind != i) {
                                //save pivoting infomation
                                temp = pvt[i];
                                pvt[i] = pvt[maxind];
                                pvt[maxind] = temp;
                                //swap rows
                                for(int j = 0; j < N; j++) {
                                        // double tempv = m[i*N+j];
                                        tempv[j] = m[i*N+j];
                                        m[i*N+j] = m[maxind*N+j];
                                        m[maxind*N+j] = tempv[j];
                                }
                        }
                        //factorization
                        for(int j = i+1; j < N; j++) {
                                m[j*N+i] = m[j*N+i]/m[i*N+i];
                                for(int k = i+1; k < N; k++) {
                                        m[j*N+k] = m[j*N+k] - m[j*N+i]*m[i*N+k];
                                }
                        }
                }

        }
        return 0;

}

int mydtrsm_forward(int N, double* A, int* pvt, double* b, double* x, double *y) {
        //forward substitution.
        double sum = 0.0;
        // double temp;

        // double* y = (double*)calloc(sizeof(double), N);
        y[0] = b[pvt[0]];
        for(int i = 1; i < N; i++) {
                for(int j = 0; j < i-1; j++) {
                        sum += y[j] * A[j*N+j];
                }
                y[i] = b[pvt[i]] - sum;
        }

        return 0;

}


int mydtrsm_back(int N, double* A, int* pvt, double* b, double* x, double* y) {
        //forward substitution.
        double sum = 0.0;

        //back substitution.
        x[N-1] = y[N-1]/A[N*N-1];
        for(int i = N-2; i >= 0; i--) {
                double sum = 0.0;
                for(int j = i+1; j < N; j++) {
                        sum += x[j]*A[i*N+j];
                }
                x[i] = (y[i] - sum)/A[i*N+i];
        }
        return 0;

}

double err(double *a, double *b, int n){
        int i,j;
        double error = 0.0;
        for(i=0; i<n; i++) {
                for(j=0; j<n; j++) {
                        if(error < abs(a[i*n+j]-b[i*n+j]))
                                error = abs(a[i*n+j]-b[i*n+j]);
                }
        }
        //printf("Error = %f\n",error);
}

int transpose(double *a, int N) {
        double temp;
        for(int i=0; i<N; i++) {
                for(int j = i; j < N; j++) {
                        temp = a[i*N+j];
                        a[i*N+j] = a[i*N+i];
                        a[j*N+i] = temp;
                }
        }
}
