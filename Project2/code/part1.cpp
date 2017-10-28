#include <stdio.h>
#include "stdlib.h"
#include "math.h"
#include <time.h>
#include <string.h>
#include "lapacke.h"
#include "blas.h"

// using namespace std;

int mydgetrf(double* m, int N, int* pvt);
int mydtrsm(int N, double* A, int* pvt, double* b, double* x);

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

                for(int i = 0; i < N; i++) {
                        B[i] = 1;
                }

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

                double max, min;
                max = N;
                min = 0.0;
                double range = 0.0+(max - min);
                double div = RAND_MAX / range;
                srand(time(NULL));

                for(int index = 0; index < N*N; index++) {
                        A[index] = min + (rand()/div);
                        A1[index] = A[index];
                }

                for(int index = 0; index < N; index++) {
                        b[index] = min + (rand()/div);
                }



                int *pvt = new int[N];

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

                clock_gettime(CLOCK_MONOTONIC, &start);

                int n = mydgetrf(A, 3, pvt);

                clock_gettime(CLOCK_MONOTONIC, &end);



                running = 1000000000L * (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec);

                printf("MY Running: %f nonoseconds.\n", running);
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

                clock_gettime(CLOCK_MONOTONIC, &start);

// LU factorization
                LAPACK_dgetrf(&N,&N,A1,&LDA,IPIV,&INFO);

                clock_gettime(CLOCK_MONOTONIC, &end);

                running = 1000000000L * (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec);

                printf("LAPACK Running: %f nonoseconds.\n", running);

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

int mydgetrf(double* m, int N, int *pvt) {
        double max;
        int maxind;
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
                                double temps = pvt[i];
                                pvt[i] = pvt[maxind];
                                pvt[maxind] = temps;
                                //swap rows
                                for(int j = 0; j < N; j++) {
                                        double tempv = m[i*N+j];
                                        m[i*N+j] = m[maxind*N+j];
                                        m[maxind*N+j] = tempv;
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

        // for(int j = 0; j < N; j++) {
        //         printf("%f ",pvt[j]);
        // }


        // for(int i = 0; i < N; i++) {
        //         for(int j = 0; j < N; j++) {
        //                 printf("%f ", m[i*N+j]);
        //         }
        //         printf("\n");
        // }

        return 0;

}

int mydtrsm(int N, double* A, int* pvt, double* b, double* x) {
        //forward substitution.
        double* y = (double*)calloc(sizeof(double), N);
        y[1] = b[pvt[1]];
        for(int i = 1; i < N; i++) {
                double sum = 0.0;
                for(int j = 0; j < i-1; j++) {
                        sum += y[j] * A[j*N+j];
                }
                y[i] = b[pvt[i]] - sum;
        }
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
