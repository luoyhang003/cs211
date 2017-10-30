#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include "lapacke.h"
//#include "blas.h"

void mydgetrf_blockversion(double* a, int* pvt, double* tempv, int N, int block);
void mydtrsm_forward(int N, double* A, int* pvt, double* b, double* x, double *y);
void mydtrsm_back(int N, double* A, int* pvt, double* b, double* x, double* y);
double err(double *a, double *b, int n);

int main(int argc, char** argv) {
        int N = 3000;
        int max_block_size = 100;

        struct timespec start, end;
        double running;
        double gflops;

        double max, min;
        max = N;
        min = 0.0;
        double range = 0.0+(max - min);
        double div = RAND_MAX / range;
        srand(time(NULL));

        double *a;

        for(int i = 10; i < max_block_size; i+=10) {
                printf("Block Size: %d\n", i);

                a = (double*)calloc(sizeof(double), N*N);

                for(int j = 0; j < N*N; j++) {
                        a[j] = min + (rand()/div);
                }

                int *pvt = (int*)calloc(sizeof(int), N);

                for(int j = 0; j < N; j++) {
                        pvt[j] = j;
                }

                double* tempv = (double*)calloc(sizeof(double), N);

                clock_gettime(CLOCK_MONOTONIC, &start);

                mydgetrf_blockversion(a, pvt, tempv, N, i);


                clock_gettime(CLOCK_MONOTONIC, &end);

                running = ((double)end.tv_sec + 1.0e-9*end.tv_nsec) -
                          ((double)start.tv_sec + 1.0e-9*start.tv_nsec);
                gflops = (2.0*powf(N,3.0))/(3.0*running*powf(10.0,9.0));

                printf("Running Time: %.6f seconds.\n", running);
                printf("GFLOPS: %f\n", gflops);


        }

}

void mydgetrf_blockversion(double* a, int* pvt, double* tempv,
                           int N, int block) {
        int maxind, end, temp;
        double max, sum;

        double *mm;

        for(int block_i=0; block_i<N; block_i+=block) {
                end = block_i+block-1;
                for(int i=block_i; i<=end; i++) {
                        maxind = i;
                        max = fabs(a[i*N+i]);
                        for(int t=i+1; t<N; t++) {
                                if(fabs(a[t*N+i])>max) {
                                        maxind = t;
                                        max = fabs(a[t*N+i]);
                                }
                        }
                        if(max==0.0) {
                                printf("LU factorization failed: coefficient matrix is singular\n");
                                return;
                        }
                        else{
                                if(maxind != i) {
                                        temp = pvt[i];
                                        pvt[i] = pvt[maxind];
                                        pvt[maxind] = temp;

                                        for(int k=0; k<N; k++) {
                                                tempv[k] = a[i*N+k];
                                                a[i*N+k] = a[maxind*N+k];
                                                a[maxind*N+k] = tempv[k];
                                        }
                                }
                        }

                        for(int j=i+1; j<N; j++) {
                                a[j*N+i] = a[j*N+i]/a[i*N+i];
                                for(int k=i+1; k<=end; k++) {
                                        a[j*N+k] = a[j*N+k] - a[j*N+i] * a[i*N+k];
                                }
                        }
                }

                mm = (double*)calloc(sizeof(double), block*block);

                int p=0, q=0;

                for(int l = block_i; l<=end; l++) {
                        for(int m=block_i; m<=end; m++) {
                                if(l>m) {
                                        mm[p*block+q] = a[l*N+m]*(-1);
                                }
                                else if(l==m) {
                                        mm[p*block+q] = 1;
                                }
                                else{
                                        mm[p*block+q] = 0;
                                }
                                q++;
                        }
                        p++;
                        q=0;
                }

                p=0;
                q=0;
                for(int j=block_i; j<=end; j++) {
                        for(int k=end+1; k<N; k++) {
                                sum = 0.0;
                                for(int m=block_i; m<=end; m++) {
                                        sum += mm[p*block+q] * a[m*N+k];
                                        q++;
                                }
                                a[j*N+k] = sum;
                                q=0;
                        }
                        p++;
                        q=0;
                }

                for(int j=end+1; j<N; j++) {
                        for(int k=end+1; k<N; k++) {
                                double gg = 0.0;
                                for(int l=block_i; l<=end; l++) {
                                        gg += a[j*N+l] * a[l*N+k];
                                }
                                a[j*N+k] -= gg;
                        }
                }
        }
}

void mydtrsm_forward(int N, double* A, int* pvt, double* b, double* x, double *y) {
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

}


void mydtrsm_back(int N, double* A, int* pvt, double* b, double* x, double* y) {
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
        return error;
        //printf("Error = %f\n",error);
}
