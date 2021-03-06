/*
 * Course: CS 211 - High Performance Computing
 * Project 1: Performance Optimization via Register and Cache Reuse
 *
 * Author: Yuanhang Luo
 * Create Date: Oct 9th, 2017
 * Revise Date: Oct 10th, 2017
 *
 */


//just for debug
#define DEBUG 0


#include "stdio.h"
#include "string.h"
#include "time.h"
#include "stdlib.h"
#include "float.h"
#include "math.h"
#include "time.h"


int main(int argc, char *argv[]) {




        if(argc > 1) {
                if(!strcmp(argv[1], "help")) {
                        printf("=========Instructions=============\n");
                }
                else if(!strcmp(argv[1], "dgemm0")) {
                        printf("=========dgemm0=================\n");

                        int N, loop;

                        if(argv[2] == NULL) {
                                for(loop = 0; loop < 6; loop++) {
                                        N = pow(2, 6+loop);
                                        double *a = NULL;
                                        double *b = NULL;
                                        double *c = NULL;

                                        a = (double *)calloc(sizeof(double), N*N);
                                        b = (double *)calloc(sizeof(double), N*N);
                                        c = (double *)calloc(sizeof(double), N*N);

                                        int index;

                                        double max, min;
                                        max = N*N*N;
                                        min = 0.0-N*N*N;
                                        double range = 0.0+(max - min);
                                        double div = RAND_MAX / range;
                                        srand(time(NULL));

                                        for(index = 0; index < N*N; index++) {
                                                a[index] = min + (rand()/div);
                                                b[index] = min + (rand()/div);
                                        }
                                        int i, j, k;

                                        clock_t begin;
                                        clock_t end;
                                        double running;

                                        begin = clock();
                                        //dgemm0: simple ijk version triple loop algorithm
                                        for(i = 0; i < N; i++) {
                                                for(j = 0; j < N; j++) {
                                                        for(k = 0; k < N; k++) {
                                                                c[i*N+j] += a[i*N+k] * b[k*N+j];
                                                        }
                                                }
                                        }
                                        end = clock();
                                        running = (double)(end - begin) / CLOCKS_PER_SEC;

                                        printf("The running time of dgemm0 with n=%d is %f s\n", N, running);

                                        free(a);
                                        free(b);
                                        free(c);
                                }
                        }else{
                                N = atoi(argv[2]);
                                double *a = NULL;
                                double *b = NULL;
                                double *c = NULL;

                                a = (double *)calloc(sizeof(double), N*N);
                                b = (double *)calloc(sizeof(double), N*N);
                                c = (double *)calloc(sizeof(double), N*N);

                                int index;

                                double max, min;
                                max = N*N*N;
                                min = 0.0-N*N*N;
                                double range = 0.0+(max - min);
                                double div = RAND_MAX / range;
                                srand(time(NULL));

                                for(index = 0; index < N*N; index++) {
                                        a[index] = min + (rand()/div);
                                        b[index] = min + (rand()/div);
                                }

                                int i, j, k;

                                clock_t begin;
                                clock_t end;
                                double running;

                                begin = clock();
                                //dgemm0: simple ijk version triple loop algorithm
                                for(i = 0; i < N; i++) {
                                        for(j = 0; j < N; j++) {
                                                for(k = 0; k < N; k++) {
                                                        c[i*N+j] += a[i*N+k] * b[k*N+j];
                                                }
                                        }
                                }
                                end = clock();

                                running = (double)(end - begin) / CLOCKS_PER_SEC;
                                printf("The running time of dgemm0 with n=%d is %f s\n", N, running);

                                free(a);
                                free(b);
                                free(c);
                        }



                }


                else if(!strcmp(argv[1], "dgemm1")) {
                        printf("=========dgemm1=================\n");

                        int N, loop;

                        if(argv[2] == NULL) {
                                for(loop = 0; loop < 6; loop++) {
                                        N = pow(2, 6+loop);
                                        double *a = NULL;
                                        double *b = NULL;
                                        double *c = NULL;

                                        a = (double *)calloc(sizeof(double), N*N);
                                        b = (double *)calloc(sizeof(double), N*N);
                                        c = (double *)calloc(sizeof(double), N*N);

                                        int index;

                                        double max, min;
                                        max = N*N*N;
                                        min = 0.0-N*N*N;
                                        double range = 0.0+(max - min);
                                        double div = RAND_MAX / range;
                                        srand(time(NULL));

                                        for(index = 0; index < N*N; index++) {
                                                a[index] = min + (rand()/div);
                                                b[index] = min + (rand()/div);
                                        }
                                        int i, j, k;

                                        clock_t begin;
                                        clock_t end;
                                        double running;

                                        begin = clock();

                                        //dgemm1: simple ijk version triple loop algorithm with register reuse
                                        for(i = 0; i < N; i++) {
                                                for(j = 0; j < N; j++) {
                                                        register double r = c[i*N+j];
                                                        for(k = 0; k < N; k++) {
                                                                r += a[i*N+k] * b[k*N+j];
                                                        }
                                                        c[i*N+j] = r;
                                                }
                                        }

                                        end = clock();
                                        running = (double)(end - begin) / CLOCKS_PER_SEC;

                                        printf("The running time of dgemm1 with n=%d is %f s\n", N, running);

                                        free(a);
                                        free(b);
                                        free(c);
                                }
                        }else{
                                N = atoi(argv[2]);
                                double *a = NULL;
                                double *b = NULL;
                                double *c = NULL;

                                a = (double *)calloc(sizeof(double), N*N);
                                b = (double *)calloc(sizeof(double), N*N);
                                c = (double *)calloc(sizeof(double), N*N);

                                int index;

                                double max, min;
                                max = N*N*N;
                                min = 0.0-N*N*N;
                                double range = 0.0+(max - min);
                                double div = RAND_MAX / range;
                                srand(time(NULL));

                                for(index = 0; index < N*N; index++) {
                                        a[index] = min + (rand()/div);
                                        b[index] = min + (rand()/div);
                                }

                                int i, j, k;

                                clock_t begin;
                                clock_t end;
                                double running;

                                begin = clock();

                                //dgemm1: simple ijk version triple loop algorithm with register reuse
                                for(i = 0; i < N; i++) {
                                        for(j = 0; j < N; j++) {
                                                register double r = c[i*N+j];
                                                for(k = 0; k < N; k++) {
                                                        r += a[i*N+k] * b[k*N+j];
                                                }
                                                c[i*N+j] = r;
                                        }
                                }

                                end = clock();

                                running = (double)(end - begin) / CLOCKS_PER_SEC;
                                printf("The running time of dgemm1 with n=%d is %f s\n", N, running);

                                free(a);
                                free(b);
                                free(c);
                        }



                }


                else if(!strcmp(argv[1], "dgemm2")) {
                        printf("=========dgemm2=================\n");

                        int N, loop;

                        if(argv[2] == NULL) {
                                for(loop = 0; loop < 6; loop++) {
                                        N = pow(2, 6+loop);
                                        double *a = NULL;
                                        double *b = NULL;
                                        double *c = NULL;

                                        a = (double *)calloc(sizeof(double), N*N);
                                        b = (double *)calloc(sizeof(double), N*N);
                                        c = (double *)calloc(sizeof(double), N*N);

                                        int index;

                                        double max, min;
                                        max = N*N*N;
                                        min = 0.0-N*N*N;
                                        double range = 0.0+(max - min);
                                        double div = RAND_MAX / range;
                                        srand(time(NULL));

                                        for(index = 0; index < N*N; index++) {
                                                a[index] = min + (rand()/div);
                                                b[index] = min + (rand()/div);
                                        }
                                        int i, j, k;

                                        clock_t begin;
                                        clock_t end;
                                        double running;

                                        begin = clock();

                                        //dgemm2: more aggresive register reuse
                                        for(i = 0; i < N; i += 2) {
                                                for(j = 0; j < N; j += 2) {
                                                        for(k = 0; k < N; k += 2) {
                                                                c[i*N + j] = a[i*N + k] * b[k*N + j] + a[i*N + k+1] * b[(k+1)*N + j] + c[i*N + j];
                                                                c[(i+1)*N + j] = a[(i+1)*N + k] * b[k*N + j] + a[(i+1)*N + k+1] * b[(k+1)*N + j] + c[(i+1)*N + j];
                                                                c[i*N + (j+1)] = a[i*N + k] * b[k*N + (j+1)] + a[i*N + k+1] * b[(k+1)*N + (j+1)] + c[i*N + (j+1)];
                                                                c[(i+1)*N + (j+1)] = a[(i+1)*N + k] * b[k*N + (j+1)] + a[(i+1)*N + k+1] * b[(k+1)*N + (j+1)] + c[(i+1)*N + (j+1)];
                                                        }
                                                }

                                        }

                                        end = clock();
                                        running = (double)(end - begin) / CLOCKS_PER_SEC;

                                        printf("The running time of dgemm2 with n=%d is %f s\n", N, running);

                                        free(a);
                                        free(b);
                                        free(c);
                                }
                        }else{
                                N = atoi(argv[2]);
                                double *a = NULL;
                                double *b = NULL;
                                double *c = NULL;

                                a = (double *)calloc(sizeof(double), N*N);
                                b = (double *)calloc(sizeof(double), N*N);
                                c = (double *)calloc(sizeof(double), N*N);

                                int index;

                                double max, min;
                                max = N*N*N;
                                min = 0.0-N*N*N;
                                double range = 0.0+(max - min);
                                double div = RAND_MAX / range;
                                srand(time(NULL));

                                for(index = 0; index < N*N; index++) {
                                        a[index] = min + (rand()/div);
                                        b[index] = min + (rand()/div);
                                }

                                int i, j, k;

                                clock_t begin;
                                clock_t end;
                                double running;

                                begin = clock();

                                //dgemm2: more aggresive register reuse
                                for(i = 0; i < N; i += 2) {
                                        for(j = 0; j < N; j += 2) {
                                                for(k = 0; k < N; k += 2) {
                                                        c[i*N + j] = a[i*N + k] * b[k*N + j] + a[i*N + k+1] * b[(k+1)*N + j] + c[i*N + j];
                                                        c[(i+1)*N + j] = a[(i+1)*N + k] * b[k*N + j] + a[(i+1)*N + k+1] * b[(k+1)*N + j] + c[(i+1)*N + j];
                                                        c[i*N + (j+1)] = a[i*N + k] * b[k*N + (j+1)] + a[i*N + k+1] * b[(k+1)*N + (j+1)] + c[i*N + (j+1)];
                                                        c[(i+1)*N + (j+1)] = a[(i+1)*N + k] * b[k*N + (j+1)] + a[(i+1)*N + k+1] * b[(k+1)*N + (j+1)] + c[(i+1)*N + (j+1)];
                                                }
                                        }

                                }

                                end = clock();

                                running = (double)(end - begin) / CLOCKS_PER_SEC;
                                printf("The running time of dgemm2 with n=%d is %f s\n", N, running);

                                free(a);
                                free(b);
                                free(c);
                        }



                }

                else if(!strcmp(argv[1], "dgemm3")) {
                        printf("=========dgemm3=================\n");

                        int N, loop;

                        if(argv[2] == NULL) {
                                for(loop = 0; loop < 6; loop++) {
                                        N = pow(2, 6+loop);
                                        double *a = NULL;
                                        double *b = NULL;
                                        double *c = NULL;

                                        a = (double *)calloc(sizeof(double), N*N);
                                        b = (double *)calloc(sizeof(double), N*N);
                                        c = (double *)calloc(sizeof(double), N*N);

                                        int index;

                                        double max, min;
                                        max = N*N*N;
                                        min = 0.0-N*N*N;
                                        double range = 0.0+(max - min);
                                        double div = RAND_MAX / range;
                                        srand(time(NULL));

                                        for(index = 0; index < N*N; index++) {
                                                a[index] = min + (rand()/div);
                                                b[index] = min + (rand()/div);
                                        }
                                        int i, j, k;

                                        clock_t begin;
                                        clock_t end;
                                        double running;

                                        begin = clock();

                                        //dgemm3
                                        for(i = 0; i < N; i += 2) {
                                                for(j = 0; j < N; j += 2)  {
                                                        register int t = i*N + j;
                                                        register int tt = t + N;
                                                        register double c00 = c[t];
                                                        register double c01 = c[t+1];
                                                        register double c10 = c[tt];
                                                        register double c11 = c[tt+1];

                                                        for(k = 0; k < N; k += 2) {
                                                                /* 2 by 2 mini matrix multiplication using registers*/
                                                                register int ta = i*N + k;
                                                                register int tta = ta + N;
                                                                register int tb = k*N + j;
                                                                register int ttb = tb + N;
                                                                register double a00 = a[ta];
                                                                register double a01 = a[ta+1];
                                                                register double a10 = a[tta];
                                                                register double a11 = a[tta+1];
                                                                register double b00 = b[tb];
                                                                register double b01 = b[tb+1];
                                                                register double b10 = b[ttb];
                                                                register double b11 = b[ttb+1];
                                                                c00 += a00*b00 + a01*b10;
                                                                c01 += a00*b01 + a01*b11;
                                                                c10 += a10*b00 + a11*b10;
                                                                c11 += a10*b01 + a11*b11;
                                                        }

                                                        c[t] = c00;
                                                        c[t+1] = c01;
                                                        c[tt] = c10;
                                                        c[tt+1] = c11;
                                                }

                                        }

                                        end = clock();
                                        running = (double)(end - begin) / CLOCKS_PER_SEC;

                                        printf("The running time of dgemm3 with n=%d is %f s\n", N, running);

                                        free(a);
                                        free(b);
                                        free(c);
                                }
                        }else{
                                N = atoi(argv[2]);
                                double *a = NULL;
                                double *b = NULL;
                                double *c = NULL;

                                a = (double *)calloc(sizeof(double), N*N);
                                b = (double *)calloc(sizeof(double), N*N);
                                c = (double *)calloc(sizeof(double), N*N);

                                int index;

                                double max, min;
                                max = N*N*N;
                                min = 0.0-N*N*N;
                                double range = 0.0+(max - min);
                                double div = RAND_MAX / range;
                                srand(time(NULL));

                                for(index = 0; index < N*N; index++) {
                                        a[index] = min + (rand()/div);
                                        b[index] = min + (rand()/div);
                                }

                                int i, j, k;

                                clock_t begin;
                                clock_t end;
                                double running;

                                begin = clock();

                                //dgemm3
                                for(i = 0; i < N; i += 2) {
                                        for(j = 0; j < N; j += 2)  {
                                                register int t = i*N + j;
                                                register int tt = t + N;
                                                register double c00 = c[t];
                                                register double c01 = c[t+1];
                                                register double c10 = c[tt];
                                                register double c11 = c[tt+1];

                                                for(k = 0; k < N; k += 2) {
                                                        /* 2 by 2 mini matrix multiplication using registers*/
                                                        register int ta = i*N + k;
                                                        register int tta = ta + N;
                                                        register int tb = k*N + j;
                                                        register int ttb = tb + N;
                                                        register double a00 = a[ta];
                                                        register double a01 = a[ta+1];
                                                        register double a10 = a[tta];
                                                        register double a11 = a[tta+1];
                                                        register double b00 = b[tb];
                                                        register double b01 = b[tb+1];
                                                        register double b10 = b[ttb];
                                                        register double b11 = b[ttb+1];
                                                        c00 += a00*b00 + a01*b10;
                                                        c01 += a00*b01 + a01*b11;
                                                        c10 += a10*b00 + a11*b10;
                                                        c11 += a10*b01 + a11*b11;
                                                }

                                                c[t] = c00;
                                                c[t+1] = c01;
                                                c[tt] = c10;
                                                c[tt+1] = c11;
                                        }

                                }
                                end = clock();

                                running = (double)(end - begin) / CLOCKS_PER_SEC;
                                printf("The running time of dgemm3 with n=%d is %f s\n", N, running);

                                free(a);
                                free(b);
                                free(c);
                        }



                }



                // else if(!strcmp(argv[1], "part1")) {
                //         printf("=========PART ONE=================\n");
                //
                //         int N, loop;
                //
                //         if(argv[2] == NULL) {
                //                 for(loop = 0; loop < 6; loop++) {
                //                         N = pow(2, 6+loop);
                //                         double *a = NULL;
                //                         double *b = NULL;
                //                         double *c = NULL;
                //
                //                         a = (double *)calloc(sizeof(double), N*N);
                //                         b = (double *)calloc(sizeof(double), N*N);
                //                         c = (double *)calloc(sizeof(double), N*N);
                //
                //                         int index;
                //
                //                         double max, min;
                //                         max = N*N*N;
                //                         min = 0.0-N*N*N;
                //                         double range = 0.0+(max - min);
                //                         double div = RAND_MAX / range;
                //                         srand(time(NULL));
                //
                //                         for(index = 0; index < N*N; index++) {
                //                                 a[index] = min + (rand()/div);
                //                                 b[index] = min + (rand()/div);
                //                         }
                //
                //                         if(DEBUG) {
                //                                 printf("the range of the matrix is from %f to %f.\n", min, max);
                //                                 printf("Range: %f, Div: %f\n", range, div);
                //                                 printf("Print first N elements of generated random matrices.\n");
                //                                 for(index = 0; index < N; index++) {
                //                                         printf("%f\t", a[index]);
                //                                         printf("%f\t", b[index]);
                //                                         printf("%f\t", c[index]);
                //                                 }
                //                                 printf("\n");
                //
                //                         }
                //
                //                         int i, j, k;
                //
                //                         clock_t begin;
                //                         clock_t end;
                //                         double running;
                //
                //                         begin = clock();
                //                         //dgemm0: simple ijk version triple loop algorithm
                //                         for(i = 0; i < N; i++) {
                //                                 for(j = 0; j < N; j++) {
                //                                         for(k = 0; k < N; k++) {
                //                                                 c[i*N+j] += a[i*N+k] * b[k*N+j];
                //                                         }
                //                                 }
                //                         }
                //                         end = clock();
                //                         running = (double)(end - begin) / CLOCKS_PER_SEC;
                //
                //                         printf("The running time of dgemm0 with n=%d is %f s\n", N, running);
                //
                //                         free(a);
                //                         free(b);
                //                         free(c);
                //
                //                         a = (double *)calloc(sizeof(double), N*N);
                //                         b = (double *)calloc(sizeof(double), N*N);
                //                         c = (double *)calloc(sizeof(double), N*N);
                //                         for(index = 0; index < N*N; index++) {
                //                                 a[index] = min + (rand()/div);
                //                                 b[index] = min + (rand()/div);
                //                         }
                //
                //                         if(DEBUG) {
                //                                 printf("(dgemm0) The first N elements in calculated matrix C:\n");
                //                                 for(index = 0; index < N; index++) {
                //                                         printf("%f\t", c[index]);
                //                                 }
                //                                 printf("\n");
                //                         }
                //
                //
                //                         begin = clock();
                //                         //dgemm1: simple ijk version triple loop algorithm with register reuse
                //                         for(i = 0; i < N; i++) {
                //                                 for(j = 0; j < N; j++) {
                //                                         register double r = c[i*N+j];
                //                                         for(k = 0; k < N; k++) {
                //                                                 r += a[i*N+k] * b[k*N+j];
                //                                         }
                //                                         c[i*N+j] = r;
                //                                 }
                //                         }
                //                         end = clock();
                //                         running = (double)(end - begin) / CLOCKS_PER_SEC;
                //                         printf("The running time of dgemm1 with n=%d is %f s, Gflops is %f \n", N, running, N*N*N*2.0/running/1000000000.0);
                //
                //                         if(DEBUG) {
                //                                 printf("(dgemm1) The first N elements in calculated matrix C:\n");
                //                                 for(index = 0; index < N; index++) {
                //                                         printf("%f\t", c[index]);
                //                                 }
                //                                 printf("\n");
                //                         }
                //                 }
                //         }else {
                //                 N = atoi(argv[2]);
                //                 double *a = NULL;
                //                 double *b = NULL;
                //                 double *c = NULL;
                //
                //                 a = (double *)calloc(sizeof(double), N*N);
                //                 b = (double *)calloc(sizeof(double), N*N);
                //                 c = (double *)calloc(sizeof(double), N*N);
                //
                //                 int index;
                //
                //                 double max, min;
                //                 max = N*N*N;
                //                 min = 0.0-N*N*N;
                //                 double range = 0.0+(max - min);
                //                 double div = RAND_MAX / range;
                //                 srand(time(NULL));
                //
                //                 for(index = 0; index < N*N; index++) {
                //                         a[index] = min + (rand()/div);
                //                         b[index] = min + (rand()/div);
                //                 }
                //
                //                 int i, j, k;
                //
                //                 clock_t begin;
                //                 clock_t end;
                //                 double running;
                //
                //                 begin = clock();
                //                 //dgemm0: simple ijk version triple loop algorithm
                //                 for(i = 0; i < N; i++) {
                //                         for(j = 0; j < N; j++) {
                //                                 for(k = 0; k < N; k++) {
                //                                         c[i*N+j] += a[i*N+k] * b[k*N+j];
                //                                 }
                //                         }
                //                 }
                //                 end = clock();
                //
                //                 running = (double)(end - begin) / CLOCKS_PER_SEC;
                //                 printf("The running time of dgemm0 with n=%d is %f s\n", N, running);
                //
                //                 free(a);
                //                 free(b);
                //                 free(c);
                //
                //                 a = (double *)calloc(sizeof(double), N*N);
                //                 b = (double *)calloc(sizeof(double), N*N);
                //                 c = (double *)calloc(sizeof(double), N*N);
                //
                //                 for(index = 0; index < N*N; index++) {
                //                         a[index] = min + (rand()/div);
                //                         b[index] = min + (rand()/div);
                //                 }
                //
                //                 begin = clock();
                //                 //dgemm1: simple ijk version triple loop algorithm with register reuse
                //                 for(i = 0; i < N; i++) {
                //                         for(j = 0; j < N; j++) {
                //                                 register double r = c[i*N+j];
                //                                 for(k = 0; k < N; k++) {
                //                                         r += a[i*N+k] * b[k*N+j];
                //                                 }
                //                                 c[i*N+j] = r;
                //                         }
                //                 }
                //                 end = clock();
                //
                //                 running = (double)(end - begin) / CLOCKS_PER_SEC;
                //                 printf("The running time of dgemm1 with n=%d is %f s\n", N, running);
                //
                //         }
                // }
                // else if(!strcmp(argv[1], "part2")) {
                //         printf("=========PART TWO===============\n");
                //
                //         int N, loop;
                //
                //         double *a = NULL;
                //         double *b = NULL;
                //         double *c = NULL;
                //
                //         for(loop = 0; loop < 6; loop++) {
                //                 N = pow(2, loop+6);
                //
                //                 a = (double*)calloc(sizeof(double), N*N);
                //                 b = (double*)calloc(sizeof(double), N*N);
                //                 c = (double*)calloc(sizeof(double), N*N);
                //
                //                 int index;
                //
                //                 double max, min;
                //                 max = N*N*N;
                //                 min = 0.0-N*N*N;
                //                 double range = 0.0+(max - min);
                //                 double div = RAND_MAX / range;
                //                 srand(time(NULL));
                //
                //                 for(index = 0; index < N*N; index++) {
                //                         a[index] = min + (rand()/div);
                //                         b[index] = min + (rand()/div);
                //                 }
                //
                //
                //                 //dgemm2: more aggresive register reuse
                //                 int i, j, k;
                //
                //                 clock_t begin;
                //                 clock_t end;
                //                 double running;
                //
                //                 begin = clock();
                //
                //                 for(i = 0; i < N; i += 2) {
                //                         for(j = 0; j < N; j += 2) {
                //                                 for(k = 0; k < N; k += 2) {
                //                                         c[i*N + j] = a[i*N + k] * b[k*N + j] + a[i*N + k+1] * b[(k+1)*N + j] + c[i*N + j];
                //                                         c[(i+1)*N + j] = a[(i+1)*N + k] * b[k*N + j] + a[(i+1)*N + k+1] * b[(k+1)*N + j] + c[(i+1)*N + j];
                //                                         c[i*N + (j+1)] = a[i*N + k] * b[k*N + (j+1)] + a[i*N + k+1] * b[(k+1)*N + (j+1)] + c[i*N + (j+1)];
                //                                         c[(i+1)*N + (j+1)] = a[(i+1)*N + k] * b[k*N + (j+1)] + a[(i+1)*N + k+1] * b[(k+1)*N + (j+1)] + c[(i+1)*N + (j+1)];
                //                                 }
                //                         }
                //
                //                 }
                //
                //                 end = clock();
                //                 running = (double)(end - begin) / CLOCKS_PER_SEC;
                //                 printf("The running time of dgemm2 with n=%d is %f s, \n", N, running);
                //
                //
                //         }
                //
                //
                // }
                // else if(!strcmp(argv[1], "part3")) {
                //         printf("=========PART THREE===============\n");
                // }
                else{
                        printf("Please input the correct parameters or use `./proj1 help` to see instructions.\n");
                }
        }
        else{
                printf("Please input the parameters of the program or use `./proj1 help` to see instructions.\n");
        }
}
