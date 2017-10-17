/*
 * Course: CS 211 - High Performance Computing
 * Project 1: Performance Optimization via Register and Cache Reuse
 *
 * Author: Yuanhang Luo
 * Create Date: Oct 14th, 2017
 * Revise Date: Oct 14th, 2017
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
    int n = atoi(argv[1]);

    double *a = NULL;
    double *b = NULL;
    double *c0 = NULL;
    double *c1 = NULL;
    double *c2 = NULL;
    double *c3 = NULL;


    a = calloc(sizeof(double), n*n);
    b = calloc(sizeof(double), n*n);
    c0 = calloc(sizeof(double), n*n);
    c1 = calloc(sizeof(double), n*n);
    c2 = calloc(sizeof(double), n*n);
    c3 = calloc(sizeof(double), n*n);

    int index;

    double max, min;
    max = n*n;
    min = 0.0-n*n;
    double range = 0.0+(max - min);
    double div = RAND_MAX / range;
    srand(time(NULL));

    for(index = 0; index < n*n; index++) {
            a[index] = min + (rand()/div);
            b[index] = min + (rand()/div);
    }


    int i, j, k;

    clock_t begin;
    clock_t end;
    double running;


    begin = clock();
    //dgemm0: simple ijk version triple loop algorithm
    for(i = 0; i < n; i++) {
            for(j = 0; j < n; j++) {
                    for(k = 0; k < n; k++) {
                            c0[i*n+j] += a[i*n+k] * b[k*n+j];
                    }
            }
    }
    end = clock();
    running = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("The running time of dgemm0 with n=%d is %fs\n", n, running);


    begin = clock();

    //dgemm1: simple ijk version triple loop algorithm with register reuse
    for(i = 0; i < n; i++) {
            for(j = 0; j < n; j++) {
                    register double r = c1[i*n+j];
                    for(k = 0; k < n; k++) {
                            r += a[i*n+k] * b[k*n+j];
                    }
                    c1[i*n+j] = r;
            }
    }

    end = clock();
    running = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("The running time of dgemm1 with n=%d is %fs\n", n, running);


    double err0 = 0.0;
    for(index=0; index<n*n; index++) {
      if(fabs(c0[index]-c0[index])>err0) {
        printf("c0[%d]:%f, c0[%d]:%f\n",index, c0[index], index,c1[index]);
        err0 = fabs(c0[index]-c1[index]);
      }
    }

    printf("err:%f\n", err0);


    begin = clock();

    //dgemm2: more aggresive register reuse
    for(i = 0; i < n; i += 2) {
            for(j = 0; j < n; j += 2) {
                    for(k = 0; k < n; k += 2) {
                            c2[i*n + j] = a[i*n + k] * b[k*n + j] + a[i*n + k+1] * b[(k+1)*n + j] + c2[i*n + j];
                            c2[(i+1)*n + j] = a[(i+1)*n + k] * b[k*n + j] + a[(i+1)*n + k+1] * b[(k+1)*n + j] + c2[(i+1)*n + j];
                            c2[i*n + (j+1)] = a[i*n + k] * b[k*n + (j+1)] + a[i*n + k+1] * b[(k+1)*n + (j+1)] + c2[i*n + (j+1)];
                            c2[(i+1)*n + (j+1)] = a[(i+1)*n + k] * b[k*n + (j+1)] + a[(i+1)*n + k+1] * b[(k+1)*n + (j+1)] + c2[(i+1)*n + (j+1)];
                    }
            }

    }

    end = clock();
    running = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("The running time of dgemm2 with n=%d is %fs\n", n, running);


    begin = clock();

    //dgemm3
    for(i = 0; i < n; i += 2) {
            for(j = 0; j < n; j += 2)  {
                    register int t = i*n + j;
                    register int tt = t + n;
                    register double c00 = c3[t];
                    register double c01 = c3[t+1];
                    register double c10 = c3[tt];
                    register double c11 = c3[tt+1];

                    for(k = 0; k < n; k += 2) {
                            /* 2 by 2 mini matrix multiplication using registers*/
                            register int ta = i*n + k;
                            register int tta = ta + n;
                            register int tb = k*n + j;
                            register int ttb = tb + n;
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

                    c3[t] = c00;
                    c3[t+1] = c01;
                    c3[tt] = c10;
                    c3[tt+1] = c11;
            }

    }

    end = clock();
    running = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("The running time of dgemm3 with n=%d is %fs\n", n, running);

    double err1 = 0.0;
    for(index=0; index<n*n; index++) {
      if(fabs(c2[index]-c3[index])>err1) {
        printf("c2[%d]:%f, c3[%d]:%f\n",index, c2[index], index,c3[index]);
        err1 = fabs(c0[index]-c1[index]);
      }
    }
    printf("err:%f\n", err1);



  }else{
    printf("Please input n.\n");
  }


}
