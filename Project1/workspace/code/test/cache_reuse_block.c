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
        int n = 2000;
        printf("******  n=%d    *******\n", n);

        double *a = NULL;
        double *b = NULL;


        double *c_ijk = NULL;
        double *c_ikj = NULL;
        double *c_jik = NULL;
        double *c_jki = NULL;
        double *c_kij = NULL;
        double *c_kji = NULL;

        a = calloc(sizeof(double), n*n);
        b = calloc(sizeof(double), n*n);


        c_ijk = calloc(sizeof(double), n*n);
        c_ikj = calloc(sizeof(double), n*n);
        c_jik = calloc(sizeof(double), n*n);
        c_jki = calloc(sizeof(double), n*n);
        c_kij = calloc(sizeof(double), n*n);
        c_kji = calloc(sizeof(double), n*n);

        int i, j, k;
        int i1, j1, k1;

        int index;

        int B = 10;
        printf("******  B=%d    *******\n", B);

        double max, min;
        max = n;
        min = 0.0;
        double range = 0.0+(max - min);
        double div = RAND_MAX / range;
        srand(time(NULL));

        for(index = 0; index < n*n; index++) {
                a[index] = min + (rand()/div);
                b[index] = min + (rand()/div);
        }

        clock_t begin;
        clock_t end;
        double running;

        //ijk
        begin = clock();

        for (i = 0; i < n; i+=B) {
                for (j = 0; j < n; j+=B) {
                        for (k = 0; k < n; k+=B) {
                                /* B x B mini matrix multiplications */
                                for (i1 = i; i1 < i+B; i1++) {
                                        for (j1 = j; j1 < j+B; j1++) {
                                                register double r=c_ijk[i1*n+j1];

                                                for (k1 = k; k1 < k+B; k1++) {
                                                        c_ijk[i1*n+j1] += a[i1*n+k1]*b[k1*n+j1];
                                                }
                                        }
                                }
                        }
                }
        }
        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of ijk with B = 10 is %f.\n", running);

        //ikj

        begin = clock();

        for (i = 0; i < n; i+=B) {
                for (k = 0; k < n; k+=B) {
                        for (j = 0; j < n; j+=B) {
                                /* B x B mini matrix multiplications */
                                for (i1 = i; i1 < i+B; i1++) {
                                        for (k1 = k; k1 < k+B; k1++) {
                                                register double r = a[i1*n+k1];
                                                for (j1 = j; j1 < j+B; j1++) {
                                                        c_ikj[i1*n+j1] += r * b[k1*n+j1];
                                                }
                                        }
                                }
                        }
                }
        }
        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of ikj with B = 10 is %f.\n", running);

        double err;
        for(index=0,err=0.0; index<n*n; index++) {
                if(fabs(c_ijk[index]-c_ikj[index])>err) {
                        printf("c_ijk[%d]:%f, c_ikj[%d]:%f\n",index, c_ijk[index], index,c_ikj[index]);
                        err = fabs(c_ijk[index]-c_ikj[index]);
                }
        }

        printf("err:%f\n", err);


        //jik
        begin = clock();

        for (j = 0; j < n; j+=B) {
                for (i = 0; i < n; i+=B) {
                        for (k = 0; k < n; k+=B) {
                                /* B x B mini matrix multiplications */
                                for (j1 = j; j1 < j+B; j1++) {
                                        for (i1 = i; i1 < i+B; i1++) {
                                                register double sum = 0.0;

                                                for (k1 = k; k1 < k+B; k1++) {
                                                        sum += a[i1*n+k1] * b[k1*n+j1];
                                                }
                                                c_jik[i1*n+j1] = sum;

                                        }
                                }
                        }
                }
        }
        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of jik with B = 10 is %f.\n", running);

        // for(index=0,err=0.0; index<n*n; index++) {
        //         if(fabs(c_ijk[index]-c_jik[index])>err) {
        //                 printf("c_ijk[%d]:%f, c_jik[%d]:%f\n",index, c_ijk[index], index,c_jik[index]);
        //                 err = fabs(c_ijk[index]-c_jik[index]);
        //         }
        // }
        //
        // printf("err:%f\n", err);


        //jki

        begin = clock();

        for (j = 0; j < n; j+=B) {
                for (k = 0; k < n; k+=B) {
                        for (i = 0; i < n; i+=B) {
                                /* B x B mini matrix multiplications */
                                for (j1 = j; j1 < j+B; j1++) {
                                        for (k1 = k; k1 < k+B; k1++) {
                                                register double r = b[k1*n + j1];

                                                for (i1 = i; i1 < i+B; i1++) {
                                                        c_jki[i1*n+j1] += a[i1*n+k1] * r;
                                                }
                                        }
                                }
                        }
                }
        }
        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of jki with B = 10 is %f.\n", running);

        for(index=0,err=0.0; index<n*n; index++) {
                if(fabs(c_ijk[index]-c_jki[index])>err) {
                        printf("c_ijk[%d]:%f, c_jki[%d]:%f\n",index, c_ijk[index], index,c_jki[index]);
                        err = fabs(c_ijk[index]-c_jki[index]);
                }
        }

        printf("err:%f\n", err);

        //kij

        begin = clock();

        for (k = 0; k < n; k+=B) {
                for (i = 0; i < n; i+=B) {
                        for (j = 0; j < n; j+=B) {
                                /* B x B mini matrix multiplications */
                                for (k1 = k; k1 < k+B; k1++) {
                                        for (i1 = i; i1 < i+B; i1++) {
                                                register double r = a[i1*n+k1];

                                                for (j1 = j; j1 < j+B; j1++) {
                                                        c_kij[i1*n+j1] += r * b[k1*n+j1];
                                                }
                                        }
                                }
                        }
                }
        }
        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of kij with B = 10 is %f.\n", running);

        for(index=0,err=0.0; index<n*n; index++) {
                if(fabs(c_ijk[index]-c_kij[index])>err) {
                        printf("c_ijk[%d]:%f, c_jki[%d]:%f\n",index, c_ijk[index], index,c_kij[index]);
                        err = fabs(c_ijk[index]-c_kij[index]);
                }
        }

        printf("err:%f\n", err);


        //kji
        begin = clock();

        for (k = 0; k < n; k+=B) {
                for (j = 0; j < n; j+=B) {
                        for (i = 0; i < n; i+=B) {
                                /* B x B mini matrix multiplications */
                                for (k1 = k; k1 < k+B; k1++) {
                                        for (j1 = j; j1 < j+B; j1++) {
                                                register double r = b[k1*n+j1];

                                                for (i1 = i; i1 < i+B; i1++) {
                                                    c_kji[i1*n+j1] += a[i1*n+k1] * r;

                                                }
                                        }
                                }
                        }
                }
        }
        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of kij with B = 10 is %f.\n", running);

        for(index=0,err=0.0; index<n*n; index++) {
                if(fabs(c_ijk[index]-c_kji[index])>err) {
                        printf("c_ijk[%d]:%f, c_jki[%d]:%f\n",index, c_ijk[index], index,c_kji[index]);
                        err = fabs(c_ijk[index]-c_kji[index]);
                }
        }

        printf("err:%f\n", err);



}
