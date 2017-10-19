/*
 * Course: CS 211 - High Performance Computing
 * Project 1: Performance Optimization via Register and Cache Reuse
 *
 * Author: Yuanhang Luo
 * Create Date: Oct 14th, 2017
 * Revise Date: Oct 14th, 2017
 *
 */

 #include "stdio.h"
 #include "string.h"
 #include "time.h"
 #include "stdlib.h"
 #include "float.h"
 #include "math.h"
 #include "time.h"


int main(int argc, char* argv[]) {

        int n = 2048;
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

        int index;

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
        for (i=0; i<n; i++) {
                for (j=0; j<n; j++) {
                        register double r=c_ijk[i*n+j];
                        for (k=0; k<n; k++) {
                                r += a[i*n+k] * b[k*n+j];
                        }
                        c_ijk[i*n+j]=r;
                }
        }
        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of ijk is %f.\n", running);

        //ikj
        begin = clock();
        for (i=0; i<n; i++) {
                for (k=0; k<n; k++) {
                        register double r = a[i*n+k];
                        for (j=0; j<n; j++)
                                c_ikj[i*n+j] += r * b[k*n+j];
                }
        }

        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of ikj is %f.\n", running);

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
        for (j=0; j<n; j++) {
                for (i=0; i<n; i++) {
                        register double sum = 0.0;
                        for (k=0; k<n; k++)
                                sum += a[i*n+k] * b[k*n+j];
                        c_jik[i*n+j] = sum;
                }
        }
        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of jik is %f.\n", running);

        for(index=0,err=0.0; index<n*n; index++) {
                if(fabs(c_ijk[index]-c_jik[index])>err) {
                        printf("c_ijk[%d]:%f, c_ikj[%d]:%f\n",index, c_ijk[index], index,c_jik[index]);
                        err = fabs(c_ijk[index]-c_jik[index]);
                }
        }

        printf("err:%f\n", err);

        //jki
        begin = clock();
        for (j=0; j<n; j++) {
                for (k=0; k<n; k++) {
                        register double r = b[k*n + j];
                        for (i=0; i<n; i++)
                                c_jki[i*n+j] += a[i*n+k] * r;
                }
        }

        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of jki is %f.\n", running);

        for(index=0,err=0.0; index<n*n; index++) {
                if(fabs(c_ijk[index]-c_jki[index])>err) {
                        printf("c_ijk[%d]:%f, c_jki[%d]:%f\n",index, c_ijk[index], index,c_jki[index]);
                        err = fabs(c_ijk[index]-c_jki[index]);
                }
        }

        printf("err:%f\n", err);


        //kij
        begin = clock();
        for (j=0; j<n; j++) {
                for (i=0; i<n; i++) {
                        register double sum = 0.0;
                        for (k=0; k<n; k++)
                                sum += a[i*n+k] * b[k*n+j];
                        c_kij[i*n+j] = sum;
                }
        }

        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of kij is %f.\n", running);

        for(index=0,err=0.0; index<n*n; index++) {
                if(fabs(c_ijk[index]-c_kij[index])>err) {
                        printf("c_ijk[%d]:%f, c_kij[%d]:%f\n",index, c_ijk[index], index,c_kij[index]);
                        err = fabs(c_ijk[index]-c_kij[index]);
                }
        }

        printf("err:%f\n", err);


        //kji

        begin = clock();
        for (k=0; k<n; k++) {
                for (j=0; j<n; j++) {
                        register double r = b[k*n+j];
                        for (i=0; i<n; i++)
                                c_kji[i*n+j] += a[i*n+k] * r;
                }
        }

        end = clock();
        running = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The running time of kji is %f.\n", running);

        for(index=0,err=0.0; index<n*n; index++) {
                if(fabs(c_ijk[index]-c_kji[index])>err) {
                        printf("c_ijk[%d]:%f, c_kji[%d]:%f\n",index, c_ijk[index], index,c_kji[index]);
                        err = fabs(c_ijk[index]-c_kji[index]);
                }
        }


}
