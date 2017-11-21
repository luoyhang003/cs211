#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MyMPI.h"

#define MIN(a,b) ((a<b) ? (a) : (b))

int main(int argc, char* argv[]) {
        int id, p;
        unsigned long long int i, n, k;
        unsigned long long int size, low_value, high_value;
        unsigned long long int proc_0, prime, first;
        unsigned long long int global_count;
        double elapse_time;
        char* marked;
        unsigned long long int index;
	int count;

        unsigned long long int local_size;
        char* local_marked;

        int success = MPI_Init(&argc, &argv);
        if(success != MPI_SUCCESS) {
                printf("ERROR!\n");
        }

        MPI_Barrier(MPI_COMM_WORLD);
        elapse_time = -MPI_Wtime();
        MPI_Comm_rank(MPI_COMM_WORLD, &id);
        MPI_Comm_size(MPI_COMM_WORLD, &p);

	global_count = 0;

        if(argc != 3) {
                if(!id) {
                        printf("Please input correct command line parameters\n");
                }
                MPI_Finalize();
                exit(1);
        }

        n = atoll(argv[1]);
        low_value = 3 + BLOCK_LOW(id, p, n-2) + BLOCK_LOW(id, p, n-2)%2;
        high_value = 3 + BLOCK_HIGH(id, p, n-2) - BLOCK_HIGH(id, p, n-2)%2;
        size = (high_value-low_value)/2 + 1;
        proc_0 = (n-2)/(p*2);

        if((3+proc_0)<(int)sqrt((double)n)) {
                if(!id) {
                        printf("Too many processes\n");
                }
                MPI_Finalize();
                exit(1);
        }

        local_size = sqrt(n)+1;

        local_marked = (char*)malloc(local_size);

        for(i=0; i<local_size; i++) {
                local_marked[i] = 0;
        }

        k=2;

        do {
                for(i=k*k; i<local_size; i++) {
                        if(i % k == 0) {
                                local_marked[i] = 1;
                        }
                }
                while(local_marked[++k]) ;
        } while(k*k <= local_size);

       // unsigned long long int count_temp = 0;
       // for(i = 2; i<local_size; i++) {
       //         if(!local_marked[i]) {
       //                 count_temp++;
       //         }
       //}
       // printf("Number of primes from: 2 to %llu is %llu.\n", local_size, count_temp);



        marked = (char*)malloc(size);


        if(marked == NULL) {
                printf("Cannot Allocate Enough Memory\n");
                MPI_Finalize();
                exit(1);
        }

        for(i=0; i<size; i++) {
                marked[i] = 0;
        }


        index = 2;


        prime = 3;
	
	if(prime*prime > low_value) {
        	first = (prime* prime - low_value)/2;
	}else {
        	if(!(low_value%prime)) {
                	first = 0;
        	} else {
                	if((low_value%prime)%2 == 0) {
                        	first = prime - (low_value%prime)/2;
                	}
                	else {
                        	first = (prime - (low_value%prime))/2;
                	}
        	}
	}

	for(i = first; i < size; i+=prime) {
        	do {
                	marked[i] = 1;

                	if(prime*prime > low_value) {
                        	first = (prime* prime - low_value)/2;
                	}else {
                        	if(!(low_value%prime)) {
                                	first = 0;
                        	} else {
                                	if((low_value%prime)%2 == 0) {
                                        	first = prime - (low_value%prime)/2;
                                	}
                                	else {
                                        	first = (prime - (low_value%prime))/2;
                                	}
                        	}
                	}

                	while(local_marked[++index]) ;
                	prime = index;

        	} while(prime*prime<=n);
	}

        count = 0;


        for(i=0; i<size; i++) {
                if(!marked[i]) {
                        count++;
                }
        }

	//printf("My count is: %llu\n", count);
        if(p>1) {
                MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        elapse_time += MPI_Wtime();

        if(!id) {
                global_count++;
                printf("The total number of prime:%llu, total time:%10.6f, total node%s\n", global_count, elapse_time, argv[2]);
        }
        MPI_Finalize();
        return 0;
}
