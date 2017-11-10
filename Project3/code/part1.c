#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MyMPI.h"

#define MIN(a,b) ((a<b)?(a):(b))

int main(int argc, char* argv[]) {
	int i, n, low_value, high_value,proc_0, id, p, size, count, global_count;
	double elapse_time;
	char* marked;
	int index, prime, first;	

	int success = MPI_Init(&argc, &argv);
	if(success != MPI_SUCCESS) {
		printf("ERROR!\n");
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	elapse_time = -MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	if(argc != 2) {
		if(!id) {
			printf("Command Line: %s\n", argv[0]);
		}
		MPI_Finalize();
		exit(1);
	}

	n = atoi(argv[1]);
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

	marked = (char*)malloc(size);
	
	if(marked == NULL){
		printf("Cannot Allocate Enough Memory\n");
		MPI_Finalize();
		exit(1);
	}
	
	for(i=0; i<size; i++) {
		marked[i] = 0;
	}

	if(!id) {
		index = 0;
	}

	prime = 3;

	do {
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
		for(i=first; i<size; i+=prime) {
			marked[i] = 1;
		}
		if(!id) {
			while(marked[++index]);
			prime = index*2+3;
		}
		if(p>1) {
			MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}
	} while(prime*prime<=n);

	count = 0;
	
	for(i=0; i<size; i++){
		if(!marked[i]) {
			count++;
		}
	}
	if(p>1) {
		MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	elapse_time += MPI_Wtime();

	if(!id) {
		global_count++;
		printf("%d primes are less than or equal to %d\n", global_count, n);
		printf("Total elapsed time: %10.6f", elapse_time);
	}
	MPI_Finalize();	
	return 0;
}
