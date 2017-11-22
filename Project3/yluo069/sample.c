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
	low_value = 2 + BLOCK_LOW(id, p, n-1);
	high_value = 2 + BLOCK_HIGH(id, p, n-1);
	size = BLOCK_SIZE(id, p, n-1);
	proc_0 = (n-1)/p;
	
	if((2+proc_0)<(int)sqrt((double)n)) {
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

	prime = 2;

	do {
		if(prime*prime > low_value) {
			first = prime* prime - low_value;
		}else {
			if(!(low_value%prime)) {
				first = 0;
			} else {
				first = prime-(low_value%prime);
			}
		}
		for(i=first; i<size; i+=prime) {
			marked[i] = 1;
		}
		if(!id) {
			while(marked[++index]);
			prime = index+2;
		}
		MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
	} while(prime*prime<=n);

	count = 0;
	
	for(i=0; i<size; i++){
		if(!marked[i]) {
			count++;
		}
	}
	
	MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	elapse_time += MPI_Wtime();

	if(!id) {
		printf("%d primes are less than or equal to %d\n", global_count, n);
		printf("Total elapsed time: %10.6f", elapse_time);
	MPI_Finalize();
	return 0;
	}
	
	return 0;
}
