#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char *argv[]) {
    int n, k, i;
    char* marked;

    n = atoi(argv[1]);

    marked = (char *) malloc(n);

    k=2;

    for(i=0; i<n; i++) {
        marked[i] = 0;
    }


    do{
        for(i=k*k; i<n; i++) {
            if(i % k == 0) {
                marked[i] = 1;
            }
        }
        while(marked[++k]);

    }while (k*k <= n);

    int count = 0;

    for(i=2; i<n; i++) {
        if(!marked[i]) {
            count++;
        }
    }

    printf("count: %d\n", count);


}
