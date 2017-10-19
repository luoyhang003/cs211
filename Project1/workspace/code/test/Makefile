all: register.o cache_reuse_simple.o cache_reuse_block.o
	gcc register.o -o register
	gcc cache_reuse_simple.o -o cache_simple
	gcc cache_reuse_block.o -o cache_block
register : register.c
	gcc -c register.c
cache : cache_reuse_simple.c cache_reuse_block.c
	gcc -c cache_reuse_simple.c
	gcc -c cache_reuse_block.c
clean :
	rm -f *.o register cache_simple cache_block
