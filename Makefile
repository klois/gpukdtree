CXX = g++
CC = gcc -fms-extensions 
CFLAGS = -Wall -g -DSOFTENING_PLUMMER -O3 -fopenmp
INCLUDE = -Iinclude
LIBS =  -lm -lrt
OPENCL = -I$(OPENCL_ROOT)/include -L$(OPENCL_ROOT)/lib/x86_64 -lOpenCL -D_POSIX_C_SOURCE=199309
CUDA = -I$(CUDA_PATH)/include -L$(CUDA_PATH)/lib64 -lcuda -lcudart -D_POSIX_C_SOURCE=199309 -DCUDA

all: gpukdtree 

cuda: gpukdtreecuda

lib_icl: src/lib_icl_ext.c src/lib_icl.c
	$(CC) $(CFLAGS) src/lib_icl.c  $(INCLUDE) $(OPENCL) -std=c99 -c -o bin/lib_icl.o 
	$(CC) $(CFLAGS) src/lib_icl_ext.c $(INCLUDE) $(OPENCL) -std=c99 -c -o bin/lib_icl_ext.o 

lib_icl_cuda: src/lib_icl_ext.c src/lib_icl_cuda.c
	$(CC) $(CFLAGS) src/lib_icl_cuda.c  $(INCLUDE) $(CUDA) -std=c99 -c -o bin/lib_icl_cuda.o 
	$(CC) $(CFLAGS) src/lib_icl_ext.c $(INCLUDE) $(CUDA) -std=c99 -c -o bin/lib_icl_ext.o 

particle: src/snapshot_io.c src/data_ops.c include/snapshot_io.h
	$(CC) $(CFLAGS) src/snapshot_io.c  $(INCLUDE) $(LIBS) -c -o bin/snapshot_io.o -DDEBUG
	$(CC) $(CFLAGS) src/data_ops.c  $(INCLUDE) $(LIBS) -c -o bin/data_ops.o
	
gpukdtree: src/gpukdtree.c src/gpuKdtreeWalk.c src/timer.c particle lib_icl
	$(CC) $(CFLAGS) src/gpukdtree.c src/gpuKdtreeWalk.c src/timer.c src/forcetest.c src/boltScan.c bin/snapshot_io.o bin/data_ops.o bin/lib_icl.o bin/lib_icl_ext.o -std=c99 $(INCLUDE) $(LIBS) $(OPENCL) -o bin/gpukdtree

gpukdtreecuda: src/gpukdtree.c src/gpuKdtreeWalk.c src/timer.c particle lib_icl_cuda
	$(CC) $(CFLAGS) src/gpukdtree.c src/gpuKdtreeWalk.c src/timer.c src/forcetest.c src/boltScan.c bin/snapshot_io.o bin/data_ops.o bin/lib_icl_cuda.o bin/lib_icl_ext.o -std=c99 $(INCLUDE) $(LIBS) $(CUDA) -o bin/gpukdtree

clean:
	rm bin/gpukdtree bin/*.o
