BINDIR=-I./bin/
BOOSTDIR=-I/home/mike/Documents/Utilities/boost_1_72_0/
EIGENDIR=-I/home/mike/Documents/Utilities/eigen-3.3.7/
IDIR=$(BINDIR) $(BOOSTDIR) $(EIGENDIR)

# GNU Compiler
CC=g++
LIBS=-lm 
MPFLAGS=-fopenmp
OFLAGS=-Ofast -march=native -DNDEBUG

# Nvidia CUDA Compiler
NCC=nvcc --disable-warnings
NLIBS=-lm
NOFLAGS= -O3 -DNDEBUG

# # intel compiler
# # you may need to run this
# # source /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux
# CC=icpc
# COMPILER_OPT=-std=c++17
# LIBS=-fPIC -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
# MPFLAGS=-qopenmp
# OFLAGS=-O3 -march=native -DNDEBUG -fp-model fast=2

generate_source_function:
	$(CC) generate_source_function.cpp $(COMPILER_OPT) $(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) -o generate_source_function.x

generate_source_function_gpu:
	$(NCC) generate_source_function_gpu.cu $(COMPILER_OPT) $(IDIR) $(NLIBS) $(NOFLAGS) -o generate_source_function_gpu.x

generate_source_function_profile:
	$(CC) generate_source_function.cpp $(COMPILER_OPT) $(IDIR) $(LIBS) $(OFLAGS) -g -o generate_source_function.x

# generate_source_function_gprof:
# 	$(CC) generate_source_function.cpp $(COMPILER_OPT) $(IDIR) $(LIBS) $(OFLAGS) -pg -o generate_source_function.x

generate_source_function_debug_warn:
	$(CC) generate_source_function.cpp $(COMPILER_OPT) $(IDIR) $(LIBS) -O0 -g -Wall -o generate_source_function.x
