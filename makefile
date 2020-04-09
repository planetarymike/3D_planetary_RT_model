BINDIR=-I./bin/
BOOSTDIR=-I/home/mike/Documents/Utilities/boost_1_72_0/
EIGENDIR=-I/home/mike/Documents/Utilities/eigen-3.3.7/
IDIR=$(BINDIR) $(BOOSTDIR) $(EIGENDIR)

# GNU Compiler
CC=g++
LIBS=-lm -fPIC
MPFLAGS=-fopenmp -fopenacc
OFLAGS=-O3 -ftree-vectorize -march=native

# # intel compiler
# # you may need to run this
# # source /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux
# CC=icpc
# COMPILER_OPT=-std=c++17 -mkl
# LIBS=-lm -fPIC
# MPFLAGS=-qopenmp
# OFLAGS=-O3 -march=native 

generate_source_function:
	$(CC) generate_source_function.cpp $(COMPILER_OPT) $(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) -o generate_source_function.x

generate_source_function_noassert:
	$(CC) generate_source_function.cpp $(COMPILER_OPT) $(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) -DNDEBUG -o generate_source_function.x

generate_source_function_profile:
	$(CC) generate_source_function.cpp $(COMPILER_OPT) $(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) -g -o generate_source_function.x

generate_source_function_gprof:
	$(CC) generate_source_function.cpp $(COMPILER_OPT) $(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) -pg -o generate_source_function.x

generate_source_function_debug_warn:
	$(CC) generate_source_function.cpp $(COMPILER_OPT) $(IDIR) $(LIBS) $(MPFLAGS) -O0 -g -Wall -o generate_source_function.x
