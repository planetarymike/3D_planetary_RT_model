BINDIR=-I./bin/
BOOSTDIR=-I/home/mike/Documents/Utilities/boost_1_72_0/
EIGENDIR=-I/home/mike/Documents/Utilities/eigen-3.3.7/
IDIR=$(BINDIR) $(BOOSTDIR) $(EIGENDIR)

CC=g++
LIBS=-lm -fPIC
MPFLAGS=-fopenmp

#Eigen and the GNU scientific library need to be installed to compile and run this code. 

generate_source_function:
	$(CC) generate_source_function.cpp $(IDIR) $(LIBS) $(MPFLAGS) -O2 -ftree-vectorize -march=native -o generate_source_function.x

generate_source_function_profile:
	$(CC) generate_source_function.cpp $(IDIR) $(LIBS) $(MPFLAGS) -Og -g -p -o generate_source_function.x

generate_source_function_debug_warn:
	$(CC) generate_source_function.cpp $(IDIR) $(LIBS) -O0 -g -Wall -o generate_source_function.x
