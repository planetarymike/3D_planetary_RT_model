#files that need compilin'
OBJDIR = ./bin
SRCDIR = src
SRCFILES = $(wildcard $(SRCDIR)/atm/*.cpp) $(wildcard $(SRCDIR)/*.cpp) 

NSRCFILES = $(SRCFILES) generate_source_function.cpp
NOBJFILES    := $(filter %.o, $(NSRCFILES:%.cpp=$(OBJDIR)/%.cuda.o)       $(NSRCFILES:%.cu=$(OBJDIR)/%.cuda.o      ))
NOBJFILESDBG := $(filter %.o, $(NSRCFILES:%.cpp=$(OBJDIR)/%.cuda.debug.o) $(NSRCFILES:%.cu=$(OBJDIR)/%.cuda.debug.o))

#include directories
BOOSTDIR=-I/home/mike/Documents/Utilities/boost_1_73_0/
EIGENDIR=-I/home/mike/Documents/Utilities/eigen-3.3.7/
IDIR=-I$(SRCDIR) $(BOOSTDIR) $(EIGENDIR)

# GNU Compiler
CC=g++
LIBS=-lm 
MPFLAGS=-fopenmp
#OFLAGS=-Ofast -march=native -DNDEBUG
OFLAGS=-O3 -march=native -DNDEBUG

# Nvidia CUDA Compiler
NCC=nvcc --disable-warnings
NIDIR=$(IDIR) -L/usr/local/cuda-10.2/lib64/ -I/home/mike/Documents/Utilities/cuda-samples/Common/
NLIBS=-lm -lcudart 
NOFLAGS= -O3 -DNDEBUG #-lineinfo

# # intel compiler
# you may need to run this
# source /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux
ICC=icpc
ICOMPILER_OPT=-std=c++17
ILIBS=-fPIC -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
IMPFLAGS=-qopenmp
IOFLAGS=-march=native -DNDEBUG -fast

generate_source_function:
	$(CC) generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) -o generate_source_function.x

generate_source_function_intel:
	$(ICC) generate_source_function.cpp $(SRCFILES) $(ICOMPILER_OPT) $(IDIR) $(ILIBS) $(IMPFLAGS) $(IOFLAGS) -o generate_source_function.x

generate_source_function_profile:
	$(CC) generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) $(OFLAGS) -g -o generate_source_function.x

generate_source_function_debug_warn:
	$(CC) generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) -O0 -g -Wall -o generate_source_function.x




generate_source_function_gpu: $(NOBJFILES) 
	@echo $(NOBJFILES)
	@echo "linking..."
	@$(NCC) $(NOBJFILES) $(NIDIR) $(NLIBS) $(NOFLAGS) -o generate_source_function_gpu.x

$(OBJDIR)/%.cuda.o: %.cpp
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) -x cu -D RT_FLOAT $(NIDIR) $(NLIBS) $(NOFLAGS) -dc $< -o $@

$(OBJDIR)/%.cuda.o: %.cu
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) -x cu -D RT_FLOAT $(NIDIR) $(NLIBS) $(NOFLAGS) -dc $< -o $@




generate_source_function_gpu_debug: $(NOBJFILESDBG) 
	@echo "linking ..."
	@$(NCC) $(NOBJFILESDBG) $(NIDIR) $(NLIBS) $(NOFLAGS) -o generate_source_function_gpu.x

$(OBJDIR)/%.cuda.debug.o: %.cpp
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) -x cu -D RT_FLOAT $(NIDIR) $(NLIBS) -O0 -g -G -dc $< -o $@

$(OBJDIR)/%.cuda.debug.o: %.cu
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) -x cu -D RT_FLOAT $(NIDIR) $(NLIBS) -O0 -g -G -dc $< -o $@


