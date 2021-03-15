MAKEFLAGS += -j20 #parallel compilation

# you may need these in your .bashrc
# export CUDA_HOME=/usr/local/cuda-11.0
# export PATH=$CUDA_HOME/bin${PATH:+:${PATH}}$ 
# export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

#files that need compilin'
OBJDIR = ./bin
SRCDIRS = src src/atm src/emission src/grid
PSRCFILES = $(foreach dir,$(SRCDIRS),$(wildcard $(dir)/*.cpp))
SRCFILES = $(filter-out src/observation_fit.cpp, $(PSRCFILES))

NSRCFILES = $(SRCFILES) generate_source_function.cpp
NOBJFILES    := $(filter %.o, $(NSRCFILES:%.cpp=$(OBJDIR)/%.cuda.o)       $(NSRCFILES:%.cu=$(OBJDIR)/%.cuda.o      ))
NOBJFILESDBG := $(filter %.o, $(NSRCFILES:%.cpp=$(OBJDIR)/%.cuda.debug.o) $(NSRCFILES:%.cu=$(OBJDIR)/%.cuda.debug.o))

#include directories
BOOSTDIR=-I/home/mike/Documents/Utilities/boost_1_73_0/
EIGENDIR=-I/home/mike/Documents/Utilities/eigen_git/
IDIR= $(foreach dir,$(SRCDIRS),-I$(dir)) $(BOOSTDIR) $(EIGENDIR)

# GNU Compiler
CC=g++ -std=c++17 #-D RT_FLOAT -Wfloat-conversion #these commands can be used to check for double literals
LIBS=-lm 
MPFLAGS=-fopenmp
OFLAGS=-O3 -march=native -DNDEBUG 

# Nvidia CUDA Compiler
NCC=nvcc -Xcudafe --display_error_number #--disable-warnings
NFLAGS=-x cu -D RT_FLOAT              -D EIGEN_NO_CUDA                -D BOOST_NO_CUDA
#            ^^^^32-bit calculation   ^^^^^ disable Eigen on device   ^^^^^ disable Boost on device
#                                           (needs Eigen git repo     (added this flag by hand as a wrapper
#                                            master > 21 Oct 2020      around BOOST_GPU_ENABLED in
#                                                                      boost/config/compiler/nvcc.hpp)
NIDIR=$(IDIR) \
      -L$(CUDA_HOME)/lib64/ \
      -I$(CUDA_HOME)/samples/common/inc/
NLIBS=-lm -lcudart -lcusolver -lcublas
NOBASEFLAGS= -O3 -DNDEBUG -lineinfo -arch sm_61 --use_fast_math #--maxrregcount 43
# if we are CUDA 11, link time optimization is possible
ifeq ($(shell nvcc --version | grep -o 'release.*' | cut -f2 -d' ' | cut -f1 -d.),11)
CUDA_DLTO=true
NOFLAGS=$(NOBASEFLAGS) -dlto
else
NOFLAGS=$(NOBASEFLAGS)
endif
NDBGFLAGS=-O0 -g -G -arch sm_61 # -lineinfo
#                ^^^ this -G sometimes changes the behavior of the code??

# # intel compiler
# you may need to run this
# source /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux
ICC=icpc
ICOMPILER_OPT=-std=c++17
ILIBS=-fPIC -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
IMPFLAGS=-qopenmp
IOFLAGS=-DNDEBUG -fast

all: generate_source_function generate_source_function_gpu

generate_source_function:
	$(CC) generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) -o generate_source_function.x

generate_source_function_intel:
	$(ICC) generate_source_function.cpp $(SRCFILES) $(ICOMPILER_OPT) $(IDIR) $(ILIBS) $(IMPFLAGS) $(IOFLAGS) -o generate_source_function.x

generate_source_function_profile:
	$(CC) generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) $(OFLAGS) -g -o generate_source_function.x

generate_source_function_debug_warn:
	$(CC) generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) -O0 -g -Wall -Wextra -o generate_source_function.x

generate_source_function_float:
	$(CC) -D RT_FLOAT generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) -o generate_source_function.x



generate_source_function_gpu: $(NOBJFILES) 
	@echo "linking..."
ifeq ($(CUDA_DLTO),true)
	$(info Using CUDA 11 link time optimization)
endif
	@$(NCC) $(NOBJFILES) $(NIDIR) $(NLIBS) $(NOFLAGS) -o generate_source_function_gpu.x

$(OBJDIR)/%.cuda.o: %.cpp
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NOFLAGS) -dc $< -o $@

$(OBJDIR)/%.cuda.o: %.cu
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NOFLAGS) -dc $< -o $@




generate_source_function_gpu_debug: $(NOBJFILESDBG) 
	@echo "linking ..."
	@$(NCC) $(NOBJFILESDBG) $(NIDIR) $(NLIBS) $(NDBGFLAGS) -o generate_source_function_gpu.x

$(OBJDIR)/%.cuda.debug.o: %.cpp
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NDBGFLAGS) -dc $< -o $@

$(OBJDIR)/%.cuda.debug.o: %.cu
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NDBGFLAGS) -dc $< -o $@


clean_gpu:
	rm -f generate_source_function_gpu.x
	test -d ./bin && find ./bin/ -type f -name '*cuda*' -delete 2>/dev/null

clean_all:
	rm -f generate_source_function.x
	rm -f generate_source_function_gpu.x
	rm -rf bin
