MAKEFLAGS += -j20 #parallel compilation

# you may need these in your .bashrc
# export CUDA_HOME=/usr/local/cuda-11.2
# export PATH=$CUDA_HOME/bin${PATH:+:${PATH}}$ 
# export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

#files that need compilin'
OBJDIR = ./bin
SRCDIR = ./src
SRCDIRS = ./src ./src/atm ./src/emission ./src/grid
PSRCFILES = $(foreach dir,$(SRCDIRS),$(wildcard $(dir)/*.cpp))
SRCFILES = $(filter-out ./src/observation_fit.cpp, $(PSRCFILES))

NSRCFILES = $(SRCFILES) generate_source_function.cpp
NOBJFILES    := $(filter %.o, $(NSRCFILES:%.cpp=$(OBJDIR)/%.cuda.o))
NOBJFILESDBG := $(filter %.o, $(NSRCFILES:%.cpp=$(OBJDIR)/%.cuda.debug.o))


PYSRCFILES = $(SRCFILES) ./src/observation_fit.cpp $(wildcard $(SRCDIR)/quemerais_IPH_model/*.cpp)
PYOBJFILES   := $(filter %.o, $(PYSRCFILES:%.cpp=$(OBJDIR)/%.o))
PYNOBJFILES  := $(filter %.o, $(PYSRCFILES:%.cpp=$(OBJDIR)/%.cuda.o))


#include directories
BOOSTDIR=-I/home/mike/Documents/Utilities/boost_1_73_0/
EIGENDIR=-I/home/mike/Documents/Utilities/eigen_git/
IDIR= $(foreach dir,$(SRCDIRS),-I$(abspath $(dir))) $(BOOSTDIR) $(EIGENDIR)

# GNU Compiler
CC=g++ -std=c++17 -fPIC #-D RT_FLOAT -Wfloat-conversion #these commands can be used to check for double literals
LIBS=-lm -lgomp
MPFLAGS=-fopenmp
OFLAGS=-O3 -march=native -DNDEBUG 

# Nvidia CUDA Compiler
#device spec
CUDA_DEVICE_CODE=$(shell $$CUDA_HOME/extras/demo_suite/deviceQuery | grep -o 'CUDA Capability Major/Minor version number:.*' | cut -f2 -d ':' | sed -r 's/\s+//g' | sed 's/\.//')

NCC=nvcc -std=c++17 -Xcompiler -fPIC -Xcudafe --display_error_number #--disable-warnings
NFLAGS=-x cu -D RT_FLOAT              -D EIGEN_NO_CUDA                -D BOOST_NO_CUDA
#            ^^^^32-bit calculation   ^^^^^ disable Eigen on device   ^^^^^ disable Boost on device
#                                           (needs Eigen git repo     (added this flag by hand as a wrapper
#                                            master > 21 Oct 2020      around BOOST_GPU_ENABLED in
#                                                                      boost/config/compiler/nvcc.hpp)
NIDIR=$(IDIR) \
      -L$(CUDA_HOME)/lib64/ \
      -I$(CUDA_HOME)/samples/common/inc/
NLIBS=-lm -lcudart -lcusolver -lcublas
NOBASEFLAGS= -O3 -DNDEBUG -lineinfo --use_fast_math #--maxrregcount 43
# if we are CUDA 11, link time optimization is possible
ifeq ($(shell nvcc --version | grep -o 'release.*' | cut -f2 -d' ' | cut -f1 -d.),11)
CUDA_DLTO=true
EXTRA_NOFLAGS = -dlto --gpu-architecture=sm_$(CUDA_DEVICE_CODE)
else
EXTRA_NOFLAGS = --gpu-architecture=sm_$(CUDA_DEVICE_CODE)
endif
NOFLAGS=$(NOBASEFLAGS) $(EXTRA_NOFLAGS)
NDBGFLAGS=-O0 -g -G -arch sm_$(CUDA_DEVICE_CODE) # -lineinfo
#                ^^^ this -G sometimes changes the behavior of the code??

generate_source_function:
	$(CC) generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) -o generate_source_function.x

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


generate_source_function_gpu_debug: $(NOBJFILESDBG) 
	@echo "linking ..."
	@$(NCC) $(NOBJFILESDBG) $(NIDIR) $(NLIBS) $(NDBGFLAGS) -o generate_source_function_gpu.x

$(OBJDIR)/%.cuda.debug.o: %.cpp
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NDBGFLAGS) -dc $< -o $@



py_corona_sim_cpu:
	@mkdir -p bin
	gfortran -fPIC -Ofast -c -std=legacy\
	  $(SRCDIR)/quemerais_IPH_model/ipbackgroundCFR_fun.f \
	  -o $(OBJDIR)/ipbackgroundCFR_fun.o

	@cd python; \
	export IDIR='$(IDIR)'; \
	export COMPILE_FLAGS='$(MPFLAGS) $(OFLAGS)' ; \
	export COMPILED_OBJECTS='$(PYOBJFILES)'; \
	export LIBRARIES='$(LIBS)'; \
	export SOURCE_FILES='$(PYSRCFILES)'; \
	export CC='g++'; \
	export CXX='g++'; \
	python setup_corona_sim.py build_ext --inplace

$(OBJDIR)/%.o: %.cpp
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(CC) -c $< $(IDIR) $(LIBS) $(OFLAGS) -o $@

py_corona_sim_gpu: 
	@mkdir -p bin
	gfortran -fPIC -Ofast -c -std=legacy\
	  $(SRCDIR)/quemerais_IPH_model/ipbackgroundCFR_fun.f \
	  -o $(OBJDIR)/ipbackgroundCFR_fun.o

	@cd python; \
	export IDIR='$(NIDIR)'; \
	export COMPILE_FLAGS='$(NFLAGS) $(NOFLAGS)' ; \
	export CUDA_DLTO='$(NOFLAGS)' ; \
	export COMPILED_OBJECTS='$(PYNOBJFILES)'; \
	export LIBRARIES='$(NLIBS)'; \
	export SOURCE_FILES='$(PYSRCFILES)'; \
	export CC='nvcc'; \
	export CXX='nvcc'; \
	python setup_corona_sim.py build_ext --inplace -RT_FLOAT



clean_gpu:
	rm -f generate_source_function_gpu.x
	test -d ./bin && find ./bin/ -type f -name '*cuda*' -delete 2>/dev/null

clean_all:
	rm -f generate_source_function.x
	rm -f generate_source_function_gpu.x
	rm -rf bin
	rm -rf python/build* python/*.cpp #python/*.so
