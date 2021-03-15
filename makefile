MAKEFLAGS += -j20 #parallel compilation

# you may need these in your .bashrc
# export CUDA_HOME=/usr/local/cuda-11.0
# export PATH=$CUDA_HOME/bin${PATH:+:${PATH}}$ 
# export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

#files that need compilin'
OBJDIR = ./bin
SRCDIR = ./src
SRCDIRS = src src/atm src/emission src/grid
PSRCFILES = $(foreach dir,$(SRCDIRS),$(wildcard $(dir)/*.cpp))
SRCFILES = $(filter-out src/observation_fit.cpp, $(PSRCFILES))

NSRCFILES = $(SRCFILES) generate_source_function.cpp
NOBJFILES    := $(filter %.o, $(NSRCFILES:%.cpp=$(OBJDIR)/%.cuda.o)       $(NSRCFILES:%.cu=$(OBJDIR)/%.cuda.o      ))
NOBJFILESDBG := $(filter %.o, $(NSRCFILES:%.cpp=$(OBJDIR)/%.cuda.debug.o) $(NSRCFILES:%.cu=$(OBJDIR)/%.cuda.debug.o))


PYSRCFILES = $(SRCFILES) src/observation_fit.cpp $(wildcard $(SRCDIR)/quemerais_IPH_model/*.cpp)
PYOBJFILES   := $(filter %.o, $(PYSRCFILES:%.cpp=$(OBJDIR)/%.o))
PYNOBJFILES  := $(filter %.o, $(PYSRCFILES:%.cpp=$(OBJDIR)/%.cuda.o))


#include directories
BOOSTDIR=-I/home/mike/Documents/Utilities/boost_1_73_0/
EIGENDIR=-I/home/mike/Documents/Utilities/eigen_git/
IDIR= $(foreach dir,$(SRCDIRS),-I$(dir)) $(BOOSTDIR) $(EIGENDIR)

# GNU Compiler
CC=g++ -std=c++17 -fPIC #-D RT_FLOAT -Wfloat-conversion #these commands can be used to check for double literals
LIBS=-lm 
MPFLAGS=-fopenmp
OFLAGS=-O3 -march=native -DNDEBUG 

# Nvidia CUDA Compiler
NCC=nvcc -Xcompiler -fPIC -Xcudafe --display_error_number #--disable-warnings
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

$(OBJDIR)/%.cuda.debug.o: %.cu
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NDBGFLAGS) -dc $< -o $@



py_corona_sim_cpu: $(PYOBJFILES)
	@mkdir -p python/build

# shared library --- python can't find this at runtime
#	$(CC) -Wall -shared -std=c++17 -fPIC \
	-Wl,-soname,libobservationfit.so \
	$(SRCFILES) \
	wrapper.cpp \
	$(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) \
	 -o build/libobservation_fit.so

	rm -f python/build/libobservation_fit.a

	ar rs python/build/libobservation_fit.a $(PYOBJFILES)

	ranlib python/build/libobservation_fit.a

	gfortran -fPIC -Ofast -c -std=legacy\
	  $(SRCDIR)/quemerais_IPH_model/ipbackgroundCFR_fun.f \
	  -o $(OBJDIR)/ipbackgroundCFR_fun.o

	cd python; python setup_corona_sim.py build_ext --inplace

$(OBJDIR)/%.o: %.cpp
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(CC) -c $< $(IDIR) $(LIBS) $(OFLAGS) -o $@

py_corona_sim_gpu: $(PYNOBJFILES)
	@mkdir -p python/build

	$(NCC) -dlink $(PYNOBJFILES) $(NIDIR) $(NLIBS) $(NOFLAGS) \
	-o python/build/observation_fit_wrapper_gpu_device.o

# shared library --- python can't find this at runtime
#	$(CC) -shared -Wl,-soname,libobservationfit.so \
	-o build/libobservation_fit.so \
	build/observation_fit_gpu_host.o \
	build/observation_fit_gpu_device.o -lc

	rm -f python/build/libobservation_fit.a

	ar rs python/build/libobservation_fit.a \
	python/build/observation_fit_wrapper_gpu_device.o \
	$(PYNOBJFILES)

	ranlib python/build/libobservation_fit.a

	gfortran -fPIC -Ofast -c -std=legacy\
	  $(SRCDIR)/quemerais_IPH_model/ipbackgroundCFR_fun.f \
	  -o $(OBJDIR)/ipbackgroundCFR_fun.o

	cd python; python setup_corona_sim.py build_ext --inplace -RT_FLOAT



clean_gpu:
	rm -f generate_source_function_gpu.x
	test -d ./bin && find ./bin/ -type f -name '*cuda*' -delete 2>/dev/null

clean_all:
	rm -f generate_source_function.x
	rm -f generate_source_function_gpu.x
	rm -rf bin
	rm -rf python/build* python/*.cpp
