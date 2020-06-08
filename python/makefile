MAKEFLAGS += -j20 #parallel compilation

#files that need compilin'
SRCDIR = ../src
SRCFILES = $(wildcard $(SRCDIR)/atm/*.cpp) $(wildcard $(SRCDIR)/*.cpp) 
OBJDIR = ../bin/src
OBJFILES  = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCFILES))
NOBJFILES = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.cuda.o,$(SRCFILES))

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
NIDIR=$(IDIR) \
      -L$(CUDA_HOME)/lib64/ \
      -I$(CUDA_HOME)/samples/common/inc/
NLIBS=-lm -lcudart -lcusolver -lcublas
NOFLAGS= -O3 -DNDEBUG #-lineinfo

py_corona_sim_cpu: $(OBJFILES)
	@mkdir -p build

# shared library --- python can't find this at runtime
#	$(CC) -Wall -shared -std=c++17 -fPIC \
	-Wl,-soname,libobservationfit.so \
	$(SRCFILES) \
	wrapper.cpp \
	$(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) \
	 -o build/libobservation_fit.so

	rm -f build/libobservation_fit.a

	ar rs build/libobservation_fit.a $(OBJFILES)
	ranlib build/libobservation_fit.a

	python setup_corona_sim.py build_ext --inplace

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(CC) -c $< -std=c++17 -fPIC $(IDIR) $(LIBS) $(OFLAGS) -o $@

py_corona_sim_gpu: $(NOBJFILES)
	@mkdir -p build

	$(NCC) -Xcompiler -fPIC \
	-D RT_FLOAT \
	$(IDIR) $(NIDIR) $(NLIBS) $(NOFLAGS) \
	-dlink $(NOBJFILES) \
	-o build/observation_fit_wrapper_gpu_device.o

# shared library --- python can't find this at runtime
#	$(CC) -shared -Wl,-soname,libobservationfit.so \
	-o build/libobservation_fit.so \
	build/observation_fit_gpu_host.o \
	build/observation_fit_gpu_device.o -lc

	rm -f build/libobservation_fit.a

	ar rs build/libobservation_fit.a \
	build/observation_fit_wrapper_gpu_device.o \
	$(NOBJFILES)

	ranlib build/libobservation_fit.a

	python setup_corona_sim.py build_ext --inplace -RT_FLOAT

$(OBJDIR)/%.cuda.o: $(SRCDIR)/%.cpp
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) -x cu -D RT_FLOAT -Xcompiler -fPIC $(NIDIR) $(NLIBS) $(NOFLAGS) -dc $< -o $@

clean_all:
	rm -rf *.cpp build* *.so ../bin/*

clean_python:
	rm -rf build* *.so *.cpp
