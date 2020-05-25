BINDIR=-I../bin/
BOOSTDIR=-I/home/mike/Documents/Utilities/boost_1_72_0/
EIGENDIR=-I/home/mike/Documents/Utilities/eigen-3.3.7/
IDIR=$(BINDIR) $(BOOSTDIR) $(EIGENDIR)

# GNU Compiler
CC=g++
LIBS=-lm 
MPFLAGS=-fopenmp
#OFLAGS=-Ofast -march=native -DNDEBUG
OFLAGS=-O3 -march=native -DNDEBUG

# Nvidia CUDA Compiler
NCC=nvcc --disable-warnings
NLIBS=-lm
NOFLAGS= -O3 -DNDEBUG #-lineinfo
CUDA_HELP_IDIR=-I/home/mike/Documents/Utilities/cuda-samples/Common/

py_corona_sim_cpu:
	python setup_corona_sim.py build_ext --inplace

py_corona_sim_gpu:
	$(NCC) -x cu -Xcompiler -fPIC \
	$(COMPILER_OPT) $(IDIR) $(CUDA_HELP_IDIR) $(NLIBS) $(NOFLAGS) \
	-dc observation_fit_gpu.h -o build/observation_fit_gpu_host.o
	$(NCC) -Xcompiler -fPIC \
	$(COMPILER_OPT) $(IDIR) $(CUDA_HELP_IDIR) $(NLIBS) $(NOFLAGS) \
	-dlink observation_fit_gpu_host.o -o build/observation_fit_gpu_device.o
	python setup_corona_sim.py build_ext --inplace
clean:
	rm -rf build* *.cpp *.so *.o
