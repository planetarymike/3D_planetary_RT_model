ifneq ($(NO_PARALLEL_COMPILATION), true)
MAKEFLAGS += -j20 #parallel compilation
endif

#files that need compilin'
OBJDIR = ./bin
SRCDIR = ./src
SRCDIRS = ./src ./src/atm ./src/emission ./src/grid
PSRCFILES = $(foreach dir,$(SRCDIRS),$(wildcard $(dir)/*.cpp))

SRCFILES = $(filter-out ./src/observation_fit.cpp, $(PSRCFILES))
OBJFILES    := $(filter %.o, $(SRCFILES:%.cpp=$(OBJDIR)/%.o))
OBJFILESDBG := $(filter %.o, $(SRCFILES:%.cpp=$(OBJDIR)/%.debug.o))


NSRCFILES = $(SRCFILES)
NOBJFILES    := $(filter %.o, $(NSRCFILES:%.cpp=$(OBJDIR)/%.cuda.o))
NOBJFILESDBG := $(filter %.o, $(NSRCFILES:%.cpp=$(OBJDIR)/%.cuda.debug.o))


PYSRCFILES = $(SRCFILES) ./src/observation_fit.cpp $(wildcard $(SRCDIR)/quemerais_IPH_model/*.cpp)
PYOBJFILES  := $(filter %.o, $(PYSRCFILES:%.cpp=$(OBJDIR)/%.o))
PYNOBJFILES := $(filter %.o, $(PYSRCFILES:%.cpp=$(OBJDIR)/%.cuda.o))


# External library dependencies
BOOST_VERSION_NUMBER = 1.83.0
BOOST_VERSION_NUMBER_ = $(subst .,_,$(BOOST_VERSION_NUMBER))
BOOSTDIR = $(shell pwd)/lib/boost_$(BOOST_VERSION_NUMBER_)
EIGEN_VERSION_NUMBER = 3.4.0
EIGENDIR = $(shell pwd)/lib/eigen-$(EIGEN_VERSION_NUMBER)

# include flags for compiler
IDIR= $(foreach dir,$(SRCDIRS),-I$(abspath $(dir))) -I$(BOOSTDIR) -I$(EIGENDIR)

# determine current git hash for versioning
UNCOMMITTED_CHANGES = $(shell git status --porcelain=v1 2>/dev/null | wc -l)
ifeq ($(UNCOMMITTED_CHANGES), 0)
	GIT_HASH = $(shell git rev-parse HEAD)
else
	GIT_HASH+=$(shell git rev-parse HEAD)-dirty
endif

#
# CPU compilation
#

# Select either GNU Compiler or clang based on env CXX
ifneq ($(USE_CLANG),)
	CCOMP = clang++-15
else
	CCOMP = g++-9
	LIBS = -lm -lgomp
	ifeq ($(shell uname),Darwin)
		ifeq ($(shell pkgutil --pkg-info=com.apple.pkg.CLTools_Executables \
			| grep version | cut -d' ' -f2 | cut -d'.' -f1),15)
			LIBS += -Wl,-ld_classic # Xcode 15 on mac needs this linker command to compile correctly
		endif
	endif
endif
CC = $(CCOMP) -std=c++17 -fPIC #-D RT_FLOAT -Wfloat-conversion # these commands can be used to check for double literals
MPFLAGS = -fopenmp
OFLAGS = -O3 -DNDEBUG -g #-march=native


#
# GPU compilation
#
NVCC := $(shell command -v nvcc 2> /dev/null)
ifneq ($(NVCC),) # if cuda toolkit is installed
	CUDA_VERSION=$(shell nvcc --version | grep -o 'release.*' | cut -f2 -d' ' | cut -f1 -d,)
	CUDA_DEVICE_CODE=$(shell $$CUDA_HOME/extras/demo_suite/deviceQuery \
			   | grep -o 'CUDA Capability Major/Minor version number:.*' \
                           | cut -f2 -d ':' | sed -r 's/\s+//g' | sed 's/\.//')

	CUDA_SAMPLES_DIR = $(shell pwd)/lib/cuda-samples-$(CUDA_VERSION)
	NIDIR=$(IDIR) -I$(CUDA_HOME)/lib64/ -I$(CUDA_SAMPLES_DIR)/Common/

	NFLAGS =-D RT_FLOAT # 32-bit calculation
	NFLAGS += -D EIGEN_NO_CUDA # disable Eigen on device
	NFLAGS += -D BOOST_NO_CUDA # disable Boost on device (see BOOST_DIR recipe)

	NLIBS=-lcudart -lcusolver -lcublas

	NOFLAGS=-O3 -DNDEBUG  #--maxrregcount 43
	NDBGFLAGS=-g

	ifneq ($(USE_CLANG),) # use clang compiler
		NCC=$(CLANGPP) -std=c++17 -fPIC --offload-new-driver
		NFLAGS +=-x cuda
		NOFLAGS += -ffast-math

		NLIBS += -ldl -lrt -pthread -L$(CUDA_HOME)/lib64
		NOFLAGS += --offload-arch=sm_$(CUDA_DEVICE_CODE)
		NDBGFLAGS += --offload-arch=sm_$(CUDA_DEVICE_CODE)$(ARCH_SM) -O0
	else # use default nvcc compiler
		NCC=nvcc -Xcompiler -fPIC -Xcudafe --display_error_number #--disable-warnings
		NFLAGS +=-x cu
		NLIBS +=-lm
		NOFLAGS += -lineinfo --use_fast_math

		ifeq ($(shell echo "$(CUDA_VERSION) > 11" |bc -l),1)
			CUDA_SM_TYPE=lto
			DLTO=-dlto
		else
			CUDA_SM_TYPE=sm
		endif

		ARCH=--generate-code arch=compute_$(CUDA_DEVICE_CODE),code=$(CUDA_SM_TYPE)_$(CUDA_DEVICE_CODE) # local machine
		ARCH +=--generate-code arch=compute_70,code=$(CUDA_SM_TYPE)_70 # AWS P3 node

		ARCH_SM=--generate-code arch=compute_$(CUDA_DEVICE_CODE),code=sm_$(CUDA_DEVICE_CODE) # local machine
		ARCH_SM +=--generate-code arch=compute_70,code=sm_70 # AWS P3 node

		NOFLAGS +=$(ARCH)
		NDBGFLAGS += -G $(ARCH_SM) -Xcompiler -O0 -Xptxas -O0
#                           ^^ this -G sometimes changes the behavior of the code??
	endif
endif

#
# recipes
#

all:
	$(error "please specify recipe")

generate_source_function: $(EIGENDIR) $(BOOSTDIR)
	@$(CC) generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) -o generate_source_function.x

generate_source_function_profile: $(EIGENDIR) $(BOOSTDIR)
	@$(CC) generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) $(OFLAGS) -g -o generate_source_function.x

# generate_source_function_debug_warn:
# 	$(CC) generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) -v -O0 -g -Wall -Wextra -Wno-unknown-pragmas -o generate_source_function.x

generate_source_function_debug_warn: $(OBJFILESDBG) $(EIGENDIR) $(BOOSTDIR)
	@echo "compiling generate_source_function.cpp..."
	@$(CC) $(IDIR) $(LIBS) -O0 -g -Wall -Wextra -Wno-unknown-pragmas -c generate_source_function.cpp -o bin/generate_source_function.debug.o
	@echo "linking ..."
	@$(CC) $(OBJFILESDBG) bin/generate_source_function.debug.o $(IDIR) $(LIBS) -O0 -g -Wall -Wextra -Wno-unknown-pragmas -o generate_source_function.x


generate_source_function_debug_warn_float: $(EIGENDIR) $(BOOSTDIR)
	@$(CC) -D RT_FLOAT generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) -O0 -g -Wall -Wextra -Wno-unknown-pragmas -o generate_source_function.x

generate_source_function_debug_warn_MP: $(EIGENDIR) $(BOOSTDIR)
	@$(CC) generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) $(MPFLAGS) -O0 -g -Wall -Wextra -o generate_source_function.x

generate_source_function_float: $(EIGENDIR) $(BOOSTDIR)
	@$(CC) -D RT_FLOAT generate_source_function.cpp $(SRCFILES) $(IDIR) $(LIBS) $(MPFLAGS) $(OFLAGS) -o generate_source_function.x



generate_source_function_gpu: $(NOBJFILES) $(EIGENDIR) $(BOOSTDIR) $(CUDA_SAMPLES_DIR) cuda_installed
	@echo "compiling generate_source_function.cpp..."
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NOFLAGS) -dc generate_source_function.cpp -o bin/generate_source_function.cuda.o
	@echo "linking..."
ifneq ($(DLTO),)
	@echo "Using CUDA link time optimization"
endif
	@$(NCC) $(NOBJFILES) bin/generate_source_function.cuda.o $(NIDIR) $(NLIBS) $(NOBASEFLAGS) $(ARCH_SM) $(DLTO) -o generate_source_function_gpu.x

$(OBJDIR)/%.cuda.o: %.cpp $(EIGENDIR) $(BOOSTDIR) $(CUDA_SAMPLES_DIR) cuda_installed
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NOFLAGS) -dc $< -o $@


generate_source_function_gpu_monolithic: $(EIGENDIR) $(BOOSTDIR) $(CUDA_SAMPLES_DIR) cuda_installed
	$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NOFLAGS) $(NSRCFILES) generate_source_function.cpp -o generate_source_function_gpu.x

generate_source_function_gpu_debug: $(NOBJFILESDBG) $(EIGENDIR) $(BOOSTDIR) $(CUDA_SAMPLES_DIR) cuda_installed
	@echo "compiling generate_source_function.cpp..."
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NDBGFLAGS) -dc generate_source_function.cpp -o bin/generate_source_function.cuda.debug.o
	@echo "linking ..."
	@$(NCC) $(NOBJFILESDBG) bin/generate_source_function.cuda.debug.o $(NIDIR) $(NLIBS) $(NDBGFLAGS) -o generate_source_function_gpu.x

$(OBJDIR)/%.cuda.debug.o: %.cpp $(EIGENDIR) $(BOOSTDIR) $(CUDA_SAMPLES_DIR) cuda_installed
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NDBGFLAGS) -dc $< -o $@



py_corona_sim_cpu: $(EIGENDIR) $(BOOSTDIR)
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
	export CC='$(CCOMP)'; \
	export CXX='$(CCOMP)'; \
	export CPP_GIT_HASH='$(GIT_HASH)'; \
	python setup_corona_sim.py build_ext --inplace

$(OBJDIR)/%.o: %.cpp $(EIGENDIR) $(BOOSTDIR)
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(CC) -c $< $(IDIR) $(LIBS) $(OFLAGS) -o $@

$(OBJDIR)/%.debug.o: %.cpp $(EIGENDIR) $(BOOSTDIR)
	@echo "compiling $<..."
	@mkdir -p '$(@D)'
	@$(CC) -c $< $(IDIR) $(LIBS)  -O0 -g -Wall -Wextra -Wno-unknown-pragmas -o $@

py_corona_sim_gpu: $(EIGENDIR) $(BOOSTDIR) $(CUDA_SAMPLES_DIR) cuda_installed
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
	export CPP_GIT_HASH='$(GIT_HASH)'; \
	python setup_corona_sim.py build_ext --inplace -RT_FLOAT -v

py_corona_sim_all: $(EIGENDIR) $(BOOSTDIR) $(CUDA_SAMPLES_DIR) cuda_installed
	make py_corona_sim_cpu && make py_corona_sim_gpu

observation_fit_cpu_test: $(EIGENDIR) $(BOOSTDIR)
	@echo "compiling observation_fit.cpp..."
	$(CC) $(IDIR) $(LIBS) -DNDEBUG -O0 -g -c src/observation_fit.cpp -o bin/src/observation_fit.debug.o
	@echo "compiling obs_fit_test.cpp..."
	@$(CC) $(IDIR) $(LIBS) -O0 -g -c python/test/obs_fit_test.cpp -o bin/obs_fit_test.debug.o

	@echo "compiling ipbackgroundCFR_fun.f..."
	@gfortran -fPIC -Ofast -c -std=legacy \
	$(SRCDIR)/quemerais_IPH_model/ipbackgroundCFR_fun.f \
	-o $(OBJDIR)/ipbackgroundCFR_fun.o
	@$(CC) $(IDIR) $(LIBS) -O0 -g -c $(SRCDIR)/quemerais_IPH_model/*.cpp -o bin/quemerais_IPH_model.debug.o


	@echo "linking ..."
	@$(CC) $(SRCFILES) \
	bin/src/observation_fit.debug.o \
	bin/obs_fit_test.debug.o \
	bin/quemerais_IPH_model.debug.o \
	$(OBJDIR)/ipbackgroundCFR_fun.o -lgfortran \
	$(IDIR) $(LIBS)  -O0 -g -o python/test/obs_fit_test.x

observation_fit_gpu_test: $(NOBJFILESDBG) $(EIGENDIR) $(BOOSTDIR) $(CUDA_SAMPLES_DIR) cuda_installed
	@echo "compiling observation_fit.cpp..."
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NDBGFLAGS) -dc src/observation_fit.cpp -o bin/src/observation_fit.cuda.debug.o
	@echo "compiling obs_fit_test.cpp..."
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NDBGFLAGS) -dc python/test/obs_fit_test.cpp -o bin/obs_fit_test.cuda.debug.o

	@echo "compiling ipbackgroundCFR_fun.f..."
	@gfortran -fPIC -Ofast -c -std=legacy \
	$(SRCDIR)/quemerais_IPH_model/ipbackgroundCFR_fun.f \
	-o $(OBJDIR)/ipbackgroundCFR_fun.o
	@$(NCC) $(NFLAGS) $(NIDIR) $(NLIBS) $(NDBGFLAGS) -dc $(SRCDIR)/quemerais_IPH_model/*.cpp -o bin/quemerais_IPH_model.cuda.debug.o


	@echo "linking ..."
	@$(NCC) $(NOBJFILESDBG) \
	bin/src/observation_fit.cuda.debug.o \
	bin/obs_fit_test.cuda.debug.o \
	bin/quemerais_IPH_model.cuda.debug.o \
	$(OBJDIR)/ipbackgroundCFR_fun.o -lgfortran \
	$(NIDIR) $(NLIBS) $(NDBGFLAGS) -o python/test/obs_fit_test.x

$(EIGENDIR):
	@echo "Downloading Eigen library..."
	@mkdir -p lib
# remove old versions
	@cd lib && find . -name "eigen-*" -type d ! -name "eigen-$(EIGEN_VERSION_NUMBER)" -exec rm -rf {} +
	@cd lib && rm -f *.bz2
# download and extract
	@cd lib && wget -nv --no-hsts https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION_NUMBER)/eigen-$(EIGEN_VERSION_NUMBER).tar.bz2
	@cd lib && tar -xf eigen-$(EIGEN_VERSION_NUMBER).tar.bz2
	@cd lib && rm eigen-$(EIGEN_VERSION_NUMBER).tar.bz2
# fix obsolete usage of diag_suppress to eliminate scores of nvcc warnings
	@sed -i'.original' -e 's/diag_suppress/nv_diag_suppress/' lib/eigen-3.4.0/Eigen/src/Core/util/DisableStupidWarnings.h
# fix improper identification of cuda in an Eigen source file
	@sed -i'.original' -e 's/#if defined(__clang__) && defined(__CUDA__)/#if defined(EIGEN_HAS_GPU_FP16) || defined(EIGEN_HAS_ARM64_FP16_SCALAR_ARITHMETIC)/' lib/eigen-3.4.0/Eigen/src/Core/arch/Default/Half.h
	@echo "... done downloading Eigen."

$(BOOSTDIR):
	@echo "Downloading Boost library..."
	@mkdir -p lib
# remove old versions
	@cd lib && find . -name "boost_*" -type d ! -name "boost_$(BOOST_VERSION_NUMBER_)" -exec rm -rf {} +
	@cd lib && rm -f *.bz2
# download and extract
	@cd lib && wget -nv --no-hsts https://boostorg.jfrog.io/artifactory/main/release/$(BOOST_VERSION_NUMBER)/source/boost_$(BOOST_VERSION_NUMBER_).tar.bz2
	@cd lib && tar -xf boost_$(BOOST_VERSION_NUMBER_).tar.bz2
	@cd lib && rm boost_$(BOOST_VERSION_NUMBER_).tar.bz2
# implement BOOST_NO_CUDA flag
	@sed -i'.original' -e 's/#define BOOST_GPU_ENABLED __host__ __device__/#ifndef BOOST_NO_CUDA\n#define BOOST_GPU_ENABLED __host__ __device__\n#endif/' lib/boost_$(BOOST_VERSION_NUMBER_)/boost/config/compiler/nvcc.hpp
	@echo "... done downloading Boost."

$(CUDA_SAMPLES_DIR):
	@echo "Downloading cuda samples library..."
	@mkdir -p lib
# remove old versions
	@cd lib && find . -name "cuda-samples*" -type d ! -name "cuda_samples_$(CUDA_VERSION)" -exec rm -rf {} +
	@cd lib && rm -f *.tar.gz
# download and extract
	@cd lib && wget -nv --no-hsts https://github.com/NVIDIA/cuda-samples/archive/refs/tags/v$(CUDA_VERSION).tar.gz
	@cd lib && tar -xf v$(CUDA_VERSION).tar.gz
	@cd lib && rm v$(CUDA_VERSION).tar.gz
	@echo "... done downloading cuda samples library."

cuda_installed:
	@if [ -z "$(CUDA_VERSION)" ] ; \
	then echo "The CUDA nvcc compiler was not found, please install CUDA to compile GPU code.\n\n"  \
	"If the CUDA toolkit is installed, you may need to add the following lines to your .bashrc:\n" \
	"  export CUDA_HOME=/usr/local/cuda/\n" \
	"  export PATH=""$$""CUDA_HOME/bin""$$""{PATH:+:""$$""{PATH}}""$$\n" \
	"  export LD_LIBRARY_PATH=""$$""CUDA_HOME/lib64""$$""{LD_LIBRARY_PATH:+:""$$""{LD_LIBRARY_PATH}}\n" \
	; false ; fi

clean_gpu:
	rm -f generate_source_function_gpu.x
	test -d ./bin && find ./bin/ -type f -name '*cuda*' -delete 2>/dev/null

clean_all:
	rm -f generate_source_function.x
	rm -f generate_source_function_gpu.x
	rm -rf bin
	rm -rf python/build* python/*.cpp #python/*.so
