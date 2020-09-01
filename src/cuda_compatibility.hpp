#ifndef __RT_cuda_compatibility
#define __RT_cuda_compatibility

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 

#ifdef __CUDACC__
#define CUDA_CONST __constant__
#else
#define CUDA_CONST
#endif 

#ifdef __CUDACC__
#include "helper_cuda.h"
#endif

#endif
