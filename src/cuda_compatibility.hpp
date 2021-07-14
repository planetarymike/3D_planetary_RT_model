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


// these routines automate creation of static const array members

#define ARRAY_NAME(name) name ## _array
#define ARRAY_NAME_GPU(name) name ## _array_gpu

#define DECLARE_STATIC_ARRAY_HPP(type, n_elements, name, ...)	\
  static constexpr type ARRAY_NAME(name)[n_elements] = __VA_ARGS__;	\
  extern CUDA_CONST type ARRAY_NAME_GPU(name)[n_elements];

#define DECLARE_STATIC_ARRAY_CPP(type, n_elements, name, ...)	\
  CUDA_CONST       type ARRAY_NAME_GPU(name)[n_elements] = __VA_ARGS__;

#ifndef __CUDA_ARCH__
#define CUDA_STATIC_ARRAY_MEMBER(nspace, type, n_elements, name) \
  static constexpr type name(const int i) {		         \
    assert(0<=i && i<n_elements && "index must be in bounds");  \
    return nspace::ARRAY_NAME(name)[i];				 \
  }
#else
#define CUDA_STATIC_ARRAY_MEMBER(nspace, type, n_elements, name) \
  __device__ static constexpr type name(const int i){		 \
    assert(0<=i && i<n_elements && "index must be in bounds");  \
    return nspace::ARRAY_NAME_GPU(name)[i];			 \
  }
#endif

#endif
