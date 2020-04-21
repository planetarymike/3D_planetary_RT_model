#ifndef __RT_cuda_compatibility
#define __RT_cuda_compatibility

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 


#endif
