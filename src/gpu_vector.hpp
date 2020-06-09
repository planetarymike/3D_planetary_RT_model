//gpu_vector.hpp --- vector to communicate obs stuff to GPU

#ifndef __gpu_vec
#define __gpu_vec

#include "cuda_compatibility.hpp"
#ifdef __CUDACC__
#include <helper_cuda.h>
#endif

template <typename T>
class gpu_vector {
public:
  int n;
  T *v = NULL;
  T *d_v = NULL;

  gpu_vector() : n(0) { }
  ~gpu_vector() {
    if (v)
      delete [] v;
#ifdef __CUDACC__
    if(d_v)
      checkCudaErrors(cudaFree(d_v));
#endif
  }
  
  gpu_vector(const gpu_vector &copy) {
    n = copy.n;
    if (copy.v) {
      resize(n);
      for (int i = 0; i < n; i++)
	v[i] = copy.v[i];
    }
    if (copy.d_v)
      to_device();
  }
  gpu_vector& operator=(gpu_vector &rhs) {
    if(this == &rhs) return *this;
    n = rhs.n;
    if (rhs.v) {
      resize(n);
      for (int i = 0; i < n; i++)
	v[i] = rhs.v[i];
    }
    if (rhs.d_v)
      to_device();
    return *this;
  }

  void resize(const int nn) {
    n = nn;
    if (v)
      delete [] v;
    v = new T[n];
  }
  void resize(const int nn, const T& val) {
    resize(nn);
    for (int i=0;i<n;i++)
      v[i] = val;
  }

  int size() const {
    return n;
  }

  CUDA_CALLABLE_MEMBER
  T & operator[](const int n) {
    return v[n];
  }
  CUDA_CALLABLE_MEMBER
  T operator[](const int n) const {
    return v[n];
  }
  CUDA_CALLABLE_MEMBER
  T & operator()(const int n) {
    return v[n];
  }
  CUDA_CALLABLE_MEMBER
  T operator()(const int n) const {
    return v[n];
  }

  void to_device(bool transfer = true) {
#ifdef __CUDACC__  
    if (d_v != NULL)
      checkCudaErrors(cudaFree(d_v));
    //allocate the host's d_v to point at device memory
    checkCudaErrors(
		    cudaMalloc((void **) &d_v,
			       n*sizeof(T))
		    );
    if (transfer)
      //copy from host to host's device pointer
      checkCudaErrors(
		      cudaMemcpy(d_v,
				 v,
				 n*sizeof(T),
				 cudaMemcpyHostToDevice)
		      );
#endif
  }
  void to_host() {
#ifdef __CUDACC__  
    checkCudaErrors(
		    cudaMemcpy(v,
			       d_v,
			       n*sizeof(T),
			       cudaMemcpyDeviceToHost)
		    );
#endif
  }
};


#endif
