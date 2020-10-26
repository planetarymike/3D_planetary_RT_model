//voxel_vector.hpp ---

//provides an Eigen VectorX and MatrixX with a pointer to the internal
//data to make it easier to write CPU/GPU interoperable code and copy
//to/from the GPU


#ifndef __VOXEL_VECTOR_H
#define __VOXEL_VECTOR_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#ifdef __CUDACC__
#include "helper_cuda.h"
#endif

template<int N_VOXELS> //template prevents allocation errors with CUDA pointers
class voxel_vector : public VectorX {
public:
  static const int n_voxels = N_VOXELS;
  Real* vec;
  Real* d_vec = NULL;//pointer to device memory for CUDA

  voxel_vector() : VectorX() {
    resize();
  }

  ~voxel_vector() {
    free_d_vec();
  }

  
  voxel_vector(const voxel_vector &copy) 
    : VectorX(copy)
  {
    assert(n_voxels == VectorX::size());
    vec = VectorX::data();
  }
  voxel_vector& operator=(const VectorX &rhs) {
    assert(n_voxels == VectorX::size());
    VectorX::operator=(rhs);
    vec = VectorX::data();
    return *this;
  }
  CUDA_CALLABLE_MEMBER
  voxel_vector& operator=(const voxel_vector<N_VOXELS> &rhs) {
#ifdef __CUDA_ARCH__
    for (int i=0;i<n_voxels;i++)
      vec[i] = rhs.vec[i];
#else
    assert(n_voxels == VectorX::size());
    VectorX::operator=(rhs.VectorX);
    vec = VectorX::data();
#endif
    return *this;
  }

  void resize() {
    VectorX::resize(n_voxels);
    vec = VectorX::data();
  }

  CUDA_CALLABLE_MEMBER
  Real & operator[](const int n) {
    return vec[n];
  }
  CUDA_CALLABLE_MEMBER
  Real operator[](const int n) const {
    return vec[n];
  }
  CUDA_CALLABLE_MEMBER
  Real & operator()(const int n) {
    return vec[n];
  }
  CUDA_CALLABLE_MEMBER
  Real operator()(const int n) const {
    return vec[n];
  }

  void free_d_vec() {
#ifdef __CUDACC__
    if(d_vec!=NULL) {
      checkCudaErrors(cudaFree(d_vec));
      d_vec = NULL;
    }
#endif
  }
  
#ifdef __CUDACC__  
  void to_device(bool transfer = true) {
    if (d_vec == NULL)
      //allocate the host's d_vec to point at device memory
      checkCudaErrors(
		      cudaMalloc((void **) &d_vec,
				 n_voxels*sizeof(Real))
		      );
    if (transfer)
      //copy from host to host's device pointer
      checkCudaErrors(
		      cudaMemcpy(d_vec,
				 vec,
				 n_voxels*sizeof(Real),
				 cudaMemcpyHostToDevice)
		      );
  }
  void to_host() {
    checkCudaErrors(
		    cudaMemcpy(vec,
			       d_vec,
			       n_voxels*sizeof(Real),
			       cudaMemcpyDeviceToHost)
		    );
  }

  __device__
  void to_shared() {
    //declare shared array
    __shared__ Real s_vec[n_voxels];

    //move to shared array using all threads in this block
    for (int i = threadIdx.x; i < n_voxels; i+=blockDim.x)
      s_vec[i] = vec[i];
    __syncthreads();

    d_vec = vec;
    vec = s_vec;    
  }
  __device__
  void from_shared() {
    vec = d_vec;
    d_vec = NULL;
    //shared memory is automatically deallocated
  }

#endif

};

template<int N_VOXELS> //template prevents allocation errors with CUDA pointers
class voxel_matrix : public MatrixX {
public:
  static const int n_voxels = N_VOXELS;
  Real* mat;
  Real* d_mat = NULL;//pointer to device memory for CUDA

  voxel_matrix() : MatrixX() {
    resize();
  }

  ~voxel_matrix() {
    free_dmat();
  }
  
  
  voxel_matrix(const voxel_matrix &copy) 
    : MatrixX(copy)
  {
    assert(MatrixX::rows() == MatrixX::cols()
	   && "voxel_matrix must be square");
    assert(n_voxels == MatrixX::rows());
    mat = MatrixX::data();
  }

  voxel_matrix& operator=(const MatrixX &rhs) {
    assert(n_voxels == MatrixX::size());
    MatrixX::operator=(rhs);
    mat = MatrixX::data();
    return *this;
  }
  CUDA_CALLABLE_MEMBER
  voxel_matrix& operator=(const voxel_matrix<N_VOXELS> &rhs) {
#ifdef __CUDA_ARCH__
    for (int i=0;i<n_voxels*n_voxels;i++)
      mat[i] = rhs.mat[i];
#else
    assert(n_voxels == MatrixX::size());
    MatrixX::operator=(rhs.MatrixX);
    mat = MatrixX::data();
#endif
    return *this;
  }

  void resize() {
    MatrixX::resize(n_voxels,n_voxels);
    mat = MatrixX::data();
  }

  //overload () to access coefficients of mat
  CUDA_CALLABLE_MEMBER
  Real & operator()(const int n, const int m) {
#ifdef EIGEN_ROWMAJOR
    const int i = n*n_voxels + m;//row major
#else
    const int i = m*n_voxels + n;//column major
#endif
    return mat[i];
  }
  CUDA_CALLABLE_MEMBER
  const Real & operator()(const int n, const int m) const {
#ifdef EIGEN_ROWMAJOR
    const int i = n*n_voxels + m;//row major
#else
    const int i = m*n_voxels + n;//column major
#endif
    return mat[i];
  }

  void free_dmat() {
#ifdef __CUDACC__
    if(d_mat!=NULL) {
      checkCudaErrors(cudaFree(d_mat));
      d_mat = NULL;
    }
#endif
  }

#ifdef __CUDACC__
  void to_device(bool transfer = true) {
    if (d_mat == NULL)
      //allocate the host's d_vec to point at device memory
      checkCudaErrors(
		      cudaMalloc((void **) &d_mat,
				 n_voxels*n_voxels*sizeof(Real))
		      );
    if (transfer)
      //copy from host to host's device pointer
      checkCudaErrors(
		      cudaMemcpy(d_mat,
				 mat,
				 n_voxels*n_voxels*sizeof(Real),
				 cudaMemcpyHostToDevice)
		      );
  }
  void to_host() {
    checkCudaErrors(
		    cudaMemcpy(mat,
			       d_mat,
			       n_voxels*n_voxels*sizeof(Real),
			       cudaMemcpyDeviceToHost)
		    );
  }
#endif
};

#endif
