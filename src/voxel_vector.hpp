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

  voxel_vector() : VectorX() { }

  ~voxel_vector() {
#ifdef __CUDACC__
    if(d_vec!=NULL)
      checkCudaErrors(cudaFree(d_vec));
#endif
  }

  
  voxel_vector(const voxel_vector &copy) 
    : VectorX(copy)
  {
    assert(n_voxels == VectorX::size());
    vec = VectorX::data();
  }
  template <typename T>
  voxel_vector& operator=(T &rhs) {
    assert(n_voxels == VectorX::size());
    VectorX::operator=(rhs);
    vec = VectorX::data();
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
#endif
};

template<int N_VOXELS> //template prevents allocation errors with CUDA pointers
class voxel_matrix : public MatrixX {
public:
  static const int n_voxels = N_VOXELS;
  Real* mat;
  Real* d_mat = NULL;//pointer to device memory for CUDA

  voxel_matrix() : MatrixX() { }

  ~voxel_matrix() {
#ifdef __CUDACC__
    if(d_mat!=NULL)
      checkCudaErrors(cudaFree(d_mat));
#endif
  }
  
  
  voxel_matrix(const voxel_matrix &copy) 
    : MatrixX(copy)
  {
    assert(MatrixX::rows() == MatrixX::cols()
	   && "voxel_matrix must be square");
    assert(n_voxels == MatrixX::rows());
    mat = MatrixX::data();
  }
  template <typename T>
  voxel_matrix& operator=(T &rhs) {
    MatrixX::operator=(rhs);
    assert(MatrixX::rows() == MatrixX::cols()
	   && "voxel_matrix must be square");
    assert(n_voxels == MatrixX::rows());
    mat = MatrixX::data();
    return *this;
  }

  void resize() {
    MatrixX::resize(n_voxels,n_voxels);
    mat = MatrixX::data();
  }

  //overload [] to access coefficients of mat
  CUDA_CALLABLE_MEMBER
  Real & operator()(const int n, const int m) {
#ifdef EIGEN_ROWMAJOR
    const int i = n*n_voxels + m;//row major
#else
    const int i = m*n_voxels + n;//column major
#endif
    return mat[i];
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
