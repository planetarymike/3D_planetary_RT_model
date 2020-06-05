//voxel_vector.hpp ---

//provides an Eigen VectorX and MatrixX with a pointer to the internal
//data to make it easier to write CPU/GPU interoperable code and copy
//to/from the GPU


#ifndef __VOXEL_VECTOR_H
#define __VOXEL_VECTOR_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"

class voxel_vector : public VectorX {
public:
  int n_voxels;
  Real* vec;
  Real* d_vec;//pointer to device memory for CUDA

  voxel_vector() : VectorX() { }
  
  voxel_vector(const voxel_vector &copy) 
    : VectorX(copy)
  {
    n_voxels = VectorX::size();
    vec = VectorX::data();
  }
  template <typename T>
  voxel_vector& operator=(T &rhs) {
    VectorX::operator=(rhs);
    n_voxels = VectorX::size();
    vec = VectorX::data();
    return *this;
  }

  void resize(int n_voxelss) {
    n_voxels = n_voxelss;
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
};
class voxel_matrix : public MatrixX {
public:
  int n_voxels;
  Real* mat;
  Real* d_mat;//pointer to device memory for CUDA

  voxel_matrix() : MatrixX() { }
  
  voxel_matrix(const voxel_matrix &copy) 
    : MatrixX(copy)
  {
    assert(MatrixX::rows() == MatrixX::cols()
	   && "voxel_matrix must be square");
    n_voxels = MatrixX::rows();
    mat = MatrixX::data();
  }
  template <typename T>
  voxel_matrix& operator=(T &rhs) {
    MatrixX::operator=(rhs);
    assert(MatrixX::rows() == MatrixX::cols()
	   && "voxel_matrix must be square");
    n_voxels = MatrixX::rows();
    mat = MatrixX::data();
    return *this;
  }

  void resize(int n_voxelss) {
    n_voxels = n_voxelss;
    MatrixX::resize(n_voxels,n_voxels);
    mat = MatrixX::data();
  }

  //overload [] to access coefficients of mat
  CUDA_CALLABLE_MEMBER
  Real & operator()(const int n, const int m) {
    const int i = n*n_voxels + m;//row major
    return mat[i];
  }
  // CUDA_CALLABLE_MEMBER
  // Real operator()(const int n, const int m) const {
  //   const int i = n*n_voxels + m;//row major
  //   return mat[i];
  // }
};

#endif
