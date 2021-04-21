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

template<int N_VOXELS, int N_STATES_PER_VOXEL>
struct voxel_array {
  //device-side array storage used in los_tracker objects
  
  static const int n_voxels = N_VOXELS;
  static const int n_states = N_STATES_PER_VOXEL;
  static const int n_elements = N_VOXELS*N_STATES_PER_VOXEL;

  //#ifndef __CUDA_ARCH__ //shared memory version
  Real vec[N_VOXELS*N_STATES_PER_VOXEL];
  // #else
  //   Real *vec; // needs to be set by kernel to the address of an
  // 	     // allocated pointer in shared memory
  // #endif

  // element fetch
  CUDA_CALLABLE_MEMBER
  static int get_element_num(const int n_voxel, const int n_state) {
    return n_voxel*n_states + n_state;
  }
  
  CUDA_CALLABLE_MEMBER
  Real & operator()(const int n_voxel, const int n_state=0) {
    return vec[get_element_num(n_voxel, n_state)];
  }
  CUDA_CALLABLE_MEMBER
  Real operator()(const int n_voxel, const int n_state = 0) const {
    return vec[get_element_num(n_voxel, n_state)];
  }

  CUDA_CALLABLE_MEMBER
  Real & operator[](const int i_el) {
    return vec[i_el];
  }
  CUDA_CALLABLE_MEMBER
  Real operator()(const int i_el) const {
    return vec[i_el];
  }
};

template <int N_VOXELS, int N_STATES_PER_VOXEL>
class voxel_vector {
  //host-side array storage with device mirror,
  // used in emission_voxels and derived classes.
  // supports interaction with Eigen
public:
  static const int n_voxels = N_VOXELS;
  static const int n_states = N_STATES_PER_VOXEL;
  static const int n_elements = N_VOXELS*N_STATES_PER_VOXEL;

  VectorX *eigen_vec;
  Real* vec;
  Real* d_vec = NULL;//pointer to device memory for CUDA
  
  CUDA_CALLABLE_MEMBER
  voxel_vector() {
#ifndef __CUDA_ARCH__
    eigen_vec = new VectorX;
    resize();
#else
    vec = new Real[n_elements];
#endif
  }

  CUDA_CALLABLE_MEMBER
  ~voxel_vector() {
#ifndef __CUDA_ARCH__
    delete eigen_vec;
    free_d_vec();
#else
    delete [] vec;
#endif
  }
  
  voxel_vector(const voxel_vector &copy) 
  {
    assert(n_elements == copy.eigen_vec->size());
    *eigen_vec = *copy.eigen_vec;
    vec = eigen_vec->data();
  }
  voxel_vector& operator=(const VectorX &rhs) {
    assert(n_elements == rhs.size());
    *eigen_vec = rhs;
    vec = eigen_vec->data();
    return *this;
  }
  voxel_vector& operator=(const voxel_vector<N_VOXELS, N_STATES_PER_VOXEL> &rhs) {
#ifdef __CUDA_ARCH__
    for (int i=0;i<n_elements;i++)
      vec[i] = rhs.vec[i];
#else
    *eigen_vec = *rhs.eigen_vec;
    vec = eigen_vec->data();
#endif
    return *this;
  }
  operator VectorX() const {
    return eigen();
  }

  void resize() {
    eigen_vec->resize(n_elements);
    vec = eigen_vec->data();
  }

  VectorX & eigen() {
    return *eigen_vec;
  }
  const VectorX eigen() const {
    return *eigen_vec;
  }

  // element fetch
  CUDA_CALLABLE_MEMBER
  int get_element_num(const int n_voxel, const int n_state) const {
    return voxel_array<N_VOXELS, N_STATES_PER_VOXEL>::get_element_num(n_voxel, n_state);
  }

  CUDA_CALLABLE_MEMBER
  Real & operator()(const int n_voxel, const int n_state=0) {
    return vec[get_element_num(n_voxel, n_state)];
  }
  CUDA_CALLABLE_MEMBER
  Real operator()(const int n_voxel, const int n_state=0) const {
    return vec[get_element_num(n_voxel, n_state)];
  }

//   CUDA_CALLABLE_MEMBER
//   Real & operator[](const int n) {
//     return vec[n];
//   }
//   CUDA_CALLABLE_MEMBER
//   Real operator[](const int n) const {
//     return vec[get_element_num(n)];
//   }



  void free_d_vec() {
#if defined(__CUDACC__) and not defined(__CUDA_ARCH__)
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
				 n_elements*sizeof(Real))
		      );
    if (transfer)
      //copy from host to host's device pointer
      checkCudaErrors(
		      cudaMemcpy(d_vec,
				 vec,
				 n_elements*sizeof(Real),
				 cudaMemcpyHostToDevice)
		      );
  }
  void to_host() {
    checkCudaErrors(
		    cudaMemcpy(vec,
			       d_vec,
			       n_elements*sizeof(Real),
			       cudaMemcpyDeviceToHost)
		    );
  }

#endif
};

template<int N_VOXELS, int N_STATES_PER_VOXEL>
class voxel_matrix {
public:
  static const int n_voxels = N_VOXELS;
  static const int n_states = N_STATES_PER_VOXEL;
  static const int n_elements = N_VOXELS*N_STATES_PER_VOXEL;

  MatrixX *eigen_mat;
  Real* mat;
  Real* d_mat = NULL;//pointer to device memory for CUDA

  CUDA_CALLABLE_MEMBER
  voxel_matrix() {
#ifndef __CUDA_ARCH__
    eigen_mat = new MatrixX;
    resize();
#endif
  }

  CUDA_CALLABLE_MEMBER
  ~voxel_matrix() {
#ifndef __CUDA_ARCH__
    delete eigen_mat;
    free_d_mat();
#endif
  }
  
  
  voxel_matrix(const voxel_matrix &copy) {
    *eigen_mat = *copy.eigen_mat;
    mat = eigen_mat->data();
  }

  voxel_matrix& operator=(const MatrixX &rhs) {
    assert(n_elements == rhs.size());
    *eigen_mat = rhs;
    mat = eigen_mat->data();
    return *this;
  }
  CUDA_CALLABLE_MEMBER
  voxel_matrix& operator=(const voxel_matrix &rhs) {
#ifdef __CUDA_ARCH__
    for (int i=0;i<n_elements*n_elements;i++)
      mat[i] = rhs.mat[i];
#else
    *eigen_mat = *rhs.eigen_mat;
    mat = eigen_mat->data();
#endif
    return *this;
  }
  operator MatrixX() const {
    return eigen();
  }


  void resize() {
    eigen_mat->resize(n_elements,n_elements);
    mat = eigen_mat->data();
  }

  MatrixX & eigen() {
    return *eigen_mat;
  }
  const MatrixX eigen() const {
    return *eigen_mat;
  }

  // element fetch
  CUDA_CALLABLE_MEMBER
  int get_element_num(const int n_voxel, const int n_state) const {
    return voxel_array<N_VOXELS, N_STATES_PER_VOXEL>::get_element_num(n_voxel, n_state);
  }
  
  //overload () to access voxels and states
  CUDA_CALLABLE_MEMBER
  Real & operator()(const int n_voxel, const int n_state,
		    const int m_voxel, const int m_state) {
    return operator()(get_element_num(n_voxel, n_state),
		      get_element_num(m_voxel, m_state));
  }
  CUDA_CALLABLE_MEMBER
  const Real & operator()(const int n_voxel, const int n_state,
			  const int m_voxel, const int m_state) const {
    return operator()(get_element_num(n_voxel, n_state),
		      get_element_num(m_voxel, m_state));
  }

  //overload () to access coefficients of mat directly
  CUDA_CALLABLE_MEMBER
  Real & operator()(const int n, const int m) {
#ifdef EIGEN_ROWMAJOR
    const int i = n*n_elements + m;//row major
#else
    const int i = m*n_elements + n;//column major
#endif
    return mat[i];
  }
  CUDA_CALLABLE_MEMBER
  const Real & operator()(const int n, const int m) const {
#ifdef EIGEN_ROWMAJOR
    const int i = n*n_elements + m;//row major
#else
    const int i = m*n_elements + n;//column major
#endif
    return mat[i];
  }

  void free_d_mat() {
#if defined(__CUDACC__) and not defined(__CUDA_ARCH__)
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
				 n_elements*n_elements*sizeof(Real))
		      );
    if (transfer)
      //copy from host to host's device pointer
      checkCudaErrors(
		      cudaMemcpy(d_mat,
				 mat,
				 n_elements*n_elements*sizeof(Real),
				 cudaMemcpyHostToDevice)
		      );
  }
  void to_host() {
    checkCudaErrors(
		    cudaMemcpy(mat,
			       d_mat,
			       n_elements*n_elements*sizeof(Real),
			       cudaMemcpyDeviceToHost)
		    );
  }
#endif
};

#endif
