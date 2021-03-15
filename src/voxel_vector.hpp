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
class voxel_vector {
public:
  static const int n_voxels = N_VOXELS;

  VectorX *eigen_vec;
  Real* vec;
  Real* d_vec = NULL;//pointer to device memory for CUDA
  
#if defined(__CUDACC__) and defined(USE_CUDA_TEXTURES)
  bool texture;
  int n_dim; //number of grid dimensions
  int dim[3]; //size of each dimension
  cudaTextureObject_t tex;
  struct cudaResourceDesc resDesc;
  struct cudaTextureDesc texDesc;
#endif

  CUDA_CALLABLE_MEMBER
  voxel_vector() {
#ifndef __CUDA_ARCH__
    eigen_vec = new VectorX;
    resize();
#endif
#if defined(__CUDACC__) and defined(USE_CUDA_TEXTURES)
    texture = false;
#endif
  }

  CUDA_CALLABLE_MEMBER
  ~voxel_vector() {
#ifndef __CUDA_ARCH__
    delete eigen_vec;
    free_d_vec();
#endif
  }

  
  voxel_vector(const voxel_vector &copy) 
  {
    assert(n_voxels == copy.eigen_vec->size());
    *eigen_vec = *copy.eigen_vec;
    vec = eigen_vec->data();
  }
  voxel_vector& operator=(const VectorX &rhs) {
    assert(n_voxels == rhs.size());
    *eigen_vec = rhs;
    vec = eigen_vec->data();
    return *this;
  }
  CUDA_CALLABLE_MEMBER
  voxel_vector& operator=(const voxel_vector<N_VOXELS> &rhs) {
#ifdef __CUDA_ARCH__
    for (int i=0;i<n_voxels;i++)
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
    eigen_vec->resize(n_voxels);
    vec = eigen_vec->data();
  }

  VectorX & eigen() {
    return *eigen_vec;
  }
  const VectorX eigen() const {
    return *eigen_vec;
  }
  
  CUDA_CALLABLE_MEMBER
  Real & operator[](const int n) {
    return vec[n];
  }
  CUDA_CALLABLE_MEMBER
  Real operator[](const int n) const {
#if defined(__CUDA_ARCH__) and defined(USE_CUDA_TEXTURES)
    return fetch_element_device(n);
#else
    return vec[n];
#endif
  }
  CUDA_CALLABLE_MEMBER
  Real & operator()(const int n) {
    return vec[n];
  }
  CUDA_CALLABLE_MEMBER
  Real operator()(const int n) const {
#if defined(__CUDA_ARCH__) and defined(USE_CUDA_TEXTURES)
    return fetch_element_device(n);
#else
    return vec[n];
#endif
  }

#if defined(__CUDACC__) and defined(USE_CUDA_TEXTURES)
  __device__
  Real fetch_element_device(const int n) const {
    if (texture && n_dim == 2) {
      int i = n/dim[1];
      int j = n%dim[1];

      // Real **matvals;
      // matvals = new Real*[dim[0]];

      // for (int ii = 0; ii < dim[0]; ii++) {
      //   matvals[ii] = new Real[dim[1]];
      //   for (int jj = 0; jj < dim[1]; jj++)
      //     matvals[ii][jj] = tex2D<Real>(tex, jj+0.5, ii+0.5);
      // }

      // if (i==dim[0]-1)
      // 	i++;i--;
      
      Real retval = tex2D<Real>(tex,j+0.5f,i+0.5f);

      // for (int ii = 0; ii < dim[0]; ii++)
      // 	delete [] matvals[ii];
      // delete [] matvals;
      
      return retval;
    } else
      return vec[n];
  }

  __device__
  Real interp(const Real (&pt)[3], const int n_dim) const {
    if (n_dim==2)
      return tex2D<Real>(tex, pt[1], pt[0]);
    else
      return -1.0f;
  }
#endif

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
  void to_device_read_only(voxel_vector<N_VOXELS> *device_voxel_vector, const int n_dim, const int *dim) {
#if !defined(USE_CUDA_TEXTURES)
    to_device(/*transfer = */true);
#else
    //allocate and assign texture objects based on dimensionality
    if (n_dim == 2) {
      texture = true;

      const int num_rows = dim[0]; 
      const int num_cols = dim[1]; 

      // Real dataIn[num_rows*num_cols];
      // for (int i = 0; i < num_rows; i++)
      // 	for (int j=0;j<num_cols; j++) {
      // 	  int n = i*num_cols+j;
      // 	  dataIn[n] = i+j/100.0;
      // 	}

      
      //copy data to device
      size_t pitch;
      cudaMallocPitch((void**)&d_vec, &pitch,  num_cols*sizeof(Real), num_rows);
      cudaMemcpy2D(d_vec, pitch, vec, num_cols*sizeof(Real), num_cols*sizeof(Real), num_rows, cudaMemcpyHostToDevice);

      //move n_dim and dim information to device
      checkCudaErrors(
		      cudaMemcpy(&(device_voxel_vector->texture),
				 &texture,
				 sizeof(bool),
				 cudaMemcpyHostToDevice)
		      );
      checkCudaErrors(
		      cudaMemcpy(&(device_voxel_vector->n_dim),
				 &n_dim,
				 sizeof(int),
				 cudaMemcpyHostToDevice)
		      );
      checkCudaErrors(
		      cudaMemcpy(&(device_voxel_vector->dim[0]),
				 &dim[0],
				 sizeof(int),
				 cudaMemcpyHostToDevice)
		      );
      checkCudaErrors(
		      cudaMemcpy(&(device_voxel_vector->dim[1]),
				 &dim[1],
				 sizeof(int),
				 cudaMemcpyHostToDevice)
		      );

      //set up resource information
      memset(&resDesc, 0, sizeof(resDesc));
      resDesc.resType = cudaResourceTypePitch2D;
      resDesc.res.pitch2D.devPtr = d_vec;
      resDesc.res.pitch2D.width = num_cols;
      resDesc.res.pitch2D.height = num_rows;
      resDesc.res.pitch2D.desc = cudaCreateChannelDesc<Real>();
      resDesc.res.pitch2D.pitchInBytes = pitch;

      //bind texture
      memset(&texDesc, 0, sizeof(texDesc));
      texDesc.normalizedCoords = false;
      texDesc.filterMode       = cudaFilterModeLinear;
      texDesc.addressMode[0] = cudaAddressModeClamp;
      texDesc.addressMode[1] = cudaAddressModeClamp;
      texDesc.readMode = cudaReadModeElementType;

      cudaCreateTextureObject(&tex, &resDesc, &texDesc, NULL);

      //move texture information to device voxel_vec pointer
      checkCudaErrors(
      		      cudaMemcpy(&(device_voxel_vector->tex),
      				 &tex,
      				 sizeof(cudaTextureObject_t),
      				 cudaMemcpyHostToDevice)
      		      );
      checkCudaErrors(
		      cudaMemcpy(&(device_voxel_vector->resDesc),
				 &resDesc,
				 sizeof(cudaResourceDesc),
				 cudaMemcpyHostToDevice)
		      );
      checkCudaErrors(
		      cudaMemcpy(&(device_voxel_vector->texDesc),
				 &texDesc,
				 sizeof(cudaTextureDesc),
				 cudaMemcpyHostToDevice)
		      );
      //cudaCreateTextureObject(&(device_voxel_vector->tex), &(device_voxel_vector->resDesc), &(device_voxel_vector->texDesc, NULL);

    } else {
      to_device();
    }
#endif
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
class voxel_matrix {
public:
  static const int n_voxels = N_VOXELS;

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
  
  
  voxel_matrix(const voxel_matrix<N_VOXELS> &copy) {
    *eigen_mat = *copy.eigen_mat;
    mat = eigen_mat->data();
  }

  voxel_matrix& operator=(const MatrixX &rhs) {
    assert(n_voxels == rhs.size());
    *eigen_mat = rhs;
    mat = eigen_mat->data();
    return *this;
  }
  CUDA_CALLABLE_MEMBER
  voxel_matrix& operator=(const voxel_matrix<N_VOXELS> &rhs) {
#ifdef __CUDA_ARCH__
    for (int i=0;i<n_voxels*n_voxels;i++)
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
    eigen_mat->resize(n_voxels,n_voxels);
    mat = eigen_mat->data();
  }

  MatrixX & eigen() {
    return *eigen_mat;
  }
  const MatrixX eigen() const {
    return *eigen_mat;
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
