//RT_gpu.cu --- defines methods of RT_grid that require the CUDA compiler
#include "RT_grid.hpp"

#include <stdio.h>
#include <helper_cuda.h>
#include <cusolverDn.h>


template <int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::RT_to_device() {
  if (d_RT == NULL) {
    //move grid to GPU
    checkCudaErrors(
		    cudaMalloc((void **)&d_RT,
			       sizeof(RT_grid_type))
		    );
  }
  checkCudaErrors(
		  cudaMemcpy(d_RT,
			     this,
			     sizeof(RT_grid_type),
			     cudaMemcpyHostToDevice)
		  );
}

template <int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::RT_to_device_brightness() {
  RT_to_device();
  //emissions must be copied seperately as these contain dynamically allocated arrays
  for (int i_emission=0;i_emission<n_emissions;i_emission++)
    emissions[i_emission].copy_to_device_brightness(&(d_RT->emissions[i_emission]));
}
template <int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::RT_to_device_influence() {
  RT_to_device();
  //emissions must be copied seperately as these contain dynamically allocated arrays
  for (int i_emission=0;i_emission<n_emissions;i_emission++)
    emissions[i_emission].copy_to_device_influence(&(d_RT->emissions[i_emission]));
}

template <int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::emissions_influence_to_host() {
  //everything we need to copy back is in emissions
  for (int i_emission=0;i_emission<n_emissions;i_emission++)
    emissions[i_emission].copy_influence_to_host();
}

template <int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::emissions_solved_to_host() {
  //everything we need to copy back is in emissions
  for (int i_emission=0;i_emission<n_emissions;i_emission++)
    emissions[i_emission].copy_solved_to_host();
}

template <int N_EMISSIONS, typename grid_type, typename influence_type>
__global__
void brightness_kernel(const RT_grid<N_EMISSIONS,grid_type,influence_type> *RT, 
		       observation<N_EMISSIONS> *obs,
		       const int n_subsamples = 5)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  
  for (int i_obs = index; i_obs < obs->size(); i_obs += stride) {
    // if (blockIdx.x==0 && threadIdx.x==0) {
    //   printf("Hello from block %d, thread %d: i_obs = %d, RT->emissions[0].species_density[0] = %d\n",
    // 	     blockIdx.x, threadIdx.x, i_obs, RT->emissions[0].species_density.vec[0]);
    // }
    RT->brightness(obs->get_vec(i_obs),obs->emission_g_factors,
    		   obs->los[i_obs],
    		   n_subsamples);
  }
}

template<int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::brightness_gpu(observation<n_emissions> &obs, const int n_subsamples/*=10*/) {
  cudaSetDevice(0);
  cudaFree(0);

  my_clock clk;
  clk.start();

  // move grid to GPU
  RT_to_device_brightness();
  
  //move observation vectors and g values to gpu
  obs.to_device();

  //run kernel on GPU
  int blockSize = 32;
  int numBlocks = (obs.size() + blockSize - 1) / blockSize;
  
  my_clock kernel_clk;
  kernel_clk.start();
  
  brightness_kernel<N_EMISSIONS,
  		    grid_type,
  		    influence_type><<<numBlocks,blockSize>>>(d_RT,
							     obs.d_obs,
  							     n_subsamples);
  
  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );
  
  kernel_clk.stop();
  kernel_clk.print_elapsed("brightness kernel execution takes ");


  //retrieve brightness from GPU
  obs.to_host();

  clk.stop();
  clk.print_elapsed("brightness memory operations take ",
		    kernel_clk.elapsed());
  
}





template <int N_EMISSIONS, typename grid_type, typename influence_type>
__global__
void influence_kernel(RT_grid<N_EMISSIONS,grid_type,influence_type> *RT)
{
  //each thread runs one line of sight for one voxel
  int i_pt = blockIdx.x;
  int i_ray = threadIdx.x;

  //initialize matrix to zero
  for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++)
    for (int j_pt_index = threadIdx.x; j_pt_index < grid_type::n_voxels; j_pt_index+=blockDim.x) {
      int j_pt = j_pt_index;
      RT->emissions[i_emission].influence_matrix(i_pt,j_pt) = 0;
    }

  //initialize objects
  atmo_vector vec;
  influence_tracker<N_EMISSIONS,grid_type::n_voxels> temp_influence;

  //get the vector for this thread and run through the grid
  vec = atmo_vector(RT->grid.pts[i_pt], RT->grid.rays[i_ray]);
  temp_influence.reset();
  RT->voxel_traverse(vec,
		     &RT_grid<N_EMISSIONS,grid_type,influence_type>::influence_update,
		     temp_influence);

  __syncthreads();
  
  // reduce temp_influence across the thread block
  for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++)
    for (int j_pt_index = threadIdx.x; j_pt_index < grid_type::n_voxels+threadIdx.x; j_pt_index++) {
      int j_pt = j_pt_index % grid_type::n_voxels;
      RT->emissions[i_emission].influence_matrix(i_pt,j_pt) += temp_influence.influence[i_emission][j_pt];
      __syncthreads();
    }

  //now compute the single scattering function:
  //only one thread needs to do this
  if (threadIdx.x == 0) {
    temp_influence.reset();
    RT->get_single_scattering(RT->grid.pts[i_pt], temp_influence);
  }
}


template<int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::generate_S_gpu() {
  cudaSetDevice(0);
  cudaFree(0);
  
  //start timing
  my_clock clk;
  clk.start();

  //move grid to GPU
  RT_to_device_influence();
  
  //run kernel on GPU
  int blockSize = grid.n_rays;
  int numBlocks = grid.n_voxels;
  
  my_clock kernel_clk;
  kernel_clk.start();
  
  influence_kernel<N_EMISSIONS,
		   grid_type,
		   influence_type><<<numBlocks,blockSize>>>(d_RT);
  
  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );
  
  kernel_clk.stop();
  kernel_clk.print_elapsed("influence matrix generation takes ");

  // //solve on CPU with Eigen
  // emissions_influence_to_host();
  // solve();

  //solve on GPU (~2.5x slower for first call)
  //much faster than CPU on subsequent calls
  solve_gpu();
  emissions_solved_to_host();
 
  // print time elapsed
  clk.stop();
  clk.print_elapsed("source function generation takes ");
  std::cout << std::endl;
  
  return;
}






template<int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::solve_emission_gpu(emission<grid_type::n_voxels> & emiss) {
  //the CUDA matrix library has a LOT more boilerplate than Eigen
  //this is adapted from the LU dense example here:
  //https://docs.nvidia.com/cuda/pdf/CUSOLVER_Library.pdf

  //before calling this function, ensure that
  // 1) emiss.influence_matrix represents the kernel, not the influence matrix
  // 2) emiss.influence_matrix is in COLUMN-MAJOR order
  // 3) emiss.sourcefn = emiss.singlescat (solution is found in-place)
  
  cusolverDnHandle_t cusolverH = NULL;

  cudaStream_t stream = NULL;
  cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
  cudaError_t cudaStat1 = cudaSuccess;
  cudaError_t cudaStat2 = cudaSuccess;
  const int m = grid_type::n_voxels;
  const int lda = m;
  const int ldb = m;

  int info = 0; /* host copy of error info */

  int *d_Ipiv = NULL; /* pivoting sequence */
  int *d_info = NULL; /* error info */
  int lwork = 0; /* size of workspace */
  Real *d_work = NULL; /* device workspace for getrf */

  /* step 1: create cusolver handle, bind a stream */
  status = cusolverDnCreate(&cusolverH);
  assert(CUSOLVER_STATUS_SUCCESS == status);
  cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
  assert(cudaSuccess == cudaStat1);
  status = cusolverDnSetStream(cusolverH, stream);
  assert(CUSOLVER_STATUS_SUCCESS == status);


  /* step 2: allocate device solver parameters */
  cudaStat1 = cudaMalloc ((void**)&d_Ipiv, sizeof(int) * m);
  cudaStat2 = cudaMalloc ((void**)&d_info, sizeof(int));
  assert(cudaSuccess == cudaStat1);
  assert(cudaSuccess == cudaStat2);

  /* step 3: determine and allocate working space of getrf */
  //DnS refers to dense single-precision, need to change this if working with doubles
  status = cusolverDnSgetrf_bufferSize(cusolverH,
				       m, m, emiss.influence_matrix.d_mat,
				       lda,
				       &lwork);
  assert(CUSOLVER_STATUS_SUCCESS == status);
  cudaStat1 = cudaMalloc((void**)&d_work, sizeof(Real)*lwork);
  assert(cudaSuccess == cudaStat1);


  /* step 4: now we can do LU factorization */
  //this does the work in place--- matrix is replaced with solution!
  //DnS refers to dense single-precision, need to change this if working with doubles
  status = cusolverDnSgetrf(cusolverH,
			    m, m, emiss.influence_matrix.d_mat,
			    lda, d_work,
			    d_Ipiv, d_info);

  cudaStat1 = cudaDeviceSynchronize();
  assert(CUSOLVER_STATUS_SUCCESS == status);
  assert(cudaSuccess == cudaStat1);

  cudaStat2 = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
  assert(cudaSuccess == cudaStat2);

  if ( 0 > info ){
    printf("%d-th parameter is wrong \n", -info);
    exit(1);
  }

  /* step 5: finally, we get the solution */
  //this is also done in-place, replacing the vector to be solved with the solution
  //DnS refers to dense single-precision, need to change this if working with doubles
  status = cusolverDnSgetrs(cusolverH,
			    CUBLAS_OP_N,
			    m,
			    1, /* nrhs */
			    emiss.influence_matrix.d_mat,
			    lda,
			    d_Ipiv,
			    emiss.sourcefn.d_vec,
			    ldb,
			    d_info);

  cudaStat1 = cudaDeviceSynchronize();
  assert(CUSOLVER_STATUS_SUCCESS == status);
  assert(cudaSuccess == cudaStat1);

  /* free resources */
  if (d_Ipiv ) cudaFree(d_Ipiv);
  if (d_info ) cudaFree(d_info);
  if (d_work ) cudaFree(d_work);
  if (cusolverH ) cusolverDnDestroy(cusolverH);
  if (stream ) cudaStreamDestroy(stream);
}

template <int N_EMISSIONS, typename grid_type, typename influence_type>
__global__
void prepare_for_solution(RT_grid<N_EMISSIONS,grid_type,influence_type> *RT)
{
  //each block prepares one row of the influence matrix
  int i_pt = blockIdx.x;

  //convert to kernel from influence (kernel = identity - influence)
  for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++)
    for (int j_pt_index = threadIdx.x; j_pt_index < grid_type::n_voxels; j_pt_index+=blockDim.x) {
      int j_pt = j_pt_index;
      RT->emissions[i_emission].influence_matrix(i_pt,j_pt) *= -1;
    }
  //add the identity matrix
  for (int i_emission=threadIdx.x; i_emission < N_EMISSIONS; i_emission+=blockDim.x)
    RT->emissions[i_emission].influence_matrix(i_pt,i_pt) += 1;

  
  //now copy singlescat to sourcefn, preparing for in-place solution
  for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++)
    for (int j_pt = threadIdx.x; j_pt < grid_type::n_voxels; j_pt += blockDim.x) {
      RT->emissions[i_emission].sourcefn[j_pt] = RT->emissions[i_emission].singlescat[j_pt];
    }
}

template <int N_EMISSIONS, typename grid_type, typename influence_type>
__global__
void get_log_sourcefn(RT_grid<N_EMISSIONS,grid_type,influence_type> *RT)
{
  //now put singlescat in place of sourcefn
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++)
    for (int j_pt = index; j_pt < grid_type::n_voxels; j_pt += stride) {
      if (RT->emissions[i_emission].sourcefn[j_pt] == 0)
	RT->emissions[i_emission].log_sourcefn[j_pt] = -1e5;
      else
	RT->emissions[i_emission].log_sourcefn[j_pt] = std::log(RT->emissions[i_emission].sourcefn[j_pt]);
    }
}

template <int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::transpose_influence_gpu() {
  //transpose the influence matrix so it's column major
  const Real unity = 1.0;
  const Real null  = 0.0;
  cublasHandle_t handle;
  //gettimeofday(&t1, NULL);
  cublasCreate(&handle);

  //allocate space for transpose;
  const int N = grid_type::n_voxels;
  Real *d_transpose;
  checkCudaErrors(
		  cudaMalloc((void **) &d_transpose,
			     N*N*sizeof(Real))
		  );
  
  for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++) {
    //S here means single precision
    //transpose into swap memory
    cublasSgeam(handle,
		CUBLAS_OP_T,
		CUBLAS_OP_N,
		N, N,
		&unity, emissions[i_emission].influence_matrix.d_mat, N,
		&null, emissions[i_emission].influence_matrix.d_mat, N,
		d_transpose, N);

    checkCudaErrors( cudaPeekAtLastError() );
    checkCudaErrors( cudaDeviceSynchronize() );

    //reassign to original location
    cublasSgeam(handle, 
		CUBLAS_OP_N, CUBLAS_OP_N,
		N, N,
		&unity, d_transpose, N,
		&null, d_transpose, N,
		emissions[i_emission].influence_matrix.d_mat, N);

    checkCudaErrors( cudaPeekAtLastError() );
    checkCudaErrors( cudaDeviceSynchronize() );
  }

  cublasDestroy(handle);
}


template<int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::solve_gpu() {

  //run kernel on GPU
  int blockSize = grid.n_rays;
  int numBlocks = grid.n_voxels;
  
  prepare_for_solution<N_EMISSIONS,
		       grid_type,
		       influence_type><<<numBlocks,blockSize>>>(d_RT);

  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );

#ifdef EIGEN_ROWMAJOR
  // we need to transpose the influence matrix from row-major to
  // column-major before the CUDA utils can solve it
  transpose_influence_gpu();
#endif
  
  //solve for each emission
  for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++) {
    //pass to the CUDA solver
    solve_emission_gpu(emissions[i_emission]);
  }
  
  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );

  //define log_sourcefn
  get_log_sourcefn<N_EMISSIONS,
		   grid_type,
		   influence_type><<<numBlocks,blockSize>>>(d_RT);

  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );

}
