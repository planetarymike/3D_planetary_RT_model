//RT_gpu.xu --- defines methods of RT_grid that require the CUDA compiler
#include <stdio.h>
#include "RT_grid.hpp"
#include "helper_cuda.h"

template <int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::RT_to_device() {
  if (d_RT == NULL) {
    //move grid to GPU
    checkCudaErrors(
		    cudaMalloc((void **)&d_RT, sizeof(RT_grid_type))
		    );
    checkCudaErrors(
		    cudaMemcpy(d_RT, this, sizeof(RT_grid_type), cudaMemcpyHostToDevice)
		    );
  }
}
template <int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::emissions_to_device_brightness() {
  for (int i_emission=0;i_emission<n_emissions;i_emission++)
    emissions[i_emission].copy_to_device_brightness(&(d_RT->emissions[i_emission]));
}
template <int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::emissions_to_device_influence() {
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
__global__
void brightness_kernel(const atmo_vector *obs_vecs, const int n_obs_vecs,
		       const RT_grid<N_EMISSIONS,grid_type,influence_type> *RT, 
		       const Real *g, brightness_tracker<N_EMISSIONS> *los,
		       const int n_subsamples = 5)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  
  for (int i_obs = index; i_obs < n_obs_vecs; i_obs += stride) {
    // if (blockIdx.x==0 && threadIdx.x==0) {
    //   printf("Hello from block %d, thread %d: i_obs = %d, RT->emissions[0].species_density[0] = %d\n",
    // 	     blockIdx.x, threadIdx.x, i_obs, RT->emissions[0].species_density.vec[0]);
    // }
    RT->brightness(obs_vecs[i_obs],g,
    		   los[i_obs],
    		   n_subsamples);
  }
}

template<int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::brightness_gpu(observation<n_emissions> &obs, const int n_subsamples/*=10*/) {
  cudaSetDevice(0);
  cudaFree(0);

  my_clock clk;
  clk.start();

  // move grid to GPU, copying only necessary members
  if (d_RT == NULL)
    RT_to_device();
  emissions_to_device_brightness();
  
  //move observation vectors and g values to gpu
  int n_obs_vecs = obs.size();
  atmo_vector *d_obs_vecs;
  checkCudaErrors(
		  cudaMalloc(&d_obs_vecs, n_obs_vecs*sizeof(atmo_vector))
		  );
  checkCudaErrors(
		  cudaMemcpy(d_obs_vecs, obs.get_vecs().data(), n_obs_vecs*sizeof(atmo_vector), cudaMemcpyHostToDevice)
		  );
  Real *d_g;
  checkCudaErrors(
		  cudaMalloc(&d_g, n_emissions*sizeof(Real))
		  );
  checkCudaErrors(
		  cudaMemcpy(d_g, obs.emission_g_factors, n_emissions*sizeof(Real), cudaMemcpyHostToDevice)
		  );

  //prepare brightness objects
  brightness_tracker<n_emissions> *d_los;
  checkCudaErrors(
		  cudaMalloc(&d_los, n_obs_vecs*sizeof(brightness_tracker<n_emissions>) )
		  );

  //run kernel on GPU
  int blockSize = 32;
  int numBlocks = (n_obs_vecs + blockSize - 1) / blockSize;
  
  my_clock kernel_clk;
  kernel_clk.start();
  
  brightness_kernel<N_EMISSIONS,
  		    grid_type,
  		    influence_type><<<numBlocks,blockSize>>>(d_obs_vecs, n_obs_vecs,
  							     d_RT,
  							     d_g, d_los,
  							     n_subsamples);
  
  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );
  
  kernel_clk.stop();
  kernel_clk.print_elapsed("brightness kernel execution takes ");


  //retrieve brightness from GPU
  brightness_tracker<n_emissions> *los;
  los = new brightness_tracker<n_emissions>[n_obs_vecs];
  checkCudaErrors(
		  cudaMemcpy(los, d_los, n_obs_vecs*sizeof(brightness_tracker<n_emissions>), cudaMemcpyDeviceToHost)
		  );

  for (int i_obs=0;i_obs<obs.size();i_obs++) {
    for (int i_emission=0;i_emission<n_emissions;i_emission++) {
      obs.brightness[i_obs][i_emission]   = los[i_obs].brightness[i_emission];
      obs.tau_species[i_obs][i_emission]  = los[i_obs].tau_species_final[i_emission];
      obs.tau_absorber[i_obs][i_emission] = los[i_obs].tau_absorber_final[i_emission];
    }
  }

  delete [] los;

  clk.stop();
  clk.print_elapsed("brightness memory operations take ",
		    kernel_clk.elapsed());
  
  checkCudaErrors(cudaFree(d_obs_vecs));
  checkCudaErrors(cudaFree(d_g));
  checkCudaErrors(cudaFree(d_los));
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
  RT_to_device();
  emissions_to_device_influence();
  
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

  //move info back to host
  emissions_influence_to_host();

  //solve for the source function
  solve();//shoule likely move this to GPU also

  // print time elapsed
  clk.stop();
  clk.print_elapsed("source function generation takes ");
  std::cout << std::endl;
    
  return;
}

