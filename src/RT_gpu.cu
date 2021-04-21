//RT_gpu.cu --- defines methods of RT_grid that require the CUDA compiler
#include "RT_grid.hpp"
#include <stdio.h>
#include <helper_cuda.h>
#include <cusolverDn.h>


template<typename emission_type, int N_EMISSIONS, typename grid_type>
void RT_grid<emission_type, N_EMISSIONS, grid_type>::RT_to_device() {
  if (d_RT == NULL) {
    //move grid to GPU
    checkCudaErrors(
		    cudaMalloc((void **)&d_RT,
			       sizeof(RT_grid_type))
		    );
    checkCudaErrors( cudaPeekAtLastError() );
    checkCudaErrors( cudaDeviceSynchronize() );
  }
  checkCudaErrors(
		  cudaMemcpy(d_RT,
			     this,
			     sizeof(RT_grid_type),
			     cudaMemcpyHostToDevice)
		  );
  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );

  for (int i_emission=0;i_emission<n_emissions;i_emission++) {
    //allocate and move the trivial members to the grid
    emissions[i_emission]->allocate_device_emission();
    emissions[i_emission]->copy_trivial_members_to_device();

    // now ensure the emission pointers are correct in the RT object
    checkCudaErrors(
		    cudaMemcpy(&(d_RT->emissions[i_emission]),
			       &(emissions[i_emission]->device_emission),
			       sizeof(emission_type*),
			       cudaMemcpyHostToDevice)
		    );
    checkCudaErrors( cudaPeekAtLastError() );
    checkCudaErrors( cudaDeviceSynchronize() );
  }
}

template<typename emission_type, int N_EMISSIONS, typename grid_type>
void RT_grid<emission_type, N_EMISSIONS, grid_type>::RT_to_device_brightness() {
  RT_to_device();
  //emissions must be copied seperately as these contain dynamically allocated arrays
  for (int i_emission=0;i_emission<n_emissions;i_emission++)
    emissions[i_emission]->copy_to_device_brightness();
}
template<typename emission_type, int N_EMISSIONS, typename grid_type>
void RT_grid<emission_type, N_EMISSIONS, grid_type>::RT_to_device_influence() {
  RT_to_device();
  //emissions must be copied seperately as these contain dynamically allocated arrays
  for (int i_emission=0;i_emission<n_emissions;i_emission++)
    emissions[i_emission]->copy_to_device_influence();
}

template<typename emission_type, int N_EMISSIONS, typename grid_type>
void RT_grid<emission_type, N_EMISSIONS, grid_type>::emissions_influence_to_host() {
  //everything we need to copy back is in emissions
  for (int i_emission=0;i_emission<n_emissions;i_emission++) {
    emissions[i_emission]->copy_influence_to_host();
  }
}

template<typename emission_type, int N_EMISSIONS, typename grid_type>
void RT_grid<emission_type, N_EMISSIONS, grid_type>::emissions_solved_to_host() {
  //everything we need to copy back is in emissions
  for (int i_emission=0;i_emission<n_emissions;i_emission++)
    emissions[i_emission]->copy_solved_to_host();
}

template<typename emission_type, int N_EMISSIONS, typename grid_type>
void RT_grid<emission_type, N_EMISSIONS, grid_type>::device_clear() {
  for (int i_emission=0;i_emission<n_emissions;i_emission++)
    emissions[i_emission]->device_clear();
  if(d_RT!=NULL)
    checkCudaErrors(cudaFree(d_RT));
  d_RT=NULL;
  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );
}


static const int brightness_blockSize = 32;

template<typename emission_type, int N_EMISSIONS, typename grid_type>
__global__
void brightness_kernel(const RT_grid<emission_type, N_EMISSIONS, grid_type> *__restrict__ RT, 
		       observation<emission_type, N_EMISSIONS> *obs,
		       const int n_subsamples = 10)
{
  if (threadIdx.x==0 && blockIdx.x==0)
    printf("size of RT grid: %i\n",(int) sizeof(*RT));

  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  //shared objects for each thread to do the calculation with lower
  //memory latency --- ~5x speedup
  __shared__ atmo_vector obs_vecs[brightness_blockSize];

  __shared__ typename emission_type::brightness_tracker los[brightness_blockSize][N_EMISSIONS];
  for (int i_emission=0; i_emission<N_EMISSIONS; i_emission++)
    los[threadIdx.x][i_emission].init();

  typename emission_type::brightness_tracker *los_ptr[N_EMISSIONS];
  
  for (int i_obs = index; i_obs < obs->size(); i_obs += stride) {
    // if (blockIdx.x==0 && threadIdx.x==0) {
    //   printf("Hello from block %d, thread %d: i_obs = %d, RT->emissions[0].species_density[0] = %d\n",
    // 	     blockIdx.x, threadIdx.x, i_obs, RT->emissions[0].species_density.vec[0]);
    // }

    obs_vecs[threadIdx.x] = obs->get_vec(i_obs);
    //los[threadIdx.x] = obs->los[i_obs];

    for (int i_emission=0;i_emission<N_EMISSIONS;i_emission++)
      los_ptr[i_emission] = &los[threadIdx.x][i_emission];

    RT->brightness(obs_vecs[threadIdx.x],
		   los_ptr,
		   n_subsamples);

    for (int i_emission=0; i_emission<N_EMISSIONS; i_emission++) {
      // profiler shows a large number of samples here but copying
      // only parts of the object back doesn't improve speed
      obs->los[i_emission][i_obs] = los[threadIdx.x][i_emission];
    }
  }
}


template<typename emission_type, int N_EMISSIONS, typename grid_type>
void RT_grid<emission_type,
	     N_EMISSIONS,
	     grid_type>::brightness_gpu(observation<emission_type, N_EMISSIONS> &obs,
					const int n_subsamples/*=5*/) {
  cudaSetDevice(0);
  cudaFree(0);
  //  checkCudaErrors(cudaDeviceSetCacheConfig(cudaFuncCachePreferShared));//doesn't help speed things up
  
  my_clock clk;
  clk.start();

  // move grid to GPU
  RT_to_device_brightness();
  
  //move observation vectors and g values to gpu
  obs.to_device();

  //std::cout << "brightness_blockSize = " << brightness_blockSize << std::endl;

  //run kernel on GPU
  int numBlocks = (obs.size() + brightness_blockSize - 1) / brightness_blockSize;
  
  my_clock kernel_clk;
  kernel_clk.start();
  
  brightness_kernel<emission_type,
		    N_EMISSIONS,
  		    grid_type><<<numBlocks,
                                 brightness_blockSize>>>(d_RT,
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

  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );
  obs.device_clear();
  device_clear();
}





template<typename emission_type, int N_EMISSIONS, typename grid_type>
__global__
void influence_kernel(RT_grid<emission_type,N_EMISSIONS,grid_type> *RT)
{
  //each block runs one voxel;
  //each thread runs one line of sight for one voxel
  int i_vox = blockIdx.x;
  int i_ray = threadIdx.x;

  // reset each emission's solution
  for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++)
    RT->emissions[i_emission]->reset_solution(i_vox);
  __syncthreads();

  //initialize objects 
  atmo_vector vec;
  typename emission_type::influence_tracker temp_influence[N_EMISSIONS];

  // __shared__ atmo_vector vec[grid_type::n_rays];
  // __shared__ typename emission_type::influence_tracker temp_influence[grid_type::n_rays][N_EMISSIONS];
  // __shared__ Real voxel_influence[N_EMISSIONS][emission_type::n_elements];
  // // voxel_influence is shared across all trackers and updated using an atomic add

  vec.ptray(RT->grid.voxels[i_vox].pt, RT->grid.rays[i_ray]);
  for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++) {
    //temp_influence[i_ray][i_emission].influence.vec = voxel_influence[i_emission];
    temp_influence[i_emission].init();
    RT->emissions[i_emission]->reset_tracker(i_vox, temp_influence[i_emission]);
  }

  //__syncthreads();
  
  //integrate
  RT->voxel_traverse(vec,
  		     &RT_grid<emission_type,N_EMISSIONS,grid_type>::influence_update,
  		     temp_influence);

  //__syncthreads();
  
  // reduce temp_influence across the thread block
  for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++) 
    RT->emissions[i_emission]->accumulate_influence(i_vox, temp_influence[i_emission], threadIdx.x);

  //__syncthreads();
  
  //now compute the single scattering function:
  //only one thread needs to do this
  if (threadIdx.x == 0) {
    for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++)
      RT->emissions[i_emission]->reset_tracker(i_vox, temp_influence[i_emission]);
    RT->get_single_scattering(RT->grid.voxels[i_vox].pt, temp_influence);
  }

  //__syncthreads();
}


template<typename emission_type, int N_EMISSIONS, typename grid_type>
void RT_grid<emission_type, N_EMISSIONS, grid_type>::generate_S_gpu() {
  cudaSetDevice(0);
  cudaFree(0);

  //start timing
  my_clock clk;
  clk.start();

  //move grid to GPU
  RT_to_device_influence();
  
  //run kernel on GPU
  int numBlocks = grid.n_voxels;
  int blockSize = grid.n_rays;
  
  my_clock kernel_clk;
  kernel_clk.start();
  
  influence_kernel<emission_type,
		   N_EMISSIONS,
  		   grid_type><<<numBlocks,blockSize>>>(d_RT);
  
  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );
  
  kernel_clk.stop();
  kernel_clk.print_elapsed("influence matrix generation takes ");

  // //solve on CPU with Eigen
  emissions_influence_to_host();
  save_influence();
  // solve();

  //solve on GPU (~2.5x slower for first call)
  //much faster than CPU on subsequent calls
  solve_gpu();
  emissions_solved_to_host();

  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );
  device_clear();

  // print time elapsed
  clk.stop();
  clk.print_elapsed("source function generation takes ");
  std::cout << std::endl;
  
  return;
}

template<typename emission_type, int N_EMISSIONS, typename grid_type>
void RT_grid<emission_type, N_EMISSIONS, grid_type>::solve_gpu() {
  //solve for each emission
  for (int i_emission=0; i_emission < N_EMISSIONS; i_emission++) {
    //pass to the CUDA solver
    emissions[i_emission]->solve_gpu();
    checkCudaErrors( cudaPeekAtLastError() );
    checkCudaErrors( cudaDeviceSynchronize() );
  }

}
