#ifndef RT_gpu
#define RT_gpu

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <stdio.h>
#include <iostream>
#include "RT_grid.h"
#include "helper_cuda.h"

template <int N_EMISSIONS, typename grid_type, typename influence_type>
__global__
void brightness_kernel(const atmo_vector *obs_vecs, const int n_obs_vecs,
		       const RT_grid<N_EMISSIONS,grid_type,influence_type> *RT, 
		       //const grid_type *grid, const emission *emissions, const influence_type *transmission,
		       const double *g, brightness_tracker<N_EMISSIONS> *los,
		       const int n_subsamples = 5)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  
  for (int i_obs = index; i_obs < n_obs_vecs; i_obs += stride) {
    // if (blockIdx.x==0 && threadIdx.x==0) {
    //   printf("Hello from block %d, thread %d: i_obs = %d, RT->transmission.Tint_lerp(10.0) = %e\n",
    // 	     blockIdx.x, threadIdx.x, i_obs, RT->transmission.Tint_lerp(10.0));
    // }
    RT->brightness(obs_vecs[i_obs],g,
    		   los[i_obs],
    		   n_subsamples);
  }
}

template<int N_EMISSIONS, typename grid_type, typename influence_type>
void RT_grid<N_EMISSIONS,grid_type,influence_type>::brightness_gpu(observation &obs, const int n_subsamples/*=10*/) {
  cudaSetDevice(0);
  cudaFree(0);

  my_clock clk;
  clk.start();

  //move grid to GPU
  typedef RT_grid<N_EMISSIONS,grid_type,influence_type> this_RT_type;
  this_RT_type *d_RT;
  checkCudaErrors(
		  cudaMalloc((void **)&d_RT, sizeof(this_RT_type))
		  );
  checkCudaErrors(
		  cudaMemcpy(d_RT, this, sizeof(this_RT_type), cudaMemcpyHostToDevice)
		  );

  // printf("Host: transmission.Tint_lerp(10.0) = %e\n",
  // 	 transmission.Tint_lerp(10.0));
  
  //move observation vectors and g values to gpu
  int n_obs_vecs = obs.size();
  atmo_vector *d_obs_vecs;
  checkCudaErrors(
		  cudaMalloc(&d_obs_vecs, n_obs_vecs*sizeof(atmo_vector))
		  );
  checkCudaErrors(
		  cudaMemcpy(d_obs_vecs, obs.get_vecs().data(), n_obs_vecs*sizeof(atmo_vector), cudaMemcpyHostToDevice)
		  );
  double *d_g;
  checkCudaErrors(
		  cudaMalloc(&d_g, n_emissions*sizeof(double))
		  );
  checkCudaErrors(
		  cudaMemcpy(d_g, obs.emission_g_factors.data(), n_emissions*sizeof(double), cudaMemcpyHostToDevice)
		  );

  //prepare brightness objects
  brightness_tracker<n_emissions> *d_los;
  checkCudaErrors(
		  cudaMalloc(&d_los, n_obs_vecs*sizeof(brightness_tracker<n_emissions>) )
		  );

  //run kernel on GPU
  int blockSize = 512;
  int numBlocks = (n_obs_vecs + blockSize - 1) / blockSize;

  my_clock kernel_clk;
  kernel_clk.start();

  brightness_kernel<N_EMISSIONS,
		    grid_type,
		    influence_type><<<numBlocks,blockSize>>>(d_obs_vecs, n_obs_vecs,
							     d_RT,
							     d_g, d_los,
							     n_subsamples);

  kernel_clk.stop();
  kernel_clk.print_elapsed("brightness kernel execution takes ");


  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );

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
  
  checkCudaErrors(cudaFree(d_RT));
  checkCudaErrors(cudaFree(d_obs_vecs));
  checkCudaErrors(cudaFree(d_g));
  checkCudaErrors(cudaFree(d_los));
}



#endif
