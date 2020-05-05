#ifndef RT_gpu
#define RT_gpu
 
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <stdio.h>

template <typename grid_type>
__global__ void traverse_kernel(atmo_vector *obs_vecs, int n_obs_vecs,
				grid_type *grid)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  
  boundary_intersection_stepper<grid_type::n_dimensions> stepper;


  //print some stuff out to check
  printf("Hello from block %d, thread %d: n_radial_boundaries = %d\n",
	 blockIdx.x, threadIdx.x, grid->n_radial_boundaries);

  printf("Hello from block %d, thread %d: test[2] = %d\n",
	 blockIdx.x, threadIdx.x, grid->test[2]);
  
  // for (int i_obs = index; i_obs < n_obs_vecs; i_obs += stride) {
  //   grid->ray_voxel_intersections(obs_vecs[i_obs], stepper);
  // }

}

template<typename grid_type, typename influence_type>
void RT_grid<grid_type,influence_type>::traverse_gpu(observation &obs, const int n_subsamples) {
  //move grid to gpu
  grid_type* d_grid;
  cudaMalloc((void **)&d_grid, sizeof(grid_type));
  cudaMemcpy(d_grid, &grid, sizeof(grid_type), cudaMemcpyHostToDevice);
  grid.copy_to_cuda(d_grid);
  
  //move observation vectors to gpu
  atmo_vector *d_obs_vecs;
  cudaMalloc(&d_obs_vecs, obs.size()*sizeof(atmo_vector));
  cudaMemcpy(d_obs_vecs, obs.get_vecs().data(), obs.size()*sizeof(atmo_vector), cudaMemcpyHostToDevice);
  int n_obs_vecs = obs.size();

  //run kernel on GPU
  int blockSize = 256;
  int numBlocks = (obs.size() + blockSize - 1) / blockSize;
  traverse_kernel<grid_type><<<numBlocks, blockSize>>>(d_obs_vecs, n_obs_vecs,
						       d_grid);

  grid.cuda_free();
  cudaFree(d_grid);
  cudaFree(d_obs_vecs);
}



#endif
