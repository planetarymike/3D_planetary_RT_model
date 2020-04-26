#ifndef __RT_spherical_azimuthally_symmetric_gpu
#define __RT_spherical_azimuthally_symmetric_gpu
 
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <stdio.h>

__global__
void traverse_kernel(int n_dimensions,
		     atmo_vector *obs_vecs, int n_obs_vecs,
		     sphere *radial_boundary_spheres, int n_spheres,
		     cone *radial_boundary_cones, int n_cones)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  //  printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);

  double distances[2];
  int n_hits;
  
  for (int i_obs = index; i_obs < n_obs_vecs; i_obs += stride)
    for (int i_sphere=0;i_sphere<n_spheres;i_sphere++) {
      radial_boundary_spheres[i_sphere].intersections(obs_vecs[i_obs],distances,n_hits);
      if (threadIdx.x==0 && i_sphere==0) {
	for (int i_hit=0;i_hit<n_hits;i_hit++) {
	  printf("Distance is: %d\n", distances[i_hit]);
	}
      }
    }
  
}

void spherical_azimuthally_symmetric_grid::traverse_gpu(observation &obs, const int n_subsamples=5) {
  //move grid to gpu
  sphere* d_radial_boundary_spheres; 
  cone* d_sza_boundary_cones; 

  cudaMalloc(&d_radial_boundary_spheres, radial_boundary_spheres.size()*sizeof(sphere)); 
  cudaMalloc(&d_sza_boundary_cones, sza_boundary_cones.size()*sizeof(cone));
  
  cudaMemcpy(d_radial_boundary_spheres, radial_boundary_spheres.data(), radial_boundary_spheres.size()*sizeof(sphere), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sza_boundary_cones, sza_boundary_cones.data(), sza_boundary_cones.size()*sizeof(cone), cudaMemcpyHostToDevice);
  
  //move observation vectors to gpu
  thrust::device_vector<atmo_vector> dobs_vecs = obs.get_vecs();
  atmo_vector* dobs_vecs_array = thrust::raw_pointer_cast(&dobs_vecs[0]);
  int n_obs_vecs = dobs_vecs.size();
  
  // Run kernel on 1M elements on the GPU
  int blockSize = 512;
  int numBlocks = (obs.size() + blockSize - 1) / blockSize;
  
  traverse_kernel<<<numBlocks, blockSize>>>(n_dimensions,
					    dobs_vecs_array, n_obs_vecs,
					    d_radial_boundary_spheres, radial_boundary_spheres.size(),
					    d_sza_boundary_cones, sza_boundary_cones.size());

  cudaFree(d_radial_boundary_spheres);
  cudaFree(d_sza_boundary_cones);
}



#endif
