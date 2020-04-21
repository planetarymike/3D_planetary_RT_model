#ifndef __RT_spherical_azimuthally_symmetric_gpu
#define __RT_spherical_azimuthally_symmetric_gpu
 
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

__global__
void traverse_kernel(int n_dimensions,
		     thrust::device_vector<atmo_vector> &obs_vecs,
		     sphere *radial_boundary_spheres, int n_spheres,
		     cone *radial_boundary_cones, int n_cones)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  thrust::device_vector<double> distances;
  
  for (int i = index; i < obs_vecs.size(); i += stride)
    for (int i_sphere=0;i_sphere<n_spheres;i_sphere++)
      radial_boundary_spheres[i_sphere].intersections(obs_vecs[i],distances);


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

  // Run kernel on 1M elements on the GPU
  int blockSize = 128;
  int numBlocks = (obs.size() + blockSize - 1) / blockSize;
  
  traverse_kernel<<<numBlocks, blockSize>>>(n_dimensions,
					    dobs_vecs,
					    d_radial_boundary_spheres, radial_boundary_spheres.size(),
					    d_sza_boundary_cones, sza_boundary_cones.size());

  cudaFree(d_radial_boundary_spheres);
  cudaFree(d_sza_boundary_cones);
}



#endif
