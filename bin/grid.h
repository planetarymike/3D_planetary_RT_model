//grid.h -- base class for atmosphere grids

#ifndef __GRID_H
#define __GRID_H

#include "cuda_compatibility.h"
#include "emission.h"

template <int NDIM, int NVOXELS, int NRAYS>
struct grid {
  static const int n_dimensions = NDIM; //dimensionality of the grid
  static const int n_voxels = NVOXELS;//number of grid voxels
  //helper functions to swap between voxel and coordinate indices
  CUDA_CALLABLE_MEMBER virtual void indices_to_voxel(const int (&/*indices*/)[n_dimensions], int & ret) const { };
  CUDA_CALLABLE_MEMBER virtual void voxel_to_indices(const int /*i_voxel*/, int (&/*indices*/)[n_dimensions]) const { };
  CUDA_CALLABLE_MEMBER virtual void point_to_indices(const atmo_point /*pt*/, int (&/*indices*/)[n_dimensions]) const { };
  
  virtual void setup_voxels() {};
  double rmax,rmin;//max and min altitudes in the atmosphere
  
  //points inside the voxels to shoot rays from
  atmo_point pts[NVOXELS];
  
  //ray info
  static const int n_rays = NRAYS;
  atmo_ray rays[NRAYS];
  virtual void setup_rays() {}; 
  
  //how to intersect rays with voxel boundaries
  CUDA_CALLABLE_MEMBER
  virtual void ray_voxel_intersections(const atmo_vector &vec,
				       boundary_intersection_stepper<n_dimensions> &stepper) const { }; 
  
  //where the sun is, for single scattering
  vector<double> sun_direction;
  
  //function to get interpolation coefs
  static const int n_interp_points = 2*n_dimensions;
  CUDA_CALLABLE_MEMBER
  virtual void interp_weights(const int &ivoxel, const atmo_point &pt,
			      int (&/*indices*/)[n_interp_points], double (&/*weights*/)[n_interp_points] ) const { };


#ifdef __CUDACC__
  template <typename G>
  void copy_to_cuda(G *d_ptr) {
    //declare, allocate, and copy all of the subelements
    atmo_point * d_pts;
    cudaMalloc(&d_pts, this->n_voxels*sizeof(atmo_point));
    cudaMemcpy(d_pts, this->pts, this->n_voxels*sizeof(atmo_point), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_ptr->pts), &d_pts, sizeof(atmo_point*), cudaMemcpyHostToDevice);

    atmo_ray* d_rays;
    cudaMalloc(&d_rays, this->n_voxels*sizeof(atmo_ray));
    cudaMemcpy(d_rays, this->rays, this->n_voxels*sizeof(atmo_ray), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_ptr->rays), &d_rays, sizeof(atmo_ray*), cudaMemcpyHostToDevice);
  }
  void cuda_free() {
    //free the pointers? not sure if they're still in scope
    
  }
#endif


  
  virtual void save_S(const string &fname, const emission *emiss) const { };
};

#endif

