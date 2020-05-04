//grid.h -- base class for atmosphere grids

#ifndef __GRID_H
#define __GRID_H

#include "cuda_compatibility.h"
#include "emission.h"

template <int NDIM>
struct grid {
  static const int n_dimensions = NDIM; //dimensionality of the grid
  int n_voxels;//number of grid voxels
  //helper functions to swap between voxel and coordinate indices
  CUDA_CALLABLE_MEMBER virtual void indices_to_voxel(const int (&/*indices*/)[n_dimensions], int & ret) const { };
  CUDA_CALLABLE_MEMBER virtual void voxel_to_indices(const int /*i_voxel*/, int (&/*indices*/)[n_dimensions]) const { };
  CUDA_CALLABLE_MEMBER virtual void point_to_indices(const atmo_point /*pt*/, int (&/*indices*/)[n_dimensions]) const { };
  
  virtual void setup_voxels() {};
  double rmax,rmin;//max and min altitudes in the atmosphere
  
  //points inside the voxels to shoot rays from
  vector<atmo_point> pts;
  
  //ray info
  int n_rays;
  vector<atmo_ray> rays;
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
  
  virtual void save_S(const string &fname, const vector<emission> &emiss) const { };
};

#endif

