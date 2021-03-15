//grid.h -- base class for atmosphere grids

#ifndef __GRID_H
#define __GRID_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "emission.hpp"
#include "boundaries.hpp"
#include "atm/atmosphere_base.hpp"

template <int NDIM, int NVOXELS, int NRAYS, int N_MAX_INTERSECTIONS, typename derived>
struct grid {
  static const int n_dimensions = NDIM; //dimensionality of the grid
  int n_pts[NDIM];
  static const int n_voxels = NVOXELS;//number of grid voxels

  //helper functions to swap between voxel and coordinate indices
  CUDA_CALLABLE_MEMBER
  void indices_to_voxel(const int (&indices)[n_dimensions], int & ret) const {
    static_cast<const derived*>(this)->indices_to_voxel(indices, ret);
  }
  CUDA_CALLABLE_MEMBER
  void voxel_to_indices(const int i_voxel, int (&indices)[n_dimensions]) const {
    static_cast<const derived*>(this)->voxel_to_indices(i_voxel, indices);
  }
  CUDA_CALLABLE_MEMBER
  void point_to_indices(const atmo_point &pt, int (&indices)[n_dimensions]) const {
    static_cast<const derived*>(this)->point_to_indices(pt, indices);
  };
  
  Real rmax,rmin;//max and min altitudes in the atmosphere

  //where the sun is, for single scattering
  const Real sun_direction[3] = {0., 0., 1.};
  
  //points inside the voxels to shoot rays from
  atmo_voxel voxels[NVOXELS];
  void setup_voxels(const atmosphere& atm) {
    static_cast<derived*>(this)->setup_voxels(atm);
  };

  //ray info
  static const int n_rays = NRAYS;
  atmo_ray rays[NRAYS];
  void setup_rays() {
    static_cast<derived*>(this)->setup_rays();
  } 
  
  //how to intersect rays with voxel boundaries
  static const int n_max_intersections = N_MAX_INTERSECTIONS;
  CUDA_CALLABLE_MEMBER
  void ray_voxel_intersections(const atmo_vector &vec,
			       boundary_intersection_stepper<n_dimensions, n_max_intersections> &stepper) const {
    static_cast<const derived*>(this)->ray_voxel_intersections(vec, stepper);
  } 
  
  //function to get interpolation coefs
  static const int n_interp_points = 2*n_dimensions;
  CUDA_CALLABLE_MEMBER
  void interp_weights(const int &ivoxel, const atmo_point &pt,
		      int (&indices)[n_interp_points], Real (&weights)[n_interp_points],
		      int (&indices_1d)[2*n_dimensions], Real (&weights_1d)[n_dimensions] ) const {
    static_cast<const derived*>(this)->interp_weights(ivoxel, pt,
						      indices, weights,
						      indices_1d, weights_1d);
  }
  
  template <typename E>
  void save_S(const string &fname, const E* const *emiss, const int n_emissions) const {
    static_cast<const derived*>(this)->save_S(fname, emiss, n_emissions);
  }
};

#endif

