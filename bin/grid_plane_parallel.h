// grid_plane_parallel.h -- plane parallel RT grid

#ifndef __grid_plane_parallel
#define __grid_plane_parallel

#include <grid.h>
#include "coordinate_generation.h"
#include "boundaries.h"
#include "atmosphere.h"
#include <fstream>
#include <cmath>
#include "intersections.h"

template <int N_RADIAL_BOUNDARIES, int N_RAYS_THETA, int N_RAYS_PHI>
struct plane_parallel_grid : grid<1,//this is a 1d grid
				  N_RADIAL_BOUNDARIES-1,
				  N_RAYS_THETA*N_RAYS_PHI,
				  N_RADIAL_BOUNDARIES> 
{

  using parent_grid = grid<1,
			   N_RADIAL_BOUNDARIES-1,
			   N_RAYS_THETA*N_RAYS_PHI,
			   N_RADIAL_BOUNDARIES>;

  static const int n_radial_boundaries = N_RADIAL_BOUNDARIES;
  Real radial_boundaries[n_radial_boundaries];
  Real pts_radii[n_radial_boundaries-1];
  plane radial_boundary_planes[n_radial_boundaries];
  
  plane_parallel_grid() {
    this->sun_direction = {0.,0.,1.};
  }
  
  void setup_voxels(atmosphere &atm)
  {
    this->rmin = atm.rmin;
    this->rmax = atm.rmax;

    get_radial_log_linear_points(radial_boundaries, 
				 n_radial_boundaries,
				 atm.rmin, atm.rexo, atm.rmax);

    for (int i=0; i<n_radial_boundaries-1; i++) {
      pts_radii[i]=sqrt(radial_boundaries[i]*radial_boundaries[i+1]);
      this->pts[i].xyz(0.,0.,pts_radii[i]);
      this->pts[i].set_voxel_index(i);
    }

    for (int i=0; i<n_radial_boundaries; i++) 
      radial_boundary_planes[i].z = radial_boundaries[i];
  }
  
  void setup_rays() {
    //in a plane-parallel grid only the angle with the vertical matters
    vector<Real> ray_theta;
    vector<Real> ray_weights;

    gauss_quadrature_points(ray_theta,ray_weights,0,pi,parent_grid::n_rays);
    
    for (int i=0;i<parent_grid::n_rays;i++) {
      this->rays[i].tp(ray_theta[i],0.0);
      this->rays[i].set_ray_index(i,ray_weights[i],2*pi);
    }
  }

  CUDA_CALLABLE_MEMBER 
  void indices_to_voxel(const int &comp_idx, int & ret) const {
    if ((comp_idx < 0) || (comp_idx > (int) n_radial_boundaries-2))
      ret = -1;
    else
      ret = comp_idx;
  }
  CUDA_CALLABLE_MEMBER 
  void indices_to_voxel(const int (&indices)[parent_grid::n_dimensions], int & ret) const {
    indices_to_voxel(indices[0], ret);
  }
  CUDA_CALLABLE_MEMBER 
  void voxel_to_indices(const int i_voxel, int (&ret)[parent_grid::n_dimensions]) const {
    if ((i_voxel < 0) || (i_voxel > parent_grid::n_voxels-1))
      ret[0]=-1;
    else
      ret[0]=i_voxel;
  }

  template <class V>
  CUDA_CALLABLE_MEMBER 
  int find_coordinate_index(const Real &pt_coord, const V &boundaries, int n_boundaries) const {
    int i;
    
    for (i=0;i<n_boundaries;i++)
      if (pt_coord < boundaries[i])
	break;

    i--;

    assert((boundaries[n_boundaries-1]<=pt_coord ||
	    pt_coord < boundaries[0] ||
	    (boundaries[i]<=pt_coord &&pt_coord<boundaries[i+1]))
	   && "we have found the appropriate point index");
    
    return i;
  }

  CUDA_CALLABLE_MEMBER 
  void point_to_indices(const atmo_point pt, int (&indices)[parent_grid::n_dimensions]) const {
    indices[0] = find_coordinate_index(pt.r,radial_boundaries,n_radial_boundaries);
  }
  
  CUDA_CALLABLE_MEMBER 
  void ray_voxel_intersections(const atmo_vector &vec,
			       boundary_intersection_stepper<parent_grid::n_dimensions,
			                                     parent_grid::n_max_intersections> &stepper) const {
    stepper.vec = vec;
    stepper.boundaries.reset();

    //define the origin
    boundary<parent_grid::n_dimensions> origin;
    origin.entering = vec.pt.i_voxel;
    if (vec.pt.i_voxel == -1)
      point_to_indices(vec.pt, origin.entering_indices);
    else
      voxel_to_indices(origin.entering, origin.entering_indices);
    origin.distance = 0.0;
    stepper.boundaries.append(origin);
    
    //do the intersections for each coordinate
    int n_hits = 0;
    Real temp_distances[2] = {-1,-1};
    for (unsigned int ir=0;ir<n_radial_boundaries;ir++) {
      radial_boundary_planes[ir].intersections(vec, temp_distances, n_hits);
      stepper.boundaries.add_intersections(vec.pt.r, 0,
					   ir, radial_boundaries[ir],
					   temp_distances, n_hits);
    }

    //sort the list of intersections by distance & trim
    stepper.boundaries.sort();
    stepper.boundaries.propagate_indices();
    stepper.boundaries.assign_voxel_indices(this);
    stepper.boundaries.trim();

    stepper.init_stepper();
  }
  
  void save_S(const string fname, const emission<parent_grid::n_voxels> *emissions, const int n_emissions) const {
    std::ofstream file(fname.c_str());
    if (file.is_open())
      {

	VectorX r_boundaries_write_out = Eigen::Map<const VectorX>(radial_boundaries,
								     n_radial_boundaries);

	file << "radial boundaries [cm]: " << r_boundaries_write_out.transpose() << "\n\n";

	VectorX r_pts_write_out = Eigen::Map<const VectorX>(pts_radii,
							      n_radial_boundaries-1);

	file << "pts radii [cm]: " << r_pts_write_out.transpose() << "\n\n";
       	
	for (int i_emission=0;i_emission<n_emissions;i_emission++)
	  file << "For " << emissions[i_emission].name << ",\n"
	       << "Species density [cm-3]: " <<	emissions[i_emission].species_density.transpose() << "\n"
	       << "Species single scattering tau: " <<	emissions[i_emission].tau_species_single_scattering.transpose() << "\n"
	       << "Species cross section [cm2]: " << emissions[i_emission].species_sigma.transpose() << "\n"
	       << "Absorber density [cm-3]: " << emissions[i_emission].absorber_density.transpose() << "\n"
	       << "Absorber single scattering tau: " <<	emissions[i_emission].tau_absorber_single_scattering.transpose() << "\n"
	       << "Source function: " << emissions[i_emission].sourcefn.transpose() << "\n\n";
      }
  }

};

#endif
