// RT_plane_parallel.h -- plane parallel RT grid

#ifndef __RT_plane_parallel
#define __RT_plane_parallel

#include <RT_grid.h>
#include "coordinate_generation.h"
#include "boundaries.h"
#include "atmosphere.h"
#include <type_traits> // std::is_scalar
#include <fstream>
#include <cmath>
#include "intersections.h"

struct plane_parallel_grid : RT_grid
{

  vector<double> radial_boundaries;
  vector<double> pts_radii;
  vector<plane> radial_boundary_planes;
  
  plane_parallel_grid(const vector<string> &emission_names,
		      holstein_approx &transmissionn)
    : RT_grid(emission_names,transmissionn)
  {
    n_dimensions = 1;
    sun_direction = {0.,0.,1.};
  }
  
  void setup_voxels(int n_radial_boundaries, atmosphere atm) {
    rmin = atm.rmin;
    rmax = atm.rmax;

    get_radial_log_linear_points(radial_boundaries, 
				 n_radial_boundaries,
				 atm.rmin, atm.rexo, atm.rmax);

    int n_voxelss = n_radial_boundaries-1;
    init_arrays(n_voxelss);

    pts.resize(n_voxels);
    pts_radii.resize(n_voxels);

    for (int i=0; i<n_radial_boundaries-1; i++) {
      pts_radii[i]=sqrt(radial_boundaries[i]*radial_boundaries[i+1]);
      pts[i].xyz(0.,0.,pts_radii[i]);
      pts[i].set_voxel_index(i);
    }


    radial_boundary_planes.resize(n_radial_boundaries);
    for (int i=0; i<n_radial_boundaries; i++) 
      radial_boundary_planes[i].z = radial_boundaries[i];

    voxels_init=true;
    if (rays_init)
      grid_init=true;
  }
  
  void setup_rays(int n_rayss) {
    //in a plane-parallel grid only the angle with the vertical matters
    n_rays = n_rayss;
    vector<double> ray_theta;
    vector<double> ray_weights;

    gauss_quadrature_points(ray_theta,ray_weights,0,pi,n_rays);
    
    rays.resize(n_rays);

    for (int i=0;i<n_rays;i++) {
      rays[i].tp(ray_theta[i],0.0);
      rays[i].set_ray_index(i,ray_weights[i],2*pi);
    }

    rays_init=true;
    if (voxels_init)
      grid_init = true;
  }

  int indices_to_voxel(const vector<int> &indices) {
    unsigned int comp_idx = indices[0];
    
    if ((comp_idx < 0) || (comp_idx > radial_boundaries.size()-2))
      return -1;
    else
      return comp_idx;
  }

  vector<int> voxel_to_indices(const int i_voxel) {
    vector<int> v(n_dimensions,-1);
    if ((i_voxel < 0) || (i_voxel > n_voxels-1))
      v[0]=-1;
    else
      v[0]=i_voxel;
    
    return v;
  }

  int find_coordinate_index(double pt_coord, vector<double> boundaries) {
    unsigned int i=-1;
    
    for (i=0;i<boundaries.size();i++)
      if (pt_coord < boundaries[i])
	break;

    if (i != boundaries.size())
      i--;

    return i;
  }

  vector<int> point_to_indices(const atmo_point pt) {
    vector<int> indices(n_dimensions,-1);
    indices[0] = find_coordinate_index(pt.r,radial_boundaries);
    return indices;
  }
  

  boundary_intersection_stepper ray_voxel_intersections(const atmo_vector &vec) {
    
    boundary_set boundaries(n_dimensions);

    //define the origin
    boundary origin(n_dimensions);
    origin.entering = vec.pt.i_voxel;
    if (vec.pt.i_voxel == -1)
      origin.entering_indices = point_to_indices(vec.pt);
    else
      origin.entering_indices = voxel_to_indices(origin.entering);
    origin.distance = 0.0;
    boundaries.append(origin);

    //do the intersections for each coordinate
    vector<double> temp_distances;
    for (unsigned int ir=0;ir<radial_boundaries.size();ir++) {
      temp_distances = radial_boundary_planes[ir].intersections(vec);
      if (temp_distances.size()>0) 
	boundaries.add_intersections(vec.pt.r, 0,
				     ir, radial_boundaries[ir], temp_distances);
    }

    //sort the list of intersections by distance & trim
    boundaries.sort();
    boundaries.propagate_indices();
    boundaries.assign_voxel_indices(this, &plane_parallel_grid::indices_to_voxel);
    boundaries.trim();

    return boundary_intersection_stepper(vec,
					 boundaries);
  }
  
  void save_S(string fname) {
    std::ofstream file(fname.c_str());
    if (file.is_open())
      {

	VectorXd r_boundaries_write_out = Eigen::Map<VectorXd>(radial_boundaries.data(),
							       radial_boundaries.size());

	file << "radial boundaries [cm]: " << r_boundaries_write_out.transpose() << "\n\n";

	VectorXd r_pts_write_out = Eigen::Map<VectorXd>(pts_radii.data(),
							pts_radii.size());

	file << "pts radii [cm]: " << r_pts_write_out.transpose() << "\n\n";
       	
	for (int i=0; i<n_emissions; i++)
	  file << "For " << emission_names[i] << ",\n"
	       << "Species density [cm-3]: " <<	species_density[i].transpose() << "\n"
	       << "Species single scattering tau: " <<	tau_species_single_scattering[i].transpose() << "\n"
	       << "Species cross section [cm2]: " << species_sigma[i].transpose() << "\n"
	       << "Absorber density [cm-3]: " << absorber_density[i].transpose() << "\n"
	       << "Absorber single scattering tau: " <<	tau_absorber_single_scattering[i].transpose() << "\n"
	       << "Source function: " << sourcefn[i].transpose() << "\n\n";
      }
  }

};




#endif
