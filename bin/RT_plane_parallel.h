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

struct plane_parallel_grid : RT_grid<1> //this is a 1d grid
{

  vector<double> radial_boundaries;
  vector<double> pts_radii;
  vector<plane> radial_boundary_planes;
  
  plane_parallel_grid(const vector<string> &emission_names,
		      holstein_approx &transmissionn)
    : RT_grid(emission_names,transmissionn)
  {
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
  }

  int indices_to_voxel(const int &comp_idx) {
    if ((comp_idx < 0) || (comp_idx > (int) radial_boundaries.size()-2))
      return -1;
    else
      return comp_idx;
  }
  int indices_to_voxel(const vector<int> &indices) {
    return indices_to_voxel(indices[0]);
  }
  int indices_to_voxel(const int (&indices)[n_dimensions]) {
    return indices_to_voxel(indices[0]);
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
  

  boundary_intersection_stepper<n_dimensions> ray_voxel_intersections(const atmo_vector &vec) {
    
    boundary_set<n_dimensions> boundaries;

    //define the origin
    boundary<n_dimensions> origin;
    origin.entering = vec.pt.i_voxel;
    if (vec.pt.i_voxel == -1)
      origin.set_entering_indices(point_to_indices(vec.pt));
    else
      origin.set_entering_indices(voxel_to_indices(origin.entering));
    origin.distance = 0.0;
    boundaries.append(origin);

    //do the intersections for each coordinate
    int n_hits = 0;
    static thread_local double temp_distances[2] = {-1,-1};
    for (unsigned int ir=0;ir<radial_boundaries.size();ir++) {
      radial_boundary_planes[ir].intersections(vec, temp_distances, n_hits);
      boundaries.add_intersections(vec.pt.r, 0,
				   ir, radial_boundaries[ir],
				   temp_distances, n_hits);
    }

    //sort the list of intersections by distance & trim
    boundaries.sort();
    boundaries.propagate_indices();
    boundaries.assign_voxel_indices(this, &plane_parallel_grid::indices_to_voxel);
    boundaries.trim();

    return boundary_intersection_stepper<n_dimensions>(vec,
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
       	
	for (auto&& emiss: emissions)
	  file << "For " << emiss.name << ",\n"
	       << "Species density [cm-3]: " <<	emiss.species_density.transpose() << "\n"
	       << "Species single scattering tau: " <<	emiss.tau_species_single_scattering.transpose() << "\n"
	       << "Species cross section [cm2]: " << emiss.species_sigma.transpose() << "\n"
	       << "Absorber density [cm-3]: " << emiss.absorber_density.transpose() << "\n"
	       << "Absorber single scattering tau: " <<	emiss.tau_absorber_single_scattering.transpose() << "\n"
	       << "Source function: " << emiss.sourcefn.transpose() << "\n\n";
      }
  }

};




#endif
