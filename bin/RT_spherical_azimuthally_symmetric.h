// RT_spherical_azimuthally_symmetric.h -- spherical grid with symmetry about the Mars-Sun line

#ifndef __RT_spherical_azimuthally_symmetric
#define __RT_spherical_azimuthally_symmetric

#include <RT_grid.h>
#include "coordinate_generation.h"
#include "boundaries.h"
#include "atmosphere.h"
//#include "geometric_tools_interface.h"
#include "intersections.h"
#include <type_traits> 
#include <fstream>
#include <cmath>

struct spherical_azimuthally_symmetric_grid : RT_grid
{

  const int r_dimension;
  int n_radial_boundaries;
  vector<double> radial_boundaries;
  vector<double> pts_radii;
  vector<sphere> radial_boundary_spheres;
  // vector<Sphere3<double>> radial_boundary_spheres;
  // ray_sphere_query ray_sphere;

  const int sza_dimension;
  int n_sza_boundaries;
  vector<double> sza_boundaries;
  vector<double> pts_sza;
  vector<cone> sza_boundary_cones;
  // vector<Cone3<double>> sza_boundary_cones;
  // ray_cone_query ray_cone;


  vector<int> n_boundaries;
  intersection_writer saver;
  bool save_intersections;

  int n_theta;
  vector<double> ray_theta;
  int n_phi;
  vector<double> ray_phi;

  
 spherical_azimuthally_symmetric_grid(int n_emissions,
				      influence &transmissionn)
   : RT_grid(n_emissions,transmissionn),
    r_dimension(0), sza_dimension(1)
    {
      n_dimensions = 2;
      n_boundaries.resize(n_dimensions);
      sun_direction = {0.,0.,1.};
      
      
      save_intersections = false;
      saver.fname = "intersections.dat";
    }
  
  void setup_voxels(int n_radial_boundariess,
		    int n_sza_boundariess,
		    atmosphere atm)    
  {
    n_radial_boundaries = n_radial_boundariess;
      
    rmin = atm.rmin;
    rmax = atm.rmax;

    get_radial_log_linear_points(radial_boundaries, 
				 n_radial_boundaries,
				 atm.rmin, atm.rexo, atm.rmax);

    pts_radii.resize(n_radial_boundaries-1);
    for (int i=0; i<n_radial_boundaries-1; i++) {
      pts_radii[i]=sqrt(radial_boundaries[i]*radial_boundaries[i+1]);
    }

    radial_boundary_spheres.resize(n_radial_boundaries);
    // array<double,3> origin = {0.,0.,0.};
    for (int i=0; i<n_radial_boundaries; i++) {
      radial_boundary_spheres[i].set_radius(radial_boundaries[i]);
      // radial_boundary_spheres[i].center = origin; 
      // radial_boundary_spheres[i].radius = radial_boundaries[i];
    }



    n_sza_boundaries = n_sza_boundariess;
    
    double sza_spacing = pi / (n_sza_boundaries - 2.);
    sza_boundaries.resize(n_sza_boundaries);
    for (int i=0;i<n_sza_boundaries;i++) {
      sza_boundaries[i]=(i-0.5)*sza_spacing;
    }

    pts_sza.resize(n_sza_boundaries-1);
    for (unsigned int i=0;i<pts_sza.size();i++) {
      pts_sza[i]=0.5*(sza_boundaries[i] + sza_boundaries[i+1]);
    }
    
    sza_boundary_cones.resize(n_sza_boundaries-2);
    for (int i=0;i<n_sza_boundaries-2;i++) {
      sza_boundary_cones[i].set_angle(sza_boundaries[i+1]);
	
      // assert (sza_boundaries[i+1] != pi/2 && "sza_boundary = pi/2; choose an even number of sza boundaries.\n");
      // if (sza_boundaries[i+1] < pi/2) {
      // 	array<double,3> sun = {0.,0.,1.};
      // 	sza_boundary_cones[i] = Cone3<double>(Ray3<double>(origin, sun), sza_boundaries[i+1]);
      // } else (sza_boundaries[i+1] > pi/2) {
      // 	array<double,3> sun = {0.,0.,-1.};
      // 	sza_boundary_cones[i] = Cone3<double>(Ray3<double>(origin, sun), pi-sza_boundaries[i+1]);
      // }
    }


    int n_voxelss = pts_radii.size()*pts_sza.size();
    n_boundaries[r_dimension] = n_radial_boundaries-1;
    n_boundaries[sza_dimension] = n_sza_boundaries-1;
    init_arrays(n_voxelss);
    
    pts.resize(n_voxels);
    int ivoxel;
    for (unsigned int i=0; i<pts_radii.size(); i++) {
      for (unsigned int j=0;j<pts_sza.size();j++) {
	ivoxel = i*(n_sza_boundaries-1) + j;
	
	pts[ivoxel].rtp(pts_radii[i],pts_sza[j],0.);
	pts[ivoxel].set_voxel_index(ivoxel);
      }
    }

    if (save_intersections)
      saver.save_coordinates(radial_boundaries,sza_boundaries);
    
    voxels_init=true;
    if (rays_init)
      grid_init=true;
  }
  
  void setup_rays(int n_thetaa, int n_phii)
  {
    n_theta = n_thetaa;
    n_phi = n_phii;

    vector<double> ray_weights_theta;
    gauss_quadrature_points(ray_theta,ray_weights_theta,0,pi,n_theta);
    
    double phi_spacing = 2*pi/n_phi;
    ray_phi.resize(n_phi);
    for (int i=0;i<n_phi;i++)
      ray_phi[i] = (i+0.5)*phi_spacing;
    
    n_rays = n_theta*n_phi;
    rays.resize(n_rays);
    int iray;
    for (int i=0;i<n_theta;i++) {
      for (int j=0;j<n_phi;j++) {
	iray = i * n_phi + j;
	rays[iray].tp(ray_theta[i],ray_phi[j]);
	rays[iray].set_ray_index(iray, ray_weights_theta[i], phi_spacing);
      }
    }
    
    rays_init=true;
    if (voxels_init)
      grid_init = true;
  }
  
  int indices_to_voxel(const vector<int> &indices) {
    unsigned int r_idx = indices[r_dimension];
    unsigned int sza_idx = indices[sza_dimension];
    
    if ((r_idx   < 0) || (  r_idx > pts_radii.size()-1) ||
	(sza_idx < 0) || (sza_idx > pts_sza.size()-1))
      return -1;
    else
      return indices[r_dimension]*(n_sza_boundaries-1)+indices[sza_dimension];
  }

  vector<int> voxel_to_indices(const int i_voxel) {
    vector<int> v(n_dimensions,-1);
    if ((i_voxel < 0) || (i_voxel > n_voxels-1)) {
      v[r_dimension]=-1;
      v[sza_dimension]=-1;
    } else {
      v[r_dimension]=i_voxel / (n_sza_boundaries-1);
      v[sza_dimension]=i_voxel % (n_sza_boundaries-1);
    }
    
    return v;
  }

  int find_coordinate_index(double pt_coord, vector<double> boundaries) {
    unsigned int i;
    
    for (i=0;i<boundaries.size();i++)
      if (pt_coord < boundaries[i])
	break;

    return i--;
  }

  vector<int> point_to_indices(const atmo_point pt) {
    vector<int> indices(n_dimensions,-1);
    indices[r_dimension] = find_coordinate_index(pt.r, radial_boundaries);
    indices[sza_dimension] = find_coordinate_index(pt.t, sza_boundaries);
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
    for (unsigned int ir=0;ir<radial_boundary_spheres.size();ir++) {
      temp_distances = radial_boundary_spheres[ir].intersections(vec);
      if (temp_distances.size()>0) 
	boundaries.add_intersections(vec.pt.r, r_dimension,
				     ir, radial_boundaries[ir], temp_distances);
      // ray_sphere(vec.GT_ray3, radial_boundary_spheres[ir]);
      // if (ray_sphere.intersect) 
      // 	boundaries.add_intersections(vec.pt.r, 0,
      // 				     ir, radial_boundaries[ir], ray_sphere.distances);
    }

    for (unsigned int isza=0;isza<sza_boundary_cones.size();isza++) {
      temp_distances = sza_boundary_cones[isza].intersections(vec);
      if (temp_distances.size() > 0)
      	boundaries.add_intersections(vec.pt.t, sza_dimension,
      				     isza+1, sza_boundaries[isza+1], temp_distances);
      // ray_cone(vec.GT_ray3, sza_boundary_cones[isza]);
      // if (ray_cone.intersect)
      // 	boundaries.add_intersections(vec.pt.t, 1,
      // 				     isza+1, sza_boundaries[isza+1], ray_cone.distances);
    }

    //sort the list of intersections by distance & trim
    boundaries.sort();
    boundaries.propagate_indices();
    boundaries.assign_voxel_indices(this, &spherical_azimuthally_symmetric_grid::indices_to_voxel);
    boundaries.trim();
    boundaries.check(n_boundaries,n_voxels);
    
    if (save_intersections)
      saver.append_intersections(vec,boundaries);

    boundary_intersection_stepper stepper(vec,
					  boundaries,
					  n_emissions);

    

    return stepper;
  }



  
  VectorXd sza_slice(VectorXd quantity, int i_sza) {
    VectorXd ret;
    vector<int> indices(n_dimensions,i_sza);
    
    ret.resize(n_radial_boundaries-1);
    for (int i=0;i<n_radial_boundaries-1;i++) {
      indices[0]=i;
      ret(i) = quantity(indices_to_voxel(indices));
    }
    
    return ret;
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
       	
	VectorXd sza_boundaries_write_out = Eigen::Map<VectorXd>(sza_boundaries.data(),
								 sza_boundaries.size());

	file << "sza boundaries [rad]: " << sza_boundaries_write_out.transpose() << "\n\n";
	
	VectorXd sza_pts_write_out = Eigen::Map<VectorXd>(pts_sza.data(),
							  pts_sza.size());

	file << "pts sza [rad]: " << sza_pts_write_out.transpose() << "\n\n";
       	
	for (int i=0; i<n_emissions; i++) {
	  file << "For " << emission_names[i] << "\n";
	  for (unsigned int j=0; j<pts_sza.size(); j++) {
	    file << "  For SZA = " << pts_sza[j] << ": \n" 
		 << "    Species density [cm-3]: "
		 <<      sza_slice(species_density[i],j).transpose() << "\n"

		 << "    Species single scattering tau: " 
		 <<	 sza_slice(tau_species_single_scattering[i],j).transpose() << "\n"

		 << "    Species cross section [cm2]: " 
		 <<      sza_slice(species_sigma[i],j).transpose() << "\n"

		 << "    Absorber density [cm-3]: " 
		 <<      sza_slice(absorber_density[i],j).transpose() << "\n"

		 << "    Absorber single scattering tau: " 
		 <<	 sza_slice(tau_absorber_single_scattering[i],j).transpose() << "\n"

		 << "    Source function: " 
		 <<      sza_slice(sourcefn[i],j).transpose() << "\n\n";
	  }
	}
      }
    file.close();
  }

};




#endif
