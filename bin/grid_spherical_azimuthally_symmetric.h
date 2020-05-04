// grid_spherical_azimuthally_symmetric.h -- spherical grid with symmetry around the Mars-Sun line

#ifndef __grid_spherical_azimuthally_symmetric
#define __grid_spherical_azimuthally_symmetric

#include "cuda_compatibility.h"
#include "grid.h"
#include "coordinate_generation.h"
#include "boundaries.h"
#include "atmosphere.h"
#include "interp.h"
#include "intersections.h"
#include <type_traits> 
#include <fstream>
#include <cmath>
#include <string>

//convert all vectors to arrays via templates

struct spherical_azimuthally_symmetric_grid : grid<2> //this is a 2D grid 
{

  static const int r_dimension = 0;
  int n_radial_boundaries;
  int rmethod;
  static const int rmethod_altitude = 0;
  static const int rmethod_lognH = 1;
    
  vector<double> radial_boundaries;
  vector<double> pts_radii;
  vector<double> log_pts_radii;
  vector<sphere> radial_boundary_spheres;

  static const int sza_dimension = 1;
  int n_sza_boundaries;
  int szamethod;
  static const int szamethod_uniform = 0;
  static const int szamethod_uniform_cos = 1;

  vector<double> sza_boundaries;
  vector<double> pts_sza;
  vector<cone> sza_boundary_cones;

  int n_boundaries[n_dimensions];
  // intersection_writer saver;
  // bool save_intersections;

  int n_theta;
  vector<double> ray_theta;
  int n_phi;
  vector<double> ray_phi;

   
  spherical_azimuthally_symmetric_grid() {
    sun_direction = {0.,0.,1.};
    rmethod = rmethod_altitude;
    szamethod = szamethod_uniform;

    // save_intersections = false;
    // saver.fname = "intersections.dat";
  }

#ifdef __CUDACC__
  void copy_to_cuda(spherical_azimuthally_symmetric_grid *d_ptr) {
    //declare, allocate, and copy all of the subelements
    
  }
  void cuda_free() {
    //free the pointers? not sure if they're still in scope
    
  }
#endif
  
  void setup_voxels(int n_radial_boundariess,
		    int n_sza_boundariess,
		    atmosphere &atm)    
  {
    n_radial_boundaries = n_radial_boundariess;
    
    rmin = atm.rmin;
    rmax = atm.rmax;

    assert((rmethod == rmethod_altitude || rmethod == rmethod_lognH)
	   && "rmethod must match a defined radial points method");
    // don't define a tau radial points method; tau < 0.1 is
    // important and max(tau) > 10; this leads to many required
    // gridpoints
    if (rmethod == rmethod_altitude)
      get_radial_log_linear_points(radial_boundaries, n_radial_boundaries,
				   atm.rmin, atm.rexo, atm.rmax);
    if (rmethod == rmethod_lognH) {
      double lognH_max = log(atm.nH(atm.rmin));
      double lognH_min = log(atm.nH(atm.rmax));
      double lognH_step = (lognH_max-lognH_min)/(n_radial_boundaries-1.);
      
      for(int i=0;i<n_radial_boundaries;i++) {
	double nH_target=exp(lognH_max-i*lognH_step);
	radial_boundaries.push_back(atm.r_from_nH(nH_target));
      }
    }
    
    pts_radii.resize(n_radial_boundaries-1);
    log_pts_radii.resize(n_radial_boundaries-1);
    for (int i=0; i<n_radial_boundaries-1; i++) {
      pts_radii[i]=sqrt(radial_boundaries[i]*radial_boundaries[i+1]);
      log_pts_radii[i]=log(pts_radii[i]);
    }

    radial_boundary_spheres.resize(n_radial_boundaries);
    for (int i=0; i<n_radial_boundaries; i++) 
      radial_boundary_spheres[i].set_radius(radial_boundaries[i]);



    n_sza_boundaries = n_sza_boundariess;

    assert((szamethod == szamethod_uniform || szamethod == szamethod_uniform_cos)
	   && "szamethod must match a defined sza points method");
    if (szamethod == szamethod_uniform) {
      double sza_spacing = pi / (n_sza_boundaries - 2.);
      sza_boundaries.resize(n_sza_boundaries);
      for (int i=0;i<n_sza_boundaries;i++) {
	sza_boundaries[i]=(i-0.5)*sza_spacing;
      }
    }
    if (szamethod == szamethod_uniform_cos) {
      double cos_sza_spacing = 2.0 / (n_sza_boundaries - 2.);
      using std::acos;
      
      sza_boundaries.resize(n_sza_boundaries);
      sza_boundaries[0] = -acos(1.0-0.5*cos_sza_spacing);
      for (int i=1;i<n_sza_boundaries-1;i++) {
	sza_boundaries[i]=acos(1.0-(i-0.5)*cos_sza_spacing);
      }
      sza_boundaries[n_sza_boundaries-1] = pi + acos(1.0-0.5*cos_sza_spacing);
    }
  
    pts_sza.resize(n_sza_boundaries-1);
    for (unsigned int i=0;i<pts_sza.size();i++) {
      pts_sza[i]=0.5*(sza_boundaries[i] + sza_boundaries[i+1]);
    }
    
    sza_boundary_cones.resize(n_sza_boundaries-2);
    for (int i=0;i<n_sza_boundaries-2;i++)
      sza_boundary_cones[i].set_angle(sza_boundaries[i+1]);


    n_voxels = pts_radii.size()*pts_sza.size();
    n_boundaries[r_dimension] = n_radial_boundaries-1;
    n_boundaries[sza_dimension] = n_sza_boundaries-1;
    
    pts.resize(n_voxels);
    int ivoxel;
    for (unsigned int i=0; i<pts_radii.size(); i++) {
      for (unsigned int j=0;j<pts_sza.size();j++) {
	ivoxel = i*(n_sza_boundaries-1) + j;
	
	pts[ivoxel].rtp(pts_radii[i],pts_sza[j],0.);
	pts[ivoxel].set_voxel_index(ivoxel);
      }
    }

    // if (save_intersections)
    //   saver.save_coordinates(radial_boundaries,sza_boundaries);
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
  }

  CUDA_CALLABLE_MEMBER 
  void indices_to_voxel(const int &r_idx, const int &sza_idx, int & vox_idx) const {
    if ((r_idx   < 0) || (  r_idx > (int) pts_radii.size()-1) ||
	(sza_idx < 0) || (sza_idx > (int) pts_sza.size()-1))
      vox_idx = -1;
    else
      vox_idx = r_idx*(n_sza_boundaries-1)+sza_idx;
  }
  CUDA_CALLABLE_MEMBER 
  void indices_to_voxel(const int (&indices)[n_dimensions], int & vox_idx) const {
    indices_to_voxel(indices[r_dimension], indices[sza_dimension], vox_idx);
  }

  CUDA_CALLABLE_MEMBER 
  void voxel_to_indices(const int &i_voxel, int (&v)[n_dimensions]) const {
    if ((i_voxel < 0) || (i_voxel > n_voxels-1)) {
      v[r_dimension]=-1;
      v[sza_dimension]=-1;
    } else {
      v[r_dimension]=i_voxel / (n_sza_boundaries-1);
      v[sza_dimension]=i_voxel % (n_sza_boundaries-1);
    }
  }


  template <class V>
  CUDA_CALLABLE_MEMBER 
  int find_coordinate_index(const double &pt_coord, const V &boundaries) const {
    unsigned int i;
    
    for (i=0;i<boundaries.size();i++)
      if (pt_coord < boundaries[i])
	break;

    i--;

    assert((boundaries.back()<=pt_coord ||
	    pt_coord < boundaries[0] ||
	    (boundaries[i]<=pt_coord &&pt_coord<boundaries[i+1]))
	   && "we have found the appropriate point index");
    
    return i;
  }

  CUDA_CALLABLE_MEMBER 
  void point_to_indices(const atmo_point &pt, int (&indices)[n_dimensions]) const {
    indices[r_dimension] = find_coordinate_index(pt.r, radial_boundaries);
    indices[sza_dimension] = find_coordinate_index(pt.t, sza_boundaries);
  }
  


  CUDA_CALLABLE_MEMBER
  void ray_voxel_intersections(const atmo_vector &vec, boundary_intersection_stepper<n_dimensions> &stepper) const {
    
    boundary_set<n_dimensions> boundaries(2*n_radial_boundaries+n_sza_boundaries);
    boundaries.reset();
    
    //define the origin
    boundary<n_dimensions> origin;
    if (vec.pt.i_voxel == -1) {
      point_to_indices(vec.pt,origin.entering_indices);
      indices_to_voxel(origin.entering_indices,origin.entering);
    } else {
      origin.entering = vec.pt.i_voxel;
      voxel_to_indices(origin.entering,origin.entering_indices);
    }
    origin.distance = 0.0;
    boundaries.append(origin);

    //do the intersections for each coordinate
    int n_hits = 0;
    double temp_distances[2] = {-1,-1};
    for (unsigned int ir=0;ir<radial_boundary_spheres.size();ir++) {
      radial_boundary_spheres[ir].intersections(vec, temp_distances, n_hits);
      boundaries.add_intersections(vec.pt.r, r_dimension,
				   ir, radial_boundaries[ir],
				   temp_distances, n_hits);
    }

    for (unsigned int isza=0;isza<sza_boundary_cones.size();isza++) {
      sza_boundary_cones[isza].intersections(vec, temp_distances, n_hits);
      boundaries.add_intersections(vec.pt.t, sza_dimension,
				   isza+1, sza_boundaries[isza+1],
				   temp_distances, n_hits);
    }

    //sort the list of intersections by distance & trim
    boundaries.sort();
    boundaries.propagate_indices();
    boundaries.assign_voxel_indices(this);
    boundaries.trim();
    assert(boundaries.check(n_boundaries,n_voxels) && "boundary checks must pass");
    
    // if (save_intersections)
    //   saver.append_intersections(vec,boundaries);

    stepper = boundary_intersection_stepper<n_dimensions>(vec, boundaries);
  }

  CUDA_CALLABLE_MEMBER 
  void interp_weights(const int &ivoxel, const atmo_point &pt,
		      int (&indices)[n_interp_points], double (&weights)[n_interp_points]) const {
    // n_interp_points = 4, because this is a 2d grid with linear interpolation

    int coord_indices[n_dimensions];
    voxel_to_indices(ivoxel, coord_indices);
    int r_idx;
    r_idx = coord_indices[r_dimension];
    int sza_idx;
    sza_idx = coord_indices[sza_dimension];

    assert(radial_boundaries[r_idx] <= pt.r &&
	   pt.r <= radial_boundaries[r_idx+1]
	   && "pt must be in identified voxel.");
    assert(sza_boundaries[sza_idx] <= pt.t &&
	   pt.t <= sza_boundaries[sza_idx+1] &&
	   "pt must be in identified voxel.");

    int r_lower_pt_idx, r_upper_pt_idx;
    double r_wt;
    if (r_idx == 0 && pt.r <= pts_radii[0]) {
      //we are below the lowest radial point in the source function grid
      r_lower_pt_idx=r_upper_pt_idx=0;
      r_wt=1.0;
    } else if (r_idx == n_radial_boundaries-2 &&  pts_radii.back() <= pt.r) {
      //we are above the highest radial point in the source function grid
      r_lower_pt_idx=r_upper_pt_idx=n_radial_boundaries-2;
      r_wt=0.0;
    } else {
      //we are inside the radial grid

      //pts at which interp quanities are defined are offset from the
      //grid boundaries. Figure out whether to go up or down
      if (pt.r < pts_radii[r_idx]) {
	r_lower_pt_idx = r_idx - 1;
	r_upper_pt_idx = r_lower_pt_idx + 1;
      } else {
	r_lower_pt_idx = r_idx;
	r_upper_pt_idx = r_lower_pt_idx + 1;
      }

      assert(r_lower_pt_idx >= 0 && r_upper_pt_idx < n_radial_boundaries-1 && "interpolation points must lie on grid.");
      
      r_wt = (log(pt.r) - log_pts_radii[r_lower_pt_idx])/(log_pts_radii[r_upper_pt_idx]-log_pts_radii[r_lower_pt_idx]);
    }

    int sza_lower_pt_idx, sza_upper_pt_idx;
    double sza_wt;
    //we are always inside the SZA grid
    assert(pts_sza[0] <= pt.t && pt.t < pts_sza.back() && "pt must be inside SZA grid.");
    //pts at which interp quanities are defined are offset from the
    //grid boundaries. Figure out whether to go up or down
    if (pt.t < pts_sza[sza_idx]) {
      sza_lower_pt_idx = sza_idx - 1;
      sza_upper_pt_idx = sza_lower_pt_idx + 1;
    } else {
      sza_lower_pt_idx = sza_idx;
      sza_upper_pt_idx = sza_lower_pt_idx + 1;
    }
    sza_wt = (pt.t-pts_sza[sza_lower_pt_idx])/(pts_sza[sza_upper_pt_idx]-pts_sza[sza_lower_pt_idx]);

    
    indices_to_voxel(r_lower_pt_idx, sza_lower_pt_idx, indices[0]);
    weights[0] =  (1.-r_wt)   *   (1-sza_wt)    ;
    
    indices_to_voxel(r_upper_pt_idx, sza_lower_pt_idx, indices[1]);
    weights[1] =      r_wt    *   (1-sza_wt)    ;
    
    indices_to_voxel(r_lower_pt_idx, sza_upper_pt_idx, indices[2]);
    weights[2] =  (1.-r_wt)   *      sza_wt     ;
    
    indices_to_voxel(r_upper_pt_idx, sza_upper_pt_idx, indices[3]);
    weights[3] =      r_wt    *      sza_wt     ;

    assert(std::abs(weights[0]+weights[1]+weights[2]+weights[3] - 1.0) < 1e-5
	   && "interpolation weights must sum to 1.");
   }

  VectorXd sza_slice(VectorXd quantity, int i_sza) const {
    VectorXd ret;
    int indices[n_dimensions];
    indices[sza_dimension] = i_sza;
    
    ret.resize(n_radial_boundaries-1);
    for (int i=0;i<n_radial_boundaries-1;i++) {
      indices[r_dimension]=i;
      int voxel;
      indices_to_voxel(indices, voxel);
      ret(i) = quantity(voxel);
    }
    
    return ret;
  }

  void save_S(const string fname, const vector<emission> &emissions) const {
    std::ofstream file(fname.c_str());
    if (file.is_open())
      {

	VectorXd r_boundaries_write_out = Eigen::Map<const VectorXd>(radial_boundaries.data(),
								     radial_boundaries.size());

	file << "radial boundaries [cm]: " << r_boundaries_write_out.transpose() << "\n\n";

	VectorXd r_pts_write_out = Eigen::Map<const VectorXd>(pts_radii.data(),
							      pts_radii.size());

	file << "pts radii [cm]: " << r_pts_write_out.transpose() << "\n\n";
       	
	VectorXd sza_boundaries_write_out = Eigen::Map<const VectorXd>(sza_boundaries.data(),
								       sza_boundaries.size());

	file << "sza boundaries [rad]: " << sza_boundaries_write_out.transpose() << "\n\n";
	
	VectorXd sza_pts_write_out = Eigen::Map<const VectorXd>(pts_sza.data(),
								pts_sza.size());

	file << "pts sza [rad]: " << sza_pts_write_out.transpose() << "\n\n";
       	
	for (auto&& emiss: emissions) {
	  file << "For " << emiss.name << "\n";
	  for (unsigned int j=0; j<pts_sza.size(); j++) {
	    file << "  For SZA = " << pts_sza[j] << ": \n" 
		 << "    Species density [cm-3]: "
		 <<      sza_slice(emiss.species_density,j).transpose() << "\n"

		 << "    Species single scattering tau: " 
		 <<	 sza_slice(emiss.tau_species_single_scattering,j).transpose() << "\n"

		 << "    Species cross section [cm2]: " 
		 <<      sza_slice(emiss.species_sigma,j).transpose() << "\n"

		 << "    Absorber density [cm-3]: " 
		 <<      sza_slice(emiss.absorber_density,j).transpose() << "\n"

		 << "    Absorber single scattering tau: " 
		 <<	 sza_slice(emiss.tau_absorber_single_scattering,j).transpose() << "\n"

		 << "    Source function: " 
		 <<      sza_slice(emiss.sourcefn,j).transpose() << "\n\n";
	  }
	}
      }
    file.close();
  }

};

#endif
