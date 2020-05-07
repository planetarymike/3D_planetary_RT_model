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
template <int N_RADIAL_BOUNDARIES, int N_SZA_BOUNDARIES, int N_RAY_THETA, int N_RAY_PHI>
struct spherical_azimuthally_symmetric_grid : grid<2, //this is a 2D grid 
						   (N_RADIAL_BOUNDARIES-1)*(N_SZA_BOUNDARIES-1),
						   N_RAY_THETA*N_RAY_PHI,
						   2*N_RADIAL_BOUNDARIES+N_SZA_BOUNDARIES>
{

  typedef grid<2,(N_RADIAL_BOUNDARIES-1)*(N_SZA_BOUNDARIES-1),
	       N_RAY_THETA*N_RAY_PHI,
	       2*N_RADIAL_BOUNDARIES+N_SZA_BOUNDARIES> parent_grid;

  static const int r_dimension = 0;
  static const int n_radial_boundaries = N_RADIAL_BOUNDARIES;
  int rmethod;
  static const int rmethod_altitude = 0;
  static const int rmethod_lognH = 1;
    
  double radial_boundaries[n_radial_boundaries];
  double pts_radii[n_radial_boundaries-1];
  double log_pts_radii[n_radial_boundaries-1];
  sphere radial_boundary_spheres[n_radial_boundaries];

  static const int sza_dimension = 1;
  static const int n_sza_boundaries = N_SZA_BOUNDARIES;
  int szamethod;
  static const int szamethod_uniform = 0;
  static const int szamethod_uniform_cos = 1;

  double sza_boundaries[n_sza_boundaries];
  double pts_sza[n_sza_boundaries-1];
  cone sza_boundary_cones[n_sza_boundaries-2];//there are extra
					      //boundaries outside of
					      //the range to put
					      //pts_sza=0,-1 on the
					      //planet-sun line

  int n_pts[parent_grid::n_dimensions] = {n_radial_boundaries-1,n_sza_boundaries-1};
  // intersection_writer saver;
  // bool save_intersections;

  static const int n_theta = N_RAY_THETA;
  double ray_theta[n_theta];
  static const int n_phi = N_RAY_PHI;
  double ray_phi[n_phi];

   
  spherical_azimuthally_symmetric_grid() {
    this->sun_direction = {0.,0.,1.};
    rmethod = rmethod_altitude;
    szamethod = szamethod_uniform;

    // save_intersections = false;
    // saver.fname = "intersections.dat";
  }

  void setup_voxels(atmosphere &atm) {
    this->rmin = atm.rmin;
    this->rmax = atm.rmax;
    
    assert((rmethod == rmethod_altitude || rmethod == rmethod_lognH)
	   && "rmethod must match a defined radial points method");
    // don't define a tau radial points method; tau < 0.1 is
    // important and max(tau) > 10; this leads to many required
    // gridpoints
    if (rmethod == rmethod_altitude) {
      vector<double> radial_boundaries_vector;
      get_radial_log_linear_points(radial_boundaries_vector, n_radial_boundaries,
				   atm.rmin, atm.rexo, atm.rmax);
      for (int i=0;i<n_radial_boundaries;i++)
	radial_boundaries[i] = radial_boundaries_vector[i];
    }
    if (rmethod == rmethod_lognH) {
      double lognH_max = log(atm.nH(atm.rmin));
      double lognH_min = log(atm.nH(atm.rmax));
      double lognH_step = (lognH_max-lognH_min)/(n_radial_boundaries-1.);
      
      for(int i=0;i<n_radial_boundaries;i++) {
	double nH_target=exp(lognH_max-i*lognH_step);
	radial_boundaries[i]=atm.r_from_nH(nH_target);
      }
    }
    
    for (int i=0; i<n_radial_boundaries-1; i++) {
      pts_radii[i]=sqrt(radial_boundaries[i]*radial_boundaries[i+1]);
      log_pts_radii[i]=log(pts_radii[i]);
    }

    for (int i=0; i<n_radial_boundaries; i++) 
      radial_boundary_spheres[i].set_radius(radial_boundaries[i]);



    assert((szamethod == szamethod_uniform || szamethod == szamethod_uniform_cos)
	   && "szamethod must match a defined sza points method");
    if (szamethod == szamethod_uniform) {
      double sza_spacing = pi / (n_sza_boundaries - 2.);
      for (int i=0;i<n_sza_boundaries;i++) {
	sza_boundaries[i]=(i-0.5)*sza_spacing;
      }
    }
    if (szamethod == szamethod_uniform_cos) {
      double cos_sza_spacing = 2.0 / (n_sza_boundaries - 2.);
      using std::acos;
      
      sza_boundaries[0] = -acos(1.0-0.5*cos_sza_spacing);
      for (int i=1;i<n_sza_boundaries-1;i++) {
	sza_boundaries[i]=acos(1.0-(i-0.5)*cos_sza_spacing);
      }
      sza_boundaries[n_sza_boundaries-1] = pi + acos(1.0-0.5*cos_sza_spacing);
    }
  
    for (unsigned int i=0;i<n_sza_boundaries-1;i++) {
      pts_sza[i]=0.5*(sza_boundaries[i] + sza_boundaries[i+1]);
    }
    
    for (int i=0;i<n_sza_boundaries-2;i++)
      sza_boundary_cones[i].set_angle(sza_boundaries[i+1]);


    int ivoxel;
    for (unsigned int i=0; i<n_radial_boundaries-1; i++) {
      for (unsigned int j=0;j<n_sza_boundaries-1;j++) {
	ivoxel = i*(n_sza_boundaries-1) + j;

	this->pts[ivoxel].rtp(pts_radii[i],pts_sza[j],0.);
	this->pts[ivoxel].set_voxel_index(ivoxel);
      }
    }

    // if (save_intersections)
    //   saver.save_coordinates(radial_boundaries,sza_boundaries);
  }

  void setup_rays()
  {
    vector<double> ray_theta_vector;
    vector<double> ray_weights_theta;
    gauss_quadrature_points(ray_theta_vector,ray_weights_theta,0,pi,n_theta);

    for (int i=0;i<n_theta;i++)
      ray_theta[i] = ray_theta_vector[i];
    
    double phi_spacing = 2*pi/n_phi;
    for (int i=0;i<n_phi;i++)
      ray_phi[i] = (i+0.5)*phi_spacing;
    
    int iray;
    for (int i=0;i<n_theta;i++) {
      for (int j=0;j<n_phi;j++) {
	iray = i * n_phi + j;
	this->rays[iray].tp(ray_theta[i],ray_phi[j]);
	this->rays[iray].set_ray_index(iray, ray_weights_theta[i], phi_spacing);
      }
    }
  }

  CUDA_CALLABLE_MEMBER 
  void indices_to_voxel(const int &r_idx, const int &sza_idx, int & vox_idx) const {
    if ((r_idx   < 0) || (  r_idx > (int) n_radial_boundaries-2) ||
	(sza_idx < 0) || (sza_idx > (int) n_sza_boundaries-2))
      vox_idx = -1;
    else
      vox_idx = r_idx*(n_sza_boundaries-1)+sza_idx;
  }
  CUDA_CALLABLE_MEMBER 
  void indices_to_voxel(const int *indices, int & vox_idx) const {
    indices_to_voxel(indices[r_dimension], indices[sza_dimension], vox_idx);
  }

  CUDA_CALLABLE_MEMBER 
  void voxel_to_indices(const int &i_voxel, int *v) const {
    if ((i_voxel < 0) || (i_voxel > parent_grid::n_voxels-1)) {
      v[r_dimension]=-1;
      v[sza_dimension]=-1;
    } else {
      v[r_dimension]=i_voxel / (n_sza_boundaries-1);
      v[sza_dimension]=i_voxel % (n_sza_boundaries-1);
    }
  }


  template <class V>
  CUDA_CALLABLE_MEMBER 
  int find_coordinate_index(const double &pt_coord, const V &boundaries, int n_boundaries) const {
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
  void point_to_indices(const atmo_point &pt, int (&indices)[parent_grid::n_dimensions]) const {
    indices[r_dimension] = find_coordinate_index(pt.r, radial_boundaries, n_radial_boundaries);
    indices[sza_dimension] = find_coordinate_index(pt.t, sza_boundaries, n_sza_boundaries);
  }
  

  CUDA_CALLABLE_MEMBER
  void ray_voxel_intersections(const atmo_vector &vec,
			       boundary_intersection_stepper<parent_grid::n_dimensions,parent_grid::n_max_intersections> &stepper) const {
    
    stepper.vec = vec;
    stepper.boundaries.reset();
    
    //define the origin
    boundary<this->n_dimensions> origin;
    if (vec.pt.i_voxel == -1) {
      point_to_indices(vec.pt,origin.entering_indices);
      indices_to_voxel(origin.entering_indices,origin.entering);
    } else {
      origin.entering = vec.pt.i_voxel;
      voxel_to_indices(origin.entering,origin.entering_indices);
    }
    origin.distance = 0.0;
    stepper.boundaries.append(origin);

    //do the intersections for each coordinate
    int n_hits = 0;
    double temp_distances[2] = {-1,-1};
    for (unsigned int ir=0;ir<n_radial_boundaries;ir++) {
      radial_boundary_spheres[ir].intersections(vec, temp_distances, n_hits);
      stepper.boundaries.add_intersections(vec.pt.r, r_dimension,
					   ir, radial_boundaries[ir],
					   temp_distances, n_hits);
    }

    for (unsigned int isza=0;isza<n_sza_boundaries-2;isza++) {
      sza_boundary_cones[isza].intersections(vec, temp_distances, n_hits);
      stepper.boundaries.add_intersections(vec.pt.t, sza_dimension,
					   isza+1, sza_boundaries[isza+1],
					   temp_distances, n_hits);
    }

    //sort the list of intersections by distance & trim
    stepper.boundaries.sort();
    stepper.boundaries.propagate_indices();
    stepper.boundaries.assign_voxel_indices(this);
    stepper.boundaries.trim();
    int tnvoxels = this->n_voxels;
    assert(stepper.boundaries.check(n_pts, tnvoxels) && "boundary checks must pass");
    
    // if (save_intersections)
    //   saver.append_intersections(vec,stepper.boundaries);

    stepper.init_stepper();
  }

  CUDA_CALLABLE_MEMBER 
  void interp_weights(const int &ivoxel, const atmo_point &pt,
		      int (&indices)[parent_grid::n_interp_points], double (&weights)[parent_grid::n_interp_points]) const {
    // n_interp_points = 4, because this is a 2d grid with linear interpolation

    int coord_indices[parent_grid::n_dimensions];
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
    } else if (r_idx == n_radial_boundaries-2 &&  pts_radii[n_radial_boundaries-2] <= pt.r) {
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
    assert(pts_sza[0] <= pt.t && pt.t < pts_sza[n_sza_boundaries-2] && "pt must be inside SZA grid.");
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
    int indices[parent_grid::n_dimensions];
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

  void save_S(const string fname, const emission<parent_grid::n_voxels> *emissions, const int n_emissions) const {
    std::ofstream file(fname.c_str());
    if (file.is_open())
      {

	VectorXd r_boundaries_write_out = Eigen::Map<const VectorXd>(radial_boundaries,
								     n_radial_boundaries);

	file << "radial boundaries [cm]: " << r_boundaries_write_out.transpose() << "\n\n";

	VectorXd r_pts_write_out = Eigen::Map<const VectorXd>(pts_radii,
							      n_radial_boundaries-1);

	file << "pts radii [cm]: " << r_pts_write_out.transpose() << "\n\n";
       	
	VectorXd sza_boundaries_write_out = Eigen::Map<const VectorXd>(sza_boundaries,
								       n_sza_boundaries);

	file << "sza boundaries [rad]: " << sza_boundaries_write_out.transpose() << "\n\n";
	
	VectorXd sza_pts_write_out = Eigen::Map<const VectorXd>(pts_sza,
								n_sza_boundaries-1);

	file << "pts sza [rad]: " << sza_pts_write_out.transpose() << "\n\n";
       	
	for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	  file << "For " << emissions[i_emission].name << "\n";
	  for (unsigned int j=0; j<n_sza_boundaries-1; j++) {
	    file << "  For SZA = " << pts_sza[j] << ": \n" 
		 << "    Species density [cm-3]: "
		 <<      sza_slice(emissions[i_emission].species_density,j).transpose() << "\n"

		 << "    Species single scattering tau: " 
		 <<	 sza_slice(emissions[i_emission].tau_species_single_scattering,j).transpose() << "\n"

		 << "    Species cross section [cm2]: " 
		 <<      sza_slice(emissions[i_emission].species_sigma,j).transpose() << "\n"

		 << "    Absorber density [cm-3]: " 
		 <<      sza_slice(emissions[i_emission].absorber_density,j).transpose() << "\n"

		 << "    Absorber single scattering tau: " 
		 <<	 sza_slice(emissions[i_emission].tau_absorber_single_scattering,j).transpose() << "\n"

		 << "    Source function: " 
		 <<      sza_slice(emissions[i_emission].sourcefn,j).transpose() << "\n\n";
	  }
	}
      }
    file.close();
  }

};

#endif
