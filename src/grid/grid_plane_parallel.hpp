// grid_plane_parallel.h -- plane parallel RT grid

#ifndef __grid_plane_parallel
#define __grid_plane_parallel

#include "Real.hpp"
#include "constants.hpp"
#include "grid.hpp"
#include "coordinate_generation.hpp"
#include "boundaries.hpp"
#include "atm/atmosphere_base.hpp"
#include <fstream>
#include <cmath>
#include <string>
#include "intersections.hpp"

template <int N_RADIAL_BOUNDARIES, int N_RAYS_THETA>
struct plane_parallel_grid : grid<1,//N_DIMENSIONS, this is a 1D grid
			   N_RADIAL_BOUNDARIES-1,//N_VOXELS
			   N_RAYS_THETA,//N_RAYS
			   N_RADIAL_BOUNDARIES,//N_MAX_INTERSECTIONS
			   plane_parallel_grid<N_RADIAL_BOUNDARIES, N_RAYS_THETA>>
{

  using parent_grid = grid<1,//N_DIMENSIONS
			   N_RADIAL_BOUNDARIES-1,//N_VOXELS
			   N_RAYS_THETA,//N_RAYS
			   N_RADIAL_BOUNDARIES,//N_MAX_INTERSECTIONS
			   plane_parallel_grid<N_RADIAL_BOUNDARIES, N_RAYS_THETA>>;

  int rmethod;
  static const int rmethod_altitude = 0;
  static const int rmethod_log_n_species= 1;
  static const int rmethod_log_n_species_tau_absorber = 2;
  static const int rmethod_log_n_species_int = 3;

  static const int n_radial_boundaries = N_RADIAL_BOUNDARIES;
  Real radial_boundaries[n_radial_boundaries];
  Real pts_radii[n_radial_boundaries-1];
  plane radial_boundary_planes[n_radial_boundaries];
  
  plane_parallel_grid() {
    this->n_pts[0] = n_radial_boundaries-1;
 }
  
  void setup_voxels(const atmosphere &atm) {
    this->rmin = atm.rmin;
    this->rmax = atm.rmax;

    assert((rmethod == rmethod_altitude
	    || rmethod == rmethod_log_n_species
	    || rmethod == rmethod_log_n_species_tau_absorber
	    || rmethod == rmethod_log_n_species_int)
	   && "rmethod must match a defined radial points method");
    // don't define a tau radial points method; tau < 0.1 is
    // important and max(tau) > 10; this leads to many required
    // gridpoints
    if (rmethod == rmethod_altitude) {
      vector<Real> radial_boundaries_vector;
      get_radial_log_linear_points(radial_boundaries_vector, n_radial_boundaries,
				   atm.rmin, atm.rexo, atm.rmax);
      for (int i=0;i<n_radial_boundaries;i++)
	radial_boundaries[i] = radial_boundaries_vector[i];
    }
    if (rmethod == rmethod_log_n_species) {
      Real log_n_species_max = log(atm.n_species(atm.rmin));
      Real log_n_species_min = log(atm.n_species(atm.rmax));
      Real log_n_species_step = (log_n_species_max-log_n_species_min)/(n_radial_boundaries-1.);
      
      for(int i=0;i<n_radial_boundaries;i++) {
	Real n_species_target=exp(log_n_species_max-i*log_n_species_step);
	radial_boundaries[i]=atm.r_from_n_species(n_species_target);
      }
    }
    if (rmethod == rmethod_log_n_species_tau_absorber) {
      // warn the user about this method
      std::cout << "Warning: using CO2 cross section at Lyman alpha to define grid points." << std::endl;
      
      //construct the integral of log_n_species * exp(-tau_absorber)
      const int n_int_steps = 1000;
      const doubReal deriv_step = 0.0001;
      
      const doubReal log_r_int_step = (
				     (log(atm.rmax-rMars) - log(atm.rmin-rMars))
				     /
				     (n_int_steps - 1.)
				     );
      
      const doubReal abs_xsec = CO2_lyman_alpha_absorption_cross_section;
      
      vector<doubReal> log_r_int;
      vector<doubReal> n_absorber_int;
      vector<doubReal> int_lognH_exp_tauabs;
      
      log_r_int.push_back(log((atm.rmax-rMars)*(1-ATMEPS)));
      n_absorber_int.push_back(0);
      int_lognH_exp_tauabs.push_back(0);
      for (int i_int=1; i_int<n_int_steps; i_int++) {
	log_r_int.push_back(log_r_int[0]-i_int*log_r_int_step);
	if (exp(log_r_int.back()) < (atm.rmin-rMars))
	  log_r_int.back() = log((atm.rmin-rMars)*(1+ATMEPS));
	
	//unscaled quantities
	doubReal r0 = exp(log_r_int[i_int-1])+rMars;
	doubReal r1 = exp(log_r_int[i_int])+rMars;
	doubReal dr = r0-r1;

	doubReal n_absorber_diff = 0.5*(atm.n_absorber(r0) + atm.n_absorber(r1))*dr;
	n_absorber_int.push_back(n_absorber_int.back() + n_absorber_diff );

	doubReal r0up = r0*(1+deriv_step);
	doubReal r0dn = r0*(1-deriv_step);
	doubReal dr0 = r0up-r0dn;

	doubReal r1up = r1*(1+deriv_step);
	doubReal r1dn = r1*(1-deriv_step);
	r1dn = r1dn < atm.rmin ? atm.rmin : r1dn;
	doubReal dr1 = r1up-r1dn;
	
	doubReal int_diff = 0.5*(
			       (
				log(atm.n_species(r0dn))
				-
				log(atm.n_species(r0up))
				)
			       /dr0
			       *exp(-abs_xsec*n_absorber_int[i_int-1])
			       +
			       (
				log(atm.n_species(r1dn))
				-
				log(atm.n_species(r1up))
				)
			       /dr1
			       *exp(-abs_xsec*n_absorber_int[i_int])
			       )*dr;

	int_lognH_exp_tauabs.push_back(int_lognH_exp_tauabs.back() + int_diff );
      }

      // now subdivide the integral and find the appropriate grid points
      doubReal int_step = int_lognH_exp_tauabs.back() / (n_radial_boundaries - 1);
      doubReal target = int_step;
      radial_boundaries[n_radial_boundaries-1] = atm.rmax;
      int boundary = n_radial_boundaries-2;
      for (int i_int=1; i_int<n_int_steps; i_int++) {
	if (int_lognH_exp_tauabs[i_int] > target) {
	  radial_boundaries[boundary] = exp(log_r_int[i_int])+rMars;
	  target += int_step;
	  boundary--;
	}
      }
      assert(boundary==0 && "we must have found all boundaries");
      radial_boundaries[0] = atm.rmin;
    }
    if (rmethod == rmethod_log_n_species_int) {
      const atmosphere_average_1d *atm_avg = dynamic_cast<const atmosphere_average_1d*>(&atm);
      assert((atm_avg != NULL) && "This radial points method only works with a class derived from atmosphere_average_1d");
      
      const doubReal logtaumax = std::log(atm_avg->n_species_int.back());
      const doubReal logtaumin = std::log(atm_avg->n_species_int[1]);
      const doubReal logtaumax_step = (logtaumax-logtaumin)/ (n_radial_boundaries-1);

      doubReal target = logtaumax_step+logtaumin;
      radial_boundaries[n_radial_boundaries-1] = atm_avg->rmax;
      int boundary = n_radial_boundaries-2;
      for (int i_int=1; i_int<(int)atm_avg->n_species_int.size(); i_int++) {
	while (std::log(atm_avg->n_species_int[i_int]) > target && boundary>-1) {
	  doubReal upper = log(atm_avg->n_species_int[i_int-1]);
	  doubReal lower = log(atm_avg->n_species_int[i_int]);
	  if (!isfinite(upper))
	    upper = lower - 10;

	  doubReal frac = ((target - upper)
			 /
			 (lower - upper));
	  doubReal altfrac = (     frac *exp(atm_avg->log_r_int[i_int])
			    + (1-frac)*exp(atm_avg->log_r_int[i_int-1]));
	  
	  radial_boundaries[boundary] = altfrac*atm_avg->r_int_scale+rMars;
	  target += logtaumax_step;
	  boundary--;
	}
      }
      assert((boundary==0 || boundary==-1) && "we must have found all boundaries");
      radial_boundaries[0] = atm_avg->rmin;
    }

    for (int i=0; i<n_radial_boundaries-1; i++) {
      this->voxels[i].rbounds[0] = radial_boundaries[i];
      this->voxels[i].rbounds[1] = radial_boundaries[i+1];
      this->voxels[i].tbounds[0] = 0;
      this->voxels[i].tbounds[1] = pi;
      this->voxels[i].pbounds[0] = 0;
      this->voxels[i].pbounds[1] = 2*pi;
      this->voxels[i].i_voxel = i;
      // this->voxels[i].init = true;

      pts_radii[i]=sqrt(radial_boundaries[i]*radial_boundaries[i+1]);
      // Real frac = 0.5;
      // pts_radii[i]=(frac*radial_boundaries[i]+(1-frac)*radial_boundaries[i+1]);
      this->voxels[i].pt.xyz(0.,0.,pts_radii[i]);
      this->voxels[i].pt.set_voxel_index(i);
    }

    for (int i=0; i<n_radial_boundaries; i++) 
      radial_boundary_planes[i].z = radial_boundaries[i];
  }
  
  void setup_rays() {
    //in a plane-parallel grid only the angle with the vertical matters
    vector<Real> ray_theta;
    vector<Real> ray_weights;

    gauss_quadrature_points(ray_theta,ray_weights,0,pi,parent_grid::n_rays);
    for (int i=0;i<parent_grid::n_rays;i++)
      ray_weights[i]*=std::sin(ray_theta[i]);
    
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
  void voxel_to_indices(const int i_voxel, int (&indices)[parent_grid::n_dimensions]) const {
    if ((i_voxel < 0) || (i_voxel > parent_grid::n_voxels-1))
      indices[0]=-1;
    else
      indices[0]=i_voxel;
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

  CUDA_CALLABLE_MEMBER 
  void interp_weights(__attribute__((unused)) const int &ivoxel, __attribute__((unused)) const atmo_point &ptt,
		      __attribute__((unused)) int (&indices)[parent_grid::n_interp_points],
		      __attribute__((unused)) Real (&weights)[parent_grid::n_interp_points],
		      __attribute__((unused)) int (&indices_1d)[2*parent_grid::n_dimensions],
		      __attribute__((unused)) Real (&weights_1d)[parent_grid::n_dimensions]) const {
    assert(false && "interp_weights not implemented in grid_plane_parallel");
  }

  template<typename E>
  void save_S(const string &fname, const E* const *emissions, const int n_emissions) const {
    std::ofstream file(fname.c_str());
    if (file.is_open())
      {

	VectorX r_boundaries_write_out = Eigen::Map<const VectorX>(radial_boundaries,
								   n_radial_boundaries);
	
	file << "radial boundaries [cm]: " << r_boundaries_write_out.transpose() << "\n\n";
	
	VectorX r_pts_write_out = Eigen::Map<const VectorX>(pts_radii,
							    n_radial_boundaries-1);

	file << "pts radii [cm]: " << r_pts_write_out.transpose() << "\n\n";

	for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	  file << "For " << emissions[i_emission]->name() << ",\n";
	  emissions[i_emission]->save(file,[](VectorX vec, __attribute__((unused)) int i) {return vec;},0);
        }
      }
  }
};

#endif
