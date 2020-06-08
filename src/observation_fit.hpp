//observation_fit.hpp -- routines to fit an atmosphere observation

#ifndef __OBSERVATION_FIT_H
#define __OBSERVATION_FIT_H

#include "Real.hpp"
#include "observation.hpp"
#include "atm/temperature.hpp"
#include "atm/chamberlain_exosphere.hpp"
#include "atm/chamb_diff_1d.hpp"
#include "RT_grid.hpp"
#include "grid_spherical_azimuthally_symmetric.hpp"

class observation_fit {
protected:
  static const int n_emissions = 1;
  const std::string emission_names[n_emissions];// = {"H Lyman alpha"};
					   // nvcc complains about
					   // inline definition, this
					   // needs to go in
					   // constructor

  static const int n_parameters = 2; // nH_exo and some T_exo type
  static const int n_pts_per_derivative = 2; // central difference

  static const int n_simulate_per_emission = n_parameters*n_pts_per_derivative + 1;
  static const int n_simulate = n_emissions*n_simulate_per_emission;
  std::string simulate_names[n_simulate];
  
  observation<n_emissions> obs;
  observation<n_simulate> obs_deriv;
  
  krasnopolsky_temperature temp;
  const Real CO2_exobase_density = 2e8;//cm-3
  chamb_diff_1d atm;//make sure to use the same exobase alt as in Tconv

  typedef holstein_approx influence_type;
  
  static const int n_radial_boundaries = 40;
  static const int n_sza_boundaries = 20;/*20 for 10 deg increments with szamethod_uniform*/
  static const int n_rays_phi = 6;
  static const int n_rays_theta = 12;
  typedef spherical_azimuthally_symmetric_grid<n_radial_boundaries,
					       n_sza_boundaries,
					       n_rays_phi,
					       n_rays_theta> grid_type;
  
  RT_grid<n_emissions,
	  grid_type,
	  influence_type> RT;

  RT_grid<n_simulate,
	  grid_type,
	  influence_type> RT_deriv;

public:
  observation_fit();

  Temp_converter Tconv;//also takes exobase alt argument

  void add_observation(const std::vector<vector<Real>> &MSO_locations,
		       const std::vector<vector<Real>> &MSO_directions);

  void add_observed_brightness(const std::vector<vector<Real>> &brightness,
			       const std::vector<vector<Real>> &sigma);
  
  void set_g_factor(Real &g);

  void generate_source_function(Real nHexo, Real Texo);
  void generate_source_function_effv(Real nHexo, Real effv_exo);
  void generate_source_function_lc(Real nHexo, Real lc_exo);
  
  std::vector<Real> brightness();  

  std::vector<Real> likelihood_and_derivatives(Real nHexo, Real Texo);
};

//might be needed to instantiate template members
//#include "observation_fit.cpp"
//observation_fit hello;

#endif
