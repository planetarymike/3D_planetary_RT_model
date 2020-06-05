//observation_fit.hpp -- routines to fit an atmosphere observation

#ifndef __OBSERVATION_FIT_H
#define __OBSERVATION_FIT_H

#include "Real.hpp"
#include "observation.hpp"
#include "atm/temperature.hpp"
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
    
  observation<n_emissions> obs;

  krasnopolsky_temperature temp;
  const Real CO2_exobase_density = 2e8;//cm-3
  chamb_diff_1d atm;

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

public:
  observation_fit();

  void add_observation(std::vector<Vector3> MSO_locations,
		       std::vector<Vector3> MSO_directions);

  void set_g_factor(Real &g);

  void generate_source_function(Real nHexo, Real Texo);
  
  std::vector<Real> brightness();  
};

//might be needed to instantiate template members
//#include "observation_fit.cpp"
//observation_fit hello;

#endif
