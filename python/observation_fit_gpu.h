//observation_fit_gpu.h -- routines to fit an atmosphere observation

#ifndef __OBSERVATION_FIT_GPU_CU
#define __OBSERVATION_FIT_GPU_CU

#include "Real_is_float.h"
#include "observation.h"
#include "atmosphere.h"
#include "RT_grid.h"
#include "grid_spherical_azimuthally_symmetric.h"

struct observation_fit {
private:
  static const int n_emissions = 1;
  const string emission_names[n_emissions];// = {"H Lyman alpha"};
    
  observation<n_emissions> obs;

  krasnopolsky_temperature temp;
  const Real CO2_exobase_density = 2e8;//cm-3
  chamb_diff_1d atm;
  typedef holstein_approx influence_function;

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
	  influence_function> RT;

public:
  observation_fit()
    :
    emission_names({"H Lyman alpha"}),
    obs(emission_names),
    atm(/*nHexo = */5e5, /*nCO2exo = */1e9, temp), //these dummy values can be overwritten later
    RT(emission_names)
  {
    RT.grid.rmethod = RT.grid.rmethod_lognH;
    RT.grid.szamethod = RT.grid.szamethod_uniform_cos;
    RT.grid.setup_rays();
  }

  void add_observation(vector<vector<Real>> MSO_locations, vector<vector<Real>> MSO_directions) {
    obs.add_MSO_observation(MSO_locations,MSO_directions);
  }

  // void add_observation(Real* MSO_locations, Real* MSO_directions, int n_obs) {
  //   obs.add_MSO_observation(MSO_locations,MSO_directions,n_obs);
  // }

  void set_g_factor(Real &g) {
    obs.emission_g_factors[0] = g;
  }

  void generate_source_function(Real nHexo, Real Texo) {
    temp = krasnopolsky_temperature(Texo);
    atm = chamb_diff_1d(nHexo,CO2_exobase_density,temp);

    RT.grid.setup_voxels(atm);

    //update the RT grid values
    RT.define_emission("H Lyman alpha",
		       1.0,
		       atm,
		       &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lya,
		       &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lya);
    
    RT.generate_S();
  }
  
  vector<Real> brightness() {
    obs.reset_output();

    RT.brightness_gpu(obs);

    vector<Real> brightness;
    brightness.resize(obs.size());

    for (int i=0;i<obs.size();i++)
      brightness[i] = obs.brightness[i][0];

    return brightness;
  }

  
};

#endif
