//fit_observation.h -- routines to fit an atmosphere observation

#ifndef __FIT_OBSERVATION_H
#define __FIT_OBSERVATION_H


#include "observation.h"
#include "atmosphere.h"
#include "RT_grid.h"
#include "RT_gpu.h"
#include "grid_spherical_azimuthally_symmetric.h"

struct observation_fit {
private:
  inline static const vector<string> emission_names = {"H Lyman alpha"};
    
  observation obs;

  krasnopolsky_temperature temp;
  const Real CO2_exobase_density = 2e8;//cm-3
  chamb_diff_1d atm;
  holstein_approx hol;

  static const int n_radial_boundaries = 40;
  static const int n_sza_boundaries = 20;/*20 for 10 deg increments with szamethod_uniform*/
  static const int n_rays_phi = 6;
  static const int n_rays_theta = 12;
  spherical_azimuthally_symmetric_grid<n_radial_boundaries,n_sza_boundaries,n_rays_phi,n_rays_theta> grid;
  RT_grid<2,typeof(grid),holstein_approx> RT(emission_names, grid, hol);

public:
  observation_fit()
    : obs(emission_names),
      atm(/*nHexo = */5e5, /*nCO2exo = */1e9, temp), //these dummy values can be overwritten later
      grid(emission_names, hol)
  {
    grid.rmethod = grid.rmethod_lognH;
    grid.szamethod = grid.szamethod_uniform_cos;
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
    RT.grid.setup_rays();

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
