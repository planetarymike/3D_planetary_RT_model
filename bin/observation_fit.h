//fit_observation.h -- routines to fit an atmosphere observation

#ifndef __FIT_OBSERVATION_H
#define __FIT_OBSERVATION_H


#include "observation.h"
#include "atmosphere.h"
#include "RT_spherical_azimuthally_symmetric.h"

struct observation_fit {
private:
  inline static const vector<string> emission_names = {"H Lyman alpha"};
    
  observation obs;

  krasnopolsky_temperature temp;
  const double CO2_exobase_density = 2e8;//cm-3
  chamb_diff_1d atm;
  holstein_approx hol;
  spherical_azimuthally_symmetric_grid grid;

public:
  observation_fit()
    : obs(emission_names),
      atm(/*nHexo = */5e5, /*nCO2exo = */1e9, temp), //these dummy values can be overwritten later
      grid(emission_names, hol)
  {
    grid.rmethod = grid.rmethod_lognH;
    grid.szamethod = grid.szamethod_uniform_cos;
  }

  void add_observation(vector<vector<double>> MSO_locations, vector<vector<double>> MSO_directions) {
    obs.add_MSO_observation(MSO_locations,MSO_directions);
  }

  // void add_observation(double* MSO_locations, double* MSO_directions, int n_obs) {
  //   obs.add_MSO_observation(MSO_locations,MSO_directions,n_obs);
  // }

  void set_g_factor(double &g) {
    obs.emission_g_factors[0] = g;
  }

  void generate_source_function(double nHexo, double Texo) {
    temp = krasnopolsky_temperature(Texo);
    atm = chamb_diff_1d(nHexo,CO2_exobase_density,temp);

    grid.setup_voxels(40, 20/*20 for 10 deg increments with sza_uniform*/, atm);
    grid.setup_rays(6, 12);

    //update the RT grid values
    grid.define_emission("H Lyman alpha",
			 1.0,
			 atm,
			 &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lya,
			 &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lya);

    grid.generate_S();
  }
  
  vector<double> brightness() {
    obs.reset_output();
    
    grid.brightness(obs);

    vector<double> brightness;
    brightness.resize(obs.size());

    for (int i=0;i<obs.size();i++)
      brightness[i] = obs.brightness[i][0];

    return brightness;
  }

  
};

#endif
