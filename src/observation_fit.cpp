//observation_fit.cpp -- routines to fit an atmosphere observation

#include "observation_fit.hpp"
using std::vector;
using std::string;

observation_fit::observation_fit()
  : emission_names({"H Lyman alpha"}),
    obs(emission_names),
    atm(/*nHexo = */5e5, /*nCO2exo = */1e9, temp), //these dummy values are overwritten later
    RT(emission_names)
{
  RT.grid.rmethod = RT.grid.rmethod_lognH;
  RT.grid.szamethod = RT.grid.szamethod_uniform_cos;
}

void observation_fit::add_observation(vector<Vector3> MSO_locations, vector<Vector3> MSO_directions) {
  obs.add_MSO_observation(MSO_locations,MSO_directions);
}

// void add_observation(Real* MSO_locations, Real* MSO_directions, int n_obs) {
//   obs.add_MSO_observation(MSO_locations,MSO_directions,n_obs);
// }

void observation_fit::set_g_factor(Real &g) {
  obs.emission_g_factors[0] = g;
}

void observation_fit::generate_source_function(Real nHexo, Real Texo) {
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
    
  //compute source function on the GPU if compiled with NVCC
#ifdef __CUDACC__
  RT.generate_S_gpu();
#else
  RT.generate_S();
#endif
}
  
vector<Real> observation_fit::brightness() {
  obs.reset_output();
    
  //compute brightness on the GPU if compiled with NVCC
#ifdef __CUDACC__
  RT.brightness_gpu(obs);
#else
  RT.brightness(obs);
#endif

  vector<Real> brightness;
  brightness.resize(obs.size());

  for (int i=0;i<obs.size();i++)
    brightness[i] = obs.los[i].brightness[0];

  return brightness;
}
