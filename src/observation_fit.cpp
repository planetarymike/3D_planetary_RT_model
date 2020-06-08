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

  for (int i_emission = 0; i_emission < n_emissions; i_emission++) {
    int i_simulate_center = i_emission * n_simulate_per_emission;
    simulate_names[i_simulate_center] = emission_names[i_emission];

    for (int i_parameter = 0; i_parameter < n_parameters; i_parameter++) { 
      for (int i_derivative = 0; i_derivative < n_pts_per_derivative; i_derivative++) {
	int i_simulate = i_simulate_center+i_parameter*n_pts_per_derivative+i_derivative+1;

	string deriv = i_derivative == 0 ? "-" : "+";

	simulate_names[i_simulate] = emission_names[i_emission] + " " + std::to_string(i_parameter) + deriv;
      }
    }
  }

  obs_deriv.set_names(simulate_names);
  RT_deriv.set_names(simulate_names);

  RT.grid.rmethod = RT.grid.rmethod_lognH;
  RT.grid.szamethod = RT.grid.szamethod_uniform_cos;

  RT_deriv.grid.rmethod = RT_deriv.grid.rmethod_lognH;
  RT_deriv.grid.szamethod = RT_deriv.grid.szamethod_uniform_cos;
}

void observation_fit::add_observation(const vector<vector<Real>> &MSO_locations, const vector<vector<Real>> &MSO_directions) {
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

void observation_fit::generate_source_function_effv(Real nHexo, Real effv_exo) {
  Real Texo = Tconv.T_from_eff(effv_exo);
  
  generate_source_function(nHexo,Texo);
}

void observation_fit::generate_source_function_lc(Real nHexo, Real lc_exo) {
  Real Texo = Tconv.T_from_lc(lc_exo);
  
  generate_source_function(nHexo,Texo);
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

void add_observed_brightness(const std::vector<vector<Real>> &brightness,
			     const std::vector<vector<Real>> &sigma) {
  //add the observed brightness to obs (not obs_deriv)
  //so we can compute log-likelihoods

}

vector<Real> likelihood_and_derivatives(Real nHexo, Real Texo) {
  //computes likelihood and derivatives along each input dimension

  //set up the atmospheres to simulate

  //simulate to retrieve model brightness

  //compute log likelihood for each simulatation

  //compute derivative of log-likelihood for each parameter

  //return log-likelihood and derivatives in the form pyMC3 prefers

}
