//observation_fit.cpp -- routines to fit an atmosphere observation

#include "observation_fit.hpp"
using std::vector;
using std::string;

observation_fit::observation_fit()
  : emission_names({"H Lyman alpha","H Lyman beta"}),
    obs(emission_names),
    atm(/*nHexo = */5e5, /*nCO2exo = */1e9, temp), //these dummy values are overwritten later
    atm_asym(/*nHexo = */5e5, /*nCO2exo = */1e9, temp),
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

void observation_fit::set_g_factor(vector<Real> &g) {
  for (int i_emission=0;i_emission<n_emissions;i_emission++)
    obs.emission_g_factors[i_emission] = g[i_emission];
}

void observation_fit::generate_source_function(const Real &nHexo, const Real &Texo) {

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
  RT.define_emission("H Lyman beta",
		     lyman_beta_branching_ratio,
		     atm,
		     &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lyb,
		     &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lyb);
  
  //compute source function on the GPU if compiled with NVCC
#ifdef __CUDACC__
  RT.generate_S_gpu();
#else
  RT.generate_S();
#endif
}

void observation_fit::generate_source_function_effv(const Real &nHexo, const Real &effv_exo) {
  Real Texo = Tconv.T_from_eff(effv_exo);
  
  generate_source_function(nHexo,Texo);
}

void observation_fit::generate_source_function_lc(const Real &nHexo, const Real &lc_exo) {
  Real Texo = Tconv.T_from_lc(lc_exo);
  
  generate_source_function(nHexo,Texo);
}


void observation_fit::generate_source_function_asym(const Real &nHexo, const Real &Texo,
						    const Real &asym) {

  temp = krasnopolsky_temperature(Texo);
  atm_asym = chamb_diff_1d_asymmetric(nHexo,CO2_exobase_density,temp);
  atm_asym.set_asymmetry(asym);

  
  RT.grid.setup_voxels(atm);
  RT.grid.setup_rays();

  //update the RT grid values
  RT.define_emission("H Lyman alpha",
		     1.0,
		     atm_asym,
		     &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lya,
		     &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lya);
  RT.define_emission("H Lyman beta",
		     lyman_beta_branching_ratio,
		     atm_asym,
		     &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lyb,
		     &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lyb);
  
  //compute source function on the GPU if compiled with NVCC
#ifdef __CUDACC__
  RT.generate_S_gpu();
#else
  RT.generate_S();
#endif
}



vector<vector<Real>> observation_fit::brightness() {
  obs.reset_output();
  
  //compute brightness on the GPU if compiled with NVCC
#ifdef __CUDACC__
  RT.brightness_gpu(obs);
#else
  RT.brightness(obs);
#endif
  
  vector<vector<Real>> brightness;
  brightness.resize(n_emissions);
  
  for (int i_emission=0;i_emission<n_emissions;i_emission++) {
    brightness[i_emission].resize(obs.size());
    
    for (int i=0;i<obs.size();i++)
      brightness[i_emission][i] = obs.los[i].brightness[i_emission];
  }
  
  return brightness;
}

void observation_fit::add_observed_brightness(const std::vector<Real> &brightness,
					      const std::vector<Real> &sigma,
					      const int emission/* = 0*/) {
  //add the observed brightness to obs (not obs_deriv)
  //so we can compute log-likelihoods
  assert(obs.size() == (int) brightness.size());
  for (int i=0;i<obs.size();i++) {
    obs.los_observed[i].brightness[emission] = brightness[i];
    obs.los_observed[i].sigma[emission]      = sigma[i];
  }
}

// vector<Real> observation_fit::likelihood_and_derivatives(const Real &nHexo, const Real &Texo) {
//   //computes likelihood and derivatives along each input dimension
//   const Real scale = 0.01;

//   //set up the atmospheres to simulate
//   Real center_parameters[n_parameters] = {nHexo, Texo};
//   Real parameters[n_parameters];

//   for (int j_parameter = 0; j_parameter < n_parameters; j_parameter++)
//     parameters[j_parameter] = center_parameter[j_parameter];
  
//   //atmosphere at center point
//   temp = krasnopolsky_temperature(parameters[1]);
//   atm = chamb_diff_1d(parameters[0],CO2_exobase_density,temp);

//   RT_deriv.grid.setup_voxels(atm);
//   RT_deriv.grid.setup_rays();

//   //set up emission 0, H Lyman alpha
//   i_emission = 0;
  
//   int i_simulate_center = i_emission * n_simulate_per_emission;
//   RT_deriv.define_emission(emission_names[i_emission],
// 			   1.0,
// 			   atm,
// 			   &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lya,
// 			   &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lya);
  
//   for (int i_parameter = 0; i_parameter < n_parameters; i_parameter++) { 
//     for (int i_derivative = 0; i_derivative < n_pts_per_derivative; i_derivative++) {
//       int i_simulate = i_simulate_center+i_parameter*n_pts_per_derivative+i_derivative+1;
      
//       for (int j_parameter = 0; j_parameter < n_parameters; j_parameter++)
// 	parameters[j_parameter] = center_parameter[j_parameter];
      
//       parameters[i_parameter] *= i_derivative == 0 ? 1.0-scale : 1.0+scale;
      
//       //atmosphere at this derivative point
//       temp = krasnopolsky_temperature(parameters[1]);
//       atm = chamb_diff_1d(parameters[0],CO2_exobase_density,temp);
      
//       RT_deriv.define_emission(simulate_names[i_simulate],
// 			       1.0,
// 			       atm,
// 			       &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lya,
// 			       &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lya);
      
//     }
//   }

//   //no more emissions to set up!
  

//   //simulate to retrieve model brightness
// #ifdef __CUDACC__
//   RT_deriv.generate_S_gpu();
// #else
//   RT_deriv.generate_S();
// #endif

//   obs_deriv.reset_output();
//   //compute brightness on the GPU if compiled with NVCC
// #ifdef __CUDACC__
//   RT_deriv.brightness_gpu(obs_deriv);
// #else
//   RT_deriv.brightness(obs_deriv);
// #endif

//   //compute log likelihood for each simulatation
// #ifdef __CUDACC__
//   logl_gpu();
// #else
//   logl();
// #endif

//   //compute derivative of log-likelihood for each parameter
//   for (int i_parameter = 0; i_parameter < n_parameters; i_parameter++) { 
//     for (int i_derivative = 0; i_derivative < n_pts_per_derivative; i_derivative++) {
//       int i_simulate = i_simulate_center+i_parameter*n_pts_per_derivative+i_derivative+1;

// 	for (int j_parameter = 0; j_parameter < n_parameters; j_parameter++)
// 	  parameters[j_parameter] = center_parameter[j_parameter];

// 	parameters[i_parameter] *= i_derivative == 0 ? 1.0-scale : 1.0+scale;
	
// 	//atmosphere at this derivative point
// 	temp = krasnopolsky_temperature(parameters[1]);
// 	atm = chamb_diff_1d(parameters[0],CO2_exobase_density,temp);
	
// 	RT_deriv.define_emission(simulate_names[i_simulate],
// 				 1.0,
// 				 atm,
// 				 &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lya,
// 				 &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lya);
	
//       }
//     }
//   }


//   //return log-likelihood and derivatives in the form pyMC3 prefers

// }
