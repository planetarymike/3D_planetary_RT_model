//observation_fit.cpp -- routines to fit an atmosphere observation

#include <string>
#include "constants.hpp"
#include "observation_fit.hpp"
#include "quemerais_IPH_model/iph_model_interface.hpp"
using std::vector;
using std::string;

observation_fit::observation_fit(const string iph_sfn_fnamee)
  : atm_tabular(),
    hydrogen_RT_pp(grid_pp, hydrogen_emissions_pp),
    hydrogen_RT(grid, hydrogen_emissions),
    hydrogen_obs(hydrogen_emissions),
    iph_sfn_fname(iph_sfn_fnamee),
    sim_iph(false),
    deuterium_RT_pp(grid_pp, deuterium_emissions_pp),
    deuterium_RT(grid, deuterium_emissions),
    deuterium_obs(deuterium_emissions),
    ly_multiplet_RT(grid, ly_multiplet_emissions),
    ly_multiplet_obs(ly_multiplet_emissions),
    ly_singlet_RT(grid, ly_singlet_emissions),
    ly_singlet_obs(ly_singlet_emissions),
    oxygen_RT(grid, oxygen_emissions),
    oxygen_obs(oxygen_emissions)
{
  //std::cout << "iph_sfn_fname = " << iph_sfn_fname << std::endl;
  
  CO2_exobase_density = default_CO2_exobase_density;

  hydrogen_RT_pp.grid.rmethod = grid_pp.rmethod_log_n_species;
  deuterium_RT_pp.grid.rmethod = grid_pp.rmethod_log_n_species;

  hydrogen_RT.grid.rmethod = grid.rmethod_altitude;//_tau_absorber;
  // fixed altitude grid gives smooth solutions for changes in input parameters,
  // other kinds of grids do not.
  hydrogen_RT.grid.szamethod = grid.szamethod_uniform_cos;
  hydrogen_RT.grid.raymethod_theta = grid.raymethod_theta_uniform;

  deuterium_RT.grid.rmethod = grid.rmethod_altitude;//_tau_absorber;
  // fixed altitude grid gives smooth solutions for changes in input parameters,
  // other kinds of grids do not.
  deuterium_RT.grid.szamethod = grid.szamethod_uniform_cos;
  deuterium_RT.grid.raymethod_theta = grid.raymethod_theta_uniform;

  ly_multiplet_RT.grid.rmethod = grid.rmethod_altitude;//_tau_absorber;
  // fixed altitude grid gives smooth solutions for changes in input parameters,
  // other kinds of grids do not.
  ly_multiplet_RT.grid.szamethod = grid.szamethod_uniform_cos;
  ly_multiplet_RT.grid.raymethod_theta = grid.raymethod_theta_uniform;

  ly_singlet_RT.grid.rmethod = grid.rmethod_altitude;//_tau_absorber;
  // fixed altitude grid gives smooth solutions for changes in input parameters,
  // other kinds of grids do not.
  ly_singlet_RT.grid.szamethod = grid.szamethod_uniform_cos;
  ly_singlet_RT.grid.raymethod_theta = grid.raymethod_theta_uniform;

  oxygen_RT.grid.rmethod = grid.rmethod_log_n_species;
  oxygen_RT.grid.szamethod = grid.szamethod_uniform_cos;
  oxygen_RT.grid.raymethod_theta = grid.raymethod_theta_uniform;

  tweak_H_density = false;
  tweak_H_temp = false;
}

void observation_fit::add_observation(const vector<vector<Real>> &MSO_locations, const vector<vector<Real>> &MSO_directions) {
  hydrogen_obs.add_MSO_observation(MSO_locations,MSO_directions);
  deuterium_obs.add_MSO_observation(MSO_locations,MSO_directions);
  ly_multiplet_obs.add_MSO_observation(MSO_locations,MSO_directions);
  ly_singlet_obs.add_MSO_observation(MSO_locations,MSO_directions);
  oxygen_obs.add_MSO_observation(MSO_locations,MSO_directions);
}

void observation_fit::set_g_factor(vector<Real> &g) {
  for (int i_emission=0;i_emission<n_hydrogen_emissions;i_emission++) {
    hydrogen_emissions[i_emission]->set_emission_g_factor(g[i_emission]);
    deuterium_emissions[i_emission]->set_emission_g_factor(g[i_emission]);
    hydrogen_emissions_pp[i_emission]->set_emission_g_factor(g[i_emission]);
  }


  std::cout << "Ly alpha solar brightness = " << g[0]/lyman_alpha_cross_section_total << std::endl;
  std::cout << "Ly beta solar brightness = " << g[1]/lyman_beta_cross_section_total << std::endl;
  
  ly_multiplet.set_solar_brightness(g[0]/lyman_alpha_cross_section_total,
				    g[1]/lyman_beta_cross_section_total);
  ly_singlet.set_solar_brightness(g[0]/lyman_alpha_cross_section_total,
				  g[1]/lyman_beta_cross_section_total);
}

void observation_fit::simulate_iph(const bool sim_iphh) {
  sim_iph = sim_iphh;
}

void observation_fit::add_observation_ra_dec(const std::vector<Real> &mars_ecliptic_coords,
					     const std::vector<Real> &RAA,
					     const std::vector<Real> &Decc) {
  simulate_iph(true);
  hydrogen_obs.add_observation_ra_dec(mars_ecliptic_coords,
				      RAA,
				      Decc);
  get_unextincted_iph();
}

void observation_fit::get_unextincted_iph() {
  //simulate the IPH brightness using Quemerais' IPH code
  vector<Real> iph_brightness_lya = quemerais_iph_model(iph_sfn_fname,
							lyman_alpha.get_emission_g_factor(),
							hydrogen_obs.mars_ecliptic_pos,
							hydrogen_obs.ra, hydrogen_obs.dec);
  
  for (int i_obs=0; i_obs < hydrogen_obs.size(); i_obs++) {
    hydrogen_obs.iph_brightness_unextincted[i_obs][0] = iph_brightness_lya[i_obs];
    if (n_hydrogen_emissions==2)
      hydrogen_obs.iph_brightness_unextincted[i_obs][1] = (lyman_beta.get_emission_g_factor()
							   / lyman_alpha.get_emission_g_factor()
							   * iph_brightness_lya[i_obs]); 
  }
}


void observation_fit::generate_source_function(const Real &nHexo, const Real &Texo,
					       const string atmosphere_fname/* = ""*/,
					       const string sourcefn_fname/* = ""*/,
					       const bool plane_parallel/* = false*/,
					       const bool deuterium/* = false */)
{
  //std::cout << "nHexo = " << nHexo << "; Texo = " << Texo << ".\n";

  temp = krasnopolsky_temperature(Texo);
  chamb_diff_1d atm(nHexo,
		    CO2_exobase_density,
		    &temp,
		    deuterium ? &D_thermosphere : &H_thermosphere);
  atm.copy_H_options(H_cross_section_options);
  
  if (atmosphere_fname != "")
    atm.save(atmosphere_fname);
  
  if (plane_parallel) {
    if (!deuterium) {
      generate_source_function_plane_parallel(atm, Texo,
					      hydrogen_RT_pp,
					      lyman_alpha_pp,
					      lyman_beta_pp,
					      sourcefn_fname);
    } else {
      generate_source_function_plane_parallel(atm, Texo,
					      deuterium_RT_pp,
					      D_lyman_alpha_pp,
					      D_lyman_beta_pp,
					      sourcefn_fname);
    }
  } else {
    if (!deuterium) {
      generate_source_function_sph_azi_sym(atm, Texo,
					   hydrogen_RT,
					   lyman_alpha,
					   lyman_beta,
					   sourcefn_fname);
    } else {
      generate_source_function_sph_azi_sym(atm, Texo,
					   deuterium_RT,
					   D_lyman_alpha,
					   D_lyman_beta,
					   sourcefn_fname);
    }
  }
}

void observation_fit::generate_source_function_effv(const Real &nHexo, const Real &effv_exo,
						    const string atmosphere_fname/* = ""*/,
						    const string sourcefn_fname/* = ""*/,
						    const bool plane_parallel/* = false*/,
						    const bool deuterium/* =false */)
{
  Real Texo = Tconv.T_from_eff(effv_exo);
  
  generate_source_function(nHexo,Texo,
			   atmosphere_fname,
			   sourcefn_fname,
			   plane_parallel,
			   deuterium);
}

void observation_fit::generate_source_function_lc(const Real &nHexo, const Real &lc_exo,
						  const string atmosphere_fname/* = ""*/,
						  const string sourcefn_fname/* = ""*/,
						  const bool plane_parallel/* = false*/,
						  const bool deuterium/* =false */)
{
  Real Texo = Tconv.T_from_lc(lc_exo);
  
  generate_source_function(nHexo,Texo,
			   atmosphere_fname,
			   sourcefn_fname,
			   plane_parallel,
			   deuterium);
}


void observation_fit::generate_source_function_variable_thermosphere(const Real &nHexo,
								     const Real &Texo,
								     const Real &nCO2rminn, //a good number is 2.6e13 (Chaufray2008)
								     const Real rexoo,
								     const Real rminn,
								     const Real rmaxx,
								     const Real rmindiffusionn,
								     //extra args for krasnopolsky_temp					  
								     const Real T_tropo,
								     const Real r_tropo,
								     const Real shape_parameter,
								     const string atmosphere_fname/* = ""*/,
								     const string sourcefn_fname/* = ""*/,
								     const bool plane_parallel/* = false*/,
								     const bool deuterium/* =false */)
{
  // std::cout << "nHexo = " << nHexo << "; Texo = " << Texo << ".\n";
  // std::cout << "nCO2rmin = " << nCO2rminn << ".\n";
  // std::cout << "T_tropo = " << T_tropo << "; z_tropo = " << (r_tropo-rMars)/1e5 << "; shape_parameter = " << shape_parameter << ".\n";

  temp = krasnopolsky_temperature(Texo, T_tropo, r_tropo, shape_parameter, false/*shape parameter is in absolute units of km*/);
  chamb_diff_1d atm(rminn,
		    rexoo,
		    rmaxx,
		    rmindiffusionn,
		    nHexo,
		    nCO2rminn,
		    &temp,
		    deuterium ? &D_thermosphere : &H_thermosphere,
		    thermosphere_exosphere::method_rmax_nCO2rmin);
  atm.copy_H_options(H_cross_section_options);
  
  if (atmosphere_fname !="")
    atm.save(atmosphere_fname);

  if (plane_parallel) {
    if (!deuterium) {
      generate_source_function_plane_parallel(atm, Texo,
					      hydrogen_RT_pp,
					      lyman_alpha_pp,
					      lyman_beta_pp,
					      sourcefn_fname);
    } else {
      generate_source_function_plane_parallel(atm, Texo,
					      deuterium_RT_pp,
					      D_lyman_alpha_pp,
					      D_lyman_beta_pp,
					      sourcefn_fname);
    }
  } else {
    if (!deuterium) {
      generate_source_function_sph_azi_sym(atm, Texo,
					   hydrogen_RT,
					   lyman_alpha,
					   lyman_beta,
					   sourcefn_fname);
    } else {
      generate_source_function_sph_azi_sym(atm, Texo,
					   deuterium_RT,
					   D_lyman_alpha,
					   D_lyman_beta,
					   sourcefn_fname);
    }
  }
}

void observation_fit::generate_source_function_nH_asym(const Real &nHexo, const Real &Texo,
						       const Real &asym,
						       const string sourcefn_fname/* = ""*/,
						       const bool deuterium/* =false */) {
  
  temp = krasnopolsky_temperature(Texo);
  chamb_diff_1d_asymmetric atm_asym(nHexo,CO2_exobase_density,&temp,&H_thermosphere);
  atm_asym.copy_H_options(H_cross_section_options);
  atm_asym.set_asymmetry(asym);

  if (!deuterium) {
    generate_source_function_sph_azi_sym(atm_asym, Texo,
					 hydrogen_RT,
					 lyman_alpha,
					 lyman_beta,
					 sourcefn_fname);
  } else {
    generate_source_function_sph_azi_sym(atm_asym, Texo,
					 deuterium_RT,
					 D_lyman_alpha,
					 D_lyman_beta,
					 sourcefn_fname);
  }
}

void observation_fit::generate_source_function_temp_asym(const Real &nHavg,
							 const Real &Tnoon, const Real &Tmidnight,
							 const string sourcefn_fname/* = ""*/,
							 const bool deuterium/* =false */) {
  generate_source_function_temp_asym(nHavg,
				     Tnoon, Tmidnight,
				     /*      nCO2rmin = */2.6e13,			       
				     /*          rexo = */rexo_typical,
				     /*          rmin = */rMars + 80e5,
				     /*         rmaxx = */rMars + 50000e5,
				     /* rmindiffusion = */rMars + 80e5,
				     //extra args for krasnopolsky_temp					  
				     /*         T_tropo = */125.0,
				     /*         r_tropo = */rMars + 90e5,
				     /* shape_parameter = */11.4,
				     /*          Tpower = */2.5,
				     sourcefn_fname,
				     deuterium);
}

void observation_fit::generate_source_function_temp_asym(const Real &nHavg,
							 const Real &Tnoon, const Real &Tmidnight,
							 const Real nCO2rminn, //a good number is 2.6e13 (Chaufray2008)
							 const Real rexoo,
							 const Real rminn,
							 const Real rmaxx,
							 const Real rmindiffusionn,
							 //extra args for krasnopolsky_temp					  
							 const Real T_tropo,
							 const Real r_tropo,
							 const Real shape_parameter,
							 //power for temperature in the expression n*T^p = const.
							 const Real Tpowerr,
							 const string sourcefn_fname/* = ""*/,
							 const bool deuterium/* =false */) {
  // std::cout << "nHavg = " << nHavg << "; Tnoon = " << Tnoon << "; Tmidnight = " << Tmidnight <<".\n";
  // std::cout << "nCO2rmin = " << nCO2rminn << ".\n";
  // std::cout << "T_tropo = " << T_tropo << "; z_tropo = " << (r_tropo-rMars)/1e5 << "; shape_parameter = " << shape_parameter << ".\n";
  // std::cout << "Tpower = " << Tpowerr << ".\n";

  chamb_diff_temp_asymmetric atm_asym(&H_thermosphere,
				      nHavg,
				      Tnoon, Tmidnight,
				      nCO2rminn,
				      rexoo,
				      rminn,
				      rmaxx,
				      rmindiffusionn,
				      //extra args for krasnopolsky_temp					  
				      T_tropo,
				      r_tropo,
				      shape_parameter,
				      Tpowerr);
  atm_asym.copy_H_options(H_cross_section_options);
  
  if (!deuterium) {
    generate_source_function_sph_azi_sym(atm_asym, Tnoon,
					 hydrogen_RT,
					 lyman_alpha,
					 lyman_beta,
					 sourcefn_fname);
  } else {
    generate_source_function_sph_azi_sym(atm_asym, Tnoon,
					 deuterium_RT,
					 D_lyman_alpha,
					 D_lyman_beta,
					 sourcefn_fname);
  }
}

void observation_fit::generate_source_function_tabular_atmosphere(const Real rmin, const Real rexo, const Real rmax,
								  const std::vector<doubReal> &alt_nH, const std::vector<doubReal> &log_nH,
								  const std::vector<doubReal> &alt_nCO2, const std::vector<doubReal> &log_nCO2,
								  const std::vector<doubReal> &alt_temp, const std::vector<doubReal> &temp,
								  const bool compute_exosphere/* = false*/,
								  const bool plane_parallel/*= false*/,
								  const bool deuterium/* =false */,
								  const string sourcefn_fname/* = ""*/) {

  tabular_1d new_atm_tabular(rmin,rexo,rmax,compute_exosphere);
  new_atm_tabular.copy_H_options(H_cross_section_options);
  atm_tabular = new_atm_tabular;

  atm_tabular.load_log_species_density(alt_nH, log_nH);
  atm_tabular.load_log_absorber_density(alt_nCO2, log_nCO2);
  atm_tabular.load_temperature(alt_temp, temp);

  Real Texo = atm_tabular.Temp(rexo);

  if (plane_parallel) {
    if (!deuterium) {
      generate_source_function_plane_parallel(atm_tabular, Texo,
					      hydrogen_RT_pp,
					      lyman_alpha_pp,
					      lyman_beta_pp,
					      sourcefn_fname);
    } else {
      generate_source_function_plane_parallel(atm_tabular, Texo,
					      deuterium_RT_pp,
					      D_lyman_alpha_pp,
					      D_lyman_beta_pp,
					      sourcefn_fname);
    }
  } else {
    if (!deuterium) {
      generate_source_function_sph_azi_sym(atm_tabular, Texo,
					   hydrogen_RT,
					   lyman_alpha,
					   lyman_beta,
					   sourcefn_fname);
    } else {
      generate_source_function_sph_azi_sym(atm_tabular, Texo,
					   deuterium_RT,
					   D_lyman_alpha,
					   D_lyman_beta,
					   sourcefn_fname);
    }
  }
}

void observation_fit::set_use_CO2_absorption(const bool use_CO2_absorption/* = true*/) {
  H_cross_section_options.no_CO2_absorption = !use_CO2_absorption;
  atm_tabular.no_CO2_absorption = !use_CO2_absorption;
  use_CO2_absorption ? ly_multiplet.set_CO2_absorption_on() : ly_multiplet.set_CO2_absorption_off();
  use_CO2_absorption ? ly_singlet.set_CO2_absorption_on() : ly_singlet.set_CO2_absorption_off();
}
void observation_fit::set_use_temp_dependent_sH(const bool use_temp_dependent_sH/* = true*/, const Real constant_temp_sH/* = -1*/) {
  H_cross_section_options.temp_dependent_sH = use_temp_dependent_sH;
  atm_tabular.temp_dependent_sH = use_temp_dependent_sH;
  
  H_cross_section_options.constant_temp_sH = constant_temp_sH;
  atm_tabular.constant_temp_sH = constant_temp_sH;

  if (use_temp_dependent_sH) {
    ly_multiplet.set_atmosphere_temp_RT();
    ly_singlet.set_atmosphere_temp_RT();
  } else {
    ly_multiplet.set_constant_temp_RT(constant_temp_sH);
    ly_singlet.set_constant_temp_RT(constant_temp_sH);
  }
}

void observation_fit::set_sza_method_uniform() {
  hydrogen_RT.grid.szamethod = hydrogen_RT.grid.szamethod_uniform;
  ly_multiplet_RT.grid.szamethod = ly_multiplet_RT.grid.szamethod_uniform;
  ly_singlet_RT.grid.szamethod = ly_singlet_RT.grid.szamethod_uniform;
}
void observation_fit::set_sza_method_uniform_cos() {
  hydrogen_RT.grid.szamethod = hydrogen_RT.grid.szamethod_uniform_cos;
  ly_multiplet_RT.grid.szamethod = ly_multiplet_RT.grid.szamethod_uniform_cos;
  ly_singlet_RT.grid.szamethod = ly_singlet_RT.grid.szamethod_uniform_cos;
}

void observation_fit::reset_H_lya_xsec_coef(const Real xsec_coef/* = lyman_alpha_line_center_cross_secion_coef*/) {
  H_cross_section_options.H_lya_xsec_coef = xsec_coef;
  atm_tabular.H_lya_xsec_coef = xsec_coef;
}
void observation_fit::reset_H_lyb_xsec_coef(const Real xsec_coef/* = lyman_beta_line_center_cross_section_coef*/) {
  H_cross_section_options.H_lyb_xsec_coef = xsec_coef;
  atm_tabular.H_lyb_xsec_coef = xsec_coef;
}
void observation_fit::reset_CO2_lya_xsec(const Real xsec/* = CO2_lyman_alpha_absorption_cross_section*/) {
  H_cross_section_options.CO2_lya_xsec = xsec;
  atm_tabular.CO2_lya_xsec = xsec;
}
void observation_fit::reset_CO2_lyb_xsec(const Real xsec/* = CO2_lyman_beta_absorption_cross_section*/) {
  H_cross_section_options.CO2_lyb_xsec = xsec;
  atm_tabular.CO2_lyb_xsec = xsec;
}


Real observation_fit::get_CO2_exobase_density() {
  return CO2_exobase_density;
}
void observation_fit::reset_CO2_exobase_density() {
  CO2_exobase_density = default_CO2_exobase_density;
}
void observation_fit::set_CO2_exobase_density(const double nCO2) {
  CO2_exobase_density = nCO2;
}

void observation_fit::set_H_density_tweak(const bool tweak_H_densityy) {
  tweak_H_density = tweak_H_densityy;
}
void observation_fit::set_H_density_tweak_values(const vector<int> voxels_to_tweak, const Real tweak_factor) {
  tweak_H_density_voxel_numbers = voxels_to_tweak;
  tweak_H_density_factor = tweak_factor;
}
void observation_fit::set_H_temp_tweak(const bool tweak_H_tempp) {
  tweak_H_temp = tweak_H_tempp;
}
void observation_fit::set_H_temp_tweak_values(const vector<int> voxels_to_tweak, const Real tweak_factor) {
  tweak_H_temp_voxel_numbers = voxels_to_tweak;
  tweak_H_temp_factor = tweak_factor;
}



vector<vector<Real>> observation_fit::brightness() {
  //compute brightness on the GPU if compiled with NVCC
#ifdef __CUDACC__
  hydrogen_RT.brightness_gpu(hydrogen_obs);
#else
  hydrogen_RT.brightness(hydrogen_obs);
#endif

  if (sim_iph)
    hydrogen_obs.update_iph_extinction();
  
  vector<vector<Real>> brightness;
  brightness.resize(n_hydrogen_emissions);
  
  for (int i_emission=0;i_emission<n_hydrogen_emissions;i_emission++) {
    brightness[i_emission].resize(hydrogen_obs.size());
    
    for (int i=0;i<hydrogen_obs.size();i++) {
      brightness[i_emission][i] = hydrogen_obs.los[i_emission][i].brightness;
      if (sim_iph)
	brightness[i_emission][i] += hydrogen_obs.iph_brightness_observed[i][i_emission];
    }
  }
  
  return brightness;
}

vector<vector<Real>> observation_fit::species_col_dens() {
  vector<vector<Real>> species_col_dens;
  species_col_dens.resize(n_hydrogen_emissions);
  
  for (int i_emission=0;i_emission<n_hydrogen_emissions;i_emission++) {
    species_col_dens[i_emission].resize(hydrogen_obs.size());
    
    for (int i=0;i<hydrogen_obs.size();i++)
      species_col_dens[i_emission][i] = hydrogen_obs.los[i_emission][i].species_col_dens;
  }
  
  return species_col_dens;
}

vector<vector<Real>> observation_fit::tau_species_final() {
  vector<vector<Real>> tau_species;
  tau_species.resize(n_hydrogen_emissions);
  
  for (int i_emission=0;i_emission<n_hydrogen_emissions;i_emission++) {
    tau_species[i_emission].resize(hydrogen_obs.size());
    
    for (int i=0;i<hydrogen_obs.size();i++)
      tau_species[i_emission][i] = hydrogen_obs.los[i_emission][i].tau_species_final;
  }
  
  return tau_species;
}


vector<vector<Real>> observation_fit::tau_absorber_final() {
  vector<vector<Real>> tau_absorber;
  tau_absorber.resize(n_hydrogen_emissions);
  
  for (int i_emission=0;i_emission<n_hydrogen_emissions;i_emission++) {
    tau_absorber[i_emission].resize(hydrogen_obs.size());
    
    for (int i=0;i<hydrogen_obs.size();i++)
      tau_absorber[i_emission][i] = hydrogen_obs.los[i_emission][i].tau_absorber_final;
  }
  
  return tau_absorber;
}


vector<vector<Real>> observation_fit::iph_brightness_observed() {
  vector<vector<Real>> iph_b;
  iph_b.resize(n_hydrogen_emissions);

  hydrogen_obs.update_iph_extinction();
  
  for (int i_emission=0;i_emission<n_hydrogen_emissions;i_emission++) {
    iph_b[i_emission].resize(hydrogen_obs.size());
    
    for (int i=0;i<hydrogen_obs.size();i++)
      iph_b[i_emission][i] =  hydrogen_obs.iph_brightness_observed[i][i_emission];
  }
  
  return iph_b;
}

vector<vector<Real>> observation_fit::iph_brightness_unextincted() {
  vector<vector<Real>> iph_b;
  iph_b.resize(n_hydrogen_emissions);

  for (int i_emission=0;i_emission<n_hydrogen_emissions;i_emission++) {
    iph_b[i_emission].resize(hydrogen_obs.size());
    
    for (int i=0;i<hydrogen_obs.size();i++)
      iph_b[i_emission][i] =  hydrogen_obs.iph_brightness_unextincted[i][i_emission];
  }
  
  return iph_b;
}

vector<vector<Real>> observation_fit::D_brightness() {
  //compute brightness on the GPU if compiled with NVCC
#ifdef __CUDACC__
  deuterium_RT.brightness_gpu(deuterium_obs);
#else
  deuterium_RT.brightness(deuterium_obs);
#endif

  if (sim_iph)
    deuterium_obs.update_iph_extinction();
  
  vector<vector<Real>> brightness;
  brightness.resize(n_hydrogen_emissions);
  
  for (int i_emission=0;i_emission<n_hydrogen_emissions;i_emission++) {
    brightness[i_emission].resize(deuterium_obs.size());
    
    for (int i=0;i<deuterium_obs.size();i++) {
      brightness[i_emission][i] = deuterium_obs.los[i_emission][i].brightness;
      if (sim_iph)
	brightness[i_emission][i] += deuterium_obs.iph_brightness_observed[i][i_emission];
    }
  }
  
  return brightness;
}

vector<vector<Real>> observation_fit::D_col_dens() {
  vector<vector<Real>> species_col_dens;
  species_col_dens.resize(n_hydrogen_emissions);
  
  for (int i_emission=0;i_emission<n_hydrogen_emissions;i_emission++) {
    species_col_dens[i_emission].resize(deuterium_obs.size());
    
    for (int i=0;i<deuterium_obs.size();i++)
      species_col_dens[i_emission][i] = deuterium_obs.los[i_emission][i].species_col_dens;
  }
  
  return species_col_dens;
}

vector<vector<Real>> observation_fit::tau_D_final() {
  vector<vector<Real>> tau_species;
  tau_species.resize(n_hydrogen_emissions);
  
  for (int i_emission=0;i_emission<n_hydrogen_emissions;i_emission++) {
    tau_species[i_emission].resize(deuterium_obs.size());
    
    for (int i=0;i<deuterium_obs.size();i++)
      tau_species[i_emission][i] = deuterium_obs.los[i_emission][i].tau_species_final;
  }
  
  return tau_species;
}

void observation_fit::save_influence_matrix(const string fname) {
   hydrogen_RT.save_influence(fname);
}

void observation_fit::save_influence_matrix_O_1026(const string fname) {
   oxygen_RT.save_influence(fname);
}


void observation_fit::O_1026_generate_source_function(const Real &nOexo,
						      const Real &Texo,
						      const Real &solar_brightness_lyman_beta, // 1.69e-3 is a good number for solar minimum
						      const string atmosphere_fname/* = ""*/,
						      const string sourcefn_fname/* = ""*/)
{
  // std::cout << "Simulating O 102.6 nm brightness.\n";
  // std::cout << "nOexo = " << nOexo << "; Texo = " << Texo << ".\n";

  temp = krasnopolsky_temperature(Texo);

  Real rmin = rMars+100e5; // 100 km
  int method = thermosphere_exosphere::method_nspmin_nCO2exo;
  
  chamb_diff_1d atm(/* rmin = */ rmin,
		    /* rexo = */ rMars+200e5,
		    /* rmaxx_or_nspmin = */ 10, // cm3
		    /* rmindifussion = */ rmin, 
		    nOexo,
		    CO2_exobase_density,
		    &temp,
		    &O_thermosphere,
		    method);

  if (atmosphere_fname !="")
    atm.save(atmosphere_fname);
  
  atm.spherical = true;

  oxygen_RT.grid.setup_voxels(atm);
  oxygen_RT.grid.setup_rays();

  //update the emission density values
  oxygen_1026.define("O_1026",
		     atm,
		     &chamb_diff_1d::n_species_voxel_avg,
		     &chamb_diff_1d::Temp_voxel_avg,
		     &chamb_diff_1d::n_absorber_voxel_avg,
		     oxygen_RT.grid.voxels);
  oxygen_1026.set_solar_brightness(solar_brightness_lyman_beta); /* ph/cm2/s/Hz, 
								    solar line center brightness at Lyman beta */
  
  //compute source function on the GPU if compiled with NVCC
#ifdef __CUDACC__
  oxygen_RT.generate_S_gpu();
#else
  oxygen_RT.generate_S();
#endif
  
  if (sourcefn_fname!="")
    oxygen_RT.save_S(sourcefn_fname);

}

vector<vector<Real>> observation_fit::O_1026_brightness() {
  //compute brightness on the GPU if compiled with NVCC
#ifdef __CUDACC__
  oxygen_RT.brightness_gpu(oxygen_obs);
#else
  oxygen_RT.brightness(oxygen_obs);
#endif
  
  vector<vector<Real>> brightness;
  brightness.resize(oxygen_1026.n_lines);
  
  for (int i_line=0;i_line<oxygen_1026.n_lines;i_line++) {
    brightness[i_line].resize(oxygen_obs.size());
    
    for (int i=0;i<oxygen_obs.size();i++)
      brightness[i_line][i] = oxygen_obs.los[0][i].brightness[i_line];
  }
  
  return brightness;
}


void observation_fit::lyman_multiplet_generate_source_function(const Real &nHexo,
							       const Real &Texo,
							       // const Real &solar_brightness_lyman_alpha,
							       // const Real &solar_brightness_lyman_beta,
							       const string atmosphere_fname/* = ""*/,
							       const string sourcefn_fname/* = ""*/)
{
  temp = krasnopolsky_temperature(Texo);
  chamb_diff_1d atm(nHexo,
		    CO2_exobase_density,
		    &temp,
		    &H_thermosphere);
  atm.copy_H_options(H_cross_section_options);

  if (atmosphere_fname !="")
    atm.save(atmosphere_fname);
  
  atm.spherical = true;

  ly_multiplet_RT.grid.setup_voxels(atm);
  ly_multiplet_RT.grid.setup_rays();

  //update the emission density values
  ly_multiplet.define("Multiplet Lyman alpha and beta",
		      atm,
		      &chamb_diff_1d::n_species_voxel_avg,
		      &chamb_diff_1d::Temp_voxel_avg,
		      &chamb_diff_1d::n_absorber_voxel_avg,
		      ly_multiplet_RT.grid.voxels);
  // ly_multiplet.set_solar_brightness(solar_brightness_lyman_alpha,
  // 				    solar_brightness_lyman_beta);
  
  //compute source function on the GPU if compiled with NVCC
#ifdef __CUDACC__
  ly_multiplet_RT.generate_S_gpu();
#else
  ly_multiplet_RT.generate_S();
#endif
  
  if (sourcefn_fname!="")
    ly_multiplet_RT.save_S(sourcefn_fname);
}

vector<vector<Real>> observation_fit::lyman_multiplet_brightness() {
  //compute brightness on the GPU if compiled with NVCC
#ifdef __CUDACC__
  ly_multiplet_RT.brightness_gpu(ly_multiplet_obs);
#else
  ly_multiplet_RT.brightness(ly_multiplet_obs);
#endif
  
  vector<vector<Real>> brightness;
  brightness.resize(2);

  brightness[0].resize(ly_multiplet_obs.size());
  brightness[1].resize(ly_multiplet_obs.size());
  
  for (int i=0;i<ly_multiplet_obs.size();i++) {
    // lyman alpha
    brightness[0][i] = (ly_multiplet_obs.los[0][i].brightness[0]
			+ly_multiplet_obs.los[0][i].brightness[1]);

    // lyman beta
    brightness[1][i] = (ly_multiplet_obs.los[0][i].brightness[2]
			+ly_multiplet_obs.los[0][i].brightness[3]);
  }
  
  return brightness;
}

void observation_fit::lyman_singlet_generate_source_function(const Real &nHexo,
							       const Real &Texo,
							       // const Real &solar_brightness_lyman_alpha,
							       // const Real &solar_brightness_lyman_beta,
							       const string atmosphere_fname/* = ""*/,
							       const string sourcefn_fname/* = ""*/)
{
  temp = krasnopolsky_temperature(Texo);
  chamb_diff_1d atm(nHexo,
		    CO2_exobase_density,
		    &temp,
		    &H_thermosphere);
  atm.copy_H_options(H_cross_section_options);

  if (atmosphere_fname !="")
    atm.save(atmosphere_fname);
  
  atm.spherical = true;

  ly_singlet_RT.grid.setup_voxels(atm);
  ly_singlet_RT.grid.setup_rays();

  //update the emission density values
  ly_singlet.define("Singlet Lyman alpha and beta",
		      atm,
		      &chamb_diff_1d::n_species_voxel_avg,
		      &chamb_diff_1d::Temp_voxel_avg,
		      &chamb_diff_1d::n_absorber_voxel_avg,
		      ly_singlet_RT.grid.voxels);
  // ly_singlet.set_solar_brightness(solar_brightness_lyman_alpha,
  // 				    solar_brightness_lyman_beta);
  
  //compute source function on the GPU if compiled with NVCC
#ifdef __CUDACC__
  ly_singlet_RT.generate_S_gpu();
#else
  ly_singlet_RT.generate_S();
#endif
  
  if (sourcefn_fname!="")
    ly_singlet_RT.save_S(sourcefn_fname);
}

vector<vector<Real>> observation_fit::lyman_singlet_brightness() {
  //compute brightness on the GPU if compiled with NVCC
#ifdef __CUDACC__
  ly_singlet_RT.brightness_gpu(ly_singlet_obs);
#else
  ly_singlet_RT.brightness(ly_singlet_obs);
#endif
  
  vector<vector<Real>> brightness;
  brightness.resize(2);

  brightness[0].resize(ly_singlet_obs.size());
  brightness[1].resize(ly_singlet_obs.size());
  
  for (int i=0;i<ly_singlet_obs.size();i++) {
    // lyman alpha
    brightness[0][i] = ly_singlet_obs.los[0][i].brightness[0];

    // lyman beta
    brightness[1][i] = ly_singlet_obs.los[0][i].brightness[1];
  }
  
  return brightness;
}
