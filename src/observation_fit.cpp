//observation_fit.cpp -- routines to fit an atmosphere observation

#include "constants.hpp"
#include "observation_fit.hpp"
#include "quemerais_IPH_model/iph_model_interface.hpp"
using std::vector;
using std::string;

observation_fit::observation_fit()
  : atm_tabular(),
    hydrogen_RT_pp(grid_pp, hydrogen_emissions_pp),
    hydrogen_RT(grid, hydrogen_emissions),
    hydrogen_obs(hydrogen_emissions),
    sim_iph(false),
    oxygen_RT(grid, oxygen_emissions),
    oxygen_obs(oxygen_emissions)
{
  hydrogen_RT_pp.grid.rmethod = grid_pp.rmethod_log_n_species;

  hydrogen_RT.grid.rmethod = grid.rmethod_altitude;//_tau_absorber;
  // fixed altitude grid gives smooth solutions for changes in input parameters,
  // other kinds of grids do not.
  hydrogen_RT.grid.szamethod = grid.szamethod_uniform_cos;
  hydrogen_RT.grid.raymethod_theta = grid.raymethod_theta_uniform;

  oxygen_RT.grid.rmethod = grid.rmethod_log_n_species;
  oxygen_RT.grid.szamethod = grid.szamethod_uniform_cos;
  oxygen_RT.grid.raymethod_theta = grid.raymethod_theta_uniform;
}

void observation_fit::add_observation(const vector<vector<Real>> &MSO_locations, const vector<vector<Real>> &MSO_directions) {
  hydrogen_obs.add_MSO_observation(MSO_locations,MSO_directions);
  oxygen_obs.add_MSO_observation(MSO_locations,MSO_directions);
}

void observation_fit::set_g_factor(vector<Real> &g) {
  for (int i_emission=0;i_emission<n_hydrogen_emissions;i_emission++) {
    hydrogen_emissions[i_emission]->set_emission_g_factor(g[i_emission]);
    hydrogen_emissions_pp[i_emission]->set_emission_g_factor(g[i_emission]);
  }
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
  vector<Real> iph_brightness_lya = quemerais_iph_model(lyman_alpha.get_emission_g_factor(),
							hydrogen_obs.mars_ecliptic_pos,
							hydrogen_obs.ra, hydrogen_obs.dec);
  
  for (int i_obs=0; i_obs < hydrogen_obs.size(); i_obs++) {
    hydrogen_obs.iph_brightness_unextincted[i_obs][0] = iph_brightness_lya[i_obs];
    if (n_hydrogen_emissions==2)
      hydrogen_obs.iph_brightness_unextincted[i_obs][1] = 0; //could maybe estimate using lyman alpha?
  }
}


void observation_fit::generate_source_function(const Real &nHexo, const Real &Texo,
					       const string atmosphere_fname/* = ""*/,
					       const string sourcefn_fname/* = ""*/,
					       bool plane_parallel/* = false*/)
{
  std::cout << "nHexo = " << nHexo << "; Texo = " << Texo << ".\n";

  temp = krasnopolsky_temperature(Texo);
  chamb_diff_1d atm(nHexo,CO2_exobase_density,&temp,&H_thermosphere);
  atm.copy_H_options(H_cross_section_options);
  
  if (atmosphere_fname != "")
    atm.save(atmosphere_fname);

  if (plane_parallel) {
    generate_source_function_plane_parallel(atm,Texo,
					    sourcefn_fname);
  } else {
    generate_source_function_sph_azi_sym(atm,Texo,
					 sourcefn_fname);
  }
}

void observation_fit::generate_source_function_effv(const Real &nHexo, const Real &effv_exo,
						    const string atmosphere_fname/* = ""*/,
						    const string sourcefn_fname/* = ""*/,
						    bool plane_parallel/* = false*/)
{
  Real Texo = Tconv.T_from_eff(effv_exo);
  
  generate_source_function(nHexo,Texo,
			   atmosphere_fname,
			   sourcefn_fname,
			   plane_parallel);
}

void observation_fit::generate_source_function_lc(const Real &nHexo, const Real &lc_exo,
						  const string atmosphere_fname/* = ""*/,
						  const string sourcefn_fname/* = ""*/,
						  bool plane_parallel/* = false*/)
{
  Real Texo = Tconv.T_from_lc(lc_exo);
  
  generate_source_function(nHexo,Texo,
			   atmosphere_fname,
			   sourcefn_fname,
			   plane_parallel);
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
								     bool plane_parallel/* = false*/)
{
  std::cout << "nHexo = " << nHexo << "; Texo = " << Texo << ".\n";
  std::cout << "nCO2rmin = " << nCO2rminn << ".\n";
  std::cout << "T_tropo = " << T_tropo << "; z_tropo = " << (r_tropo-rMars)/1e5 << "; shape_parameter = " << shape_parameter << ".\n";

  temp = krasnopolsky_temperature(Texo, T_tropo, r_tropo, shape_parameter, false/*shape parameter is in absolute units of km*/);
  chamb_diff_1d atm(rminn,
		    rexoo,
		    rmaxx,
		    rmindiffusionn,
		    nHexo,
		    nCO2rminn,
		    &temp,
		    &H_thermosphere,
		    thermosphere_exosphere::method_rmax_nCO2rmin);
  atm.copy_H_options(H_cross_section_options);
  
  if (atmosphere_fname !="")
    atm.save(atmosphere_fname);

  if (plane_parallel) {
    generate_source_function_plane_parallel(atm,Texo,
					    sourcefn_fname);
  } else {
    generate_source_function_sph_azi_sym(atm,Texo,
					 sourcefn_fname);
  }
}



template <typename A>
void observation_fit::generate_source_function_plane_parallel(A &atmm, const Real &Texo,
							      const string sourcefn_fname/* = ""*/)
{
  bool atmm_spherical = atmm.spherical;
  atmm.spherical = false;
  
  hydrogen_RT_pp.grid.setup_voxels(atmm);
  hydrogen_RT_pp.grid.setup_rays();

  //update the emission density values
  lyman_alpha_pp.define("H Lyman alpha",
			1.0,
			Texo, atmm.sH_lya(Texo),
			atmm,
			&A::n_species_voxel_avg,   &A::Temp_voxel_avg,
			&A::n_absorber_voxel_avg,  &A::sCO2_lya,
			hydrogen_RT_pp.grid.voxels);
  lyman_beta_pp.define("H Lyman beta",
  		       lyman_beta_branching_ratio,
  		       Texo, atmm.sH_lyb(Texo),
  		       atmm,
		       &A::n_species_voxel_avg,   &A::Temp_voxel_avg,
		       &A::n_absorber_voxel_avg,  &A::sCO2_lyb,
  		       hydrogen_RT_pp.grid.voxels);

  atmm.spherical = atmm_spherical;  

  //compute source function on the GPU if compiled with NVCC
#ifdef __CUDACC__
  hydrogen_RT_pp.generate_S_gpu();
#else
  hydrogen_RT_pp.generate_S();
#endif

  if (sourcefn_fname!="")
    hydrogen_RT_pp.save_S(sourcefn_fname);
}

template <typename A>
void observation_fit::generate_source_function_sph_azi_sym(A &atmm, const Real &Texo,
							   const string sourcefn_fname/* = ""*/)
{
  bool change_spherical = false;
  if (atmm.spherical != true) {
    change_spherical = true;
    atmm.spherical = true;
  }

  hydrogen_RT.grid.setup_voxels(atmm);
  hydrogen_RT.grid.setup_rays();


  //update the emission density values
  lyman_alpha.define("H Lyman alpha",
		     1.0,
		     Texo, atmm.sH_lya(Texo),
		     atmm,
		     &A::n_species_voxel_avg,   &A::Temp_voxel_avg,
		     &A::n_absorber_voxel_avg,  &A::sCO2_lya,
		     hydrogen_RT.grid.voxels);
  lyman_beta.define("H Lyman beta",
  		    lyman_beta_branching_ratio,
  		    Texo, atmm.sH_lyb(Texo),
  		    atmm,
		    &A::n_species_voxel_avg,   &A::Temp_voxel_avg,
		    &A::n_absorber_voxel_avg,  &A::sCO2_lyb,
  		    hydrogen_RT.grid.voxels);

  if (change_spherical)
    atmm.spherical = false;    

  //compute source function on the GPU if compiled with NVCC
#ifdef __CUDACC__
  hydrogen_RT.generate_S_gpu();
#else
  hydrogen_RT.generate_S();
#endif

  if (sourcefn_fname!="")
    hydrogen_RT.save_S(sourcefn_fname);
}


void observation_fit::generate_source_function_nH_asym(const Real &nHexo, const Real &Texo,
						       const Real &asym,
						       const string sourcefn_fname/* = ""*/) {
  
  temp = krasnopolsky_temperature(Texo);
  chamb_diff_1d_asymmetric atm_asym(nHexo,CO2_exobase_density,&temp,&H_thermosphere);
  atm_asym.copy_H_options(H_cross_section_options);
  atm_asym.set_asymmetry(asym);

  generate_source_function_sph_azi_sym(atm_asym,Texo,
				       sourcefn_fname);

}

void observation_fit::generate_source_function_temp_asym(const Real &nHavg,
							 const Real &Tnoon, const Real &Tmidnight,
							 const string sourcefn_fname/* = ""*/) {
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
				     /*          Tpower = */2.5);
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
							 const string sourcefn_fname/* = ""*/) {
  std::cout << "nHavg = " << nHavg << "; Tnoon = " << Tnoon << "; Tmidnight = " << Tmidnight <<".\n";
  std::cout << "nCO2rmin = " << nCO2rminn << ".\n";
  std::cout << "T_tropo = " << T_tropo << "; z_tropo = " << (r_tropo-rMars)/1e5 << "; shape_parameter = " << shape_parameter << ".\n";
  std::cout << "Tpower = " << Tpowerr << ".\n";

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
  
  generate_source_function_sph_azi_sym(atm_asym,Tnoon,
				       sourcefn_fname);  
}

void observation_fit::generate_source_function_tabular_atmosphere(const Real rmin, const Real rexo, const Real rmax,
								  const std::vector<doubReal> &alt_nH, const std::vector<doubReal> &log_nH,
								  const std::vector<doubReal> &alt_nCO2, const std::vector<doubReal> &log_nCO2,
								  const std::vector<doubReal> &alt_temp, const std::vector<doubReal> &temp,
								  const bool compute_exosphere/* = false*/,
								  const bool plane_parallel/*= false*/,
								  const string sourcefn_fname/* = ""*/) {

  tabular_1d new_atm_tabular(rmin,rexo,rmax,compute_exosphere);
  new_atm_tabular.copy_H_options(H_cross_section_options);
  atm_tabular = new_atm_tabular;

  atm_tabular.load_log_species_density(alt_nH, log_nH);
  atm_tabular.load_log_absorber_density(alt_nCO2, log_nCO2);
  atm_tabular.load_temperature(alt_temp, temp);

  Real Texo = atm_tabular.Temp(rexo);

  if (plane_parallel) {
    generate_source_function_plane_parallel(atm_tabular, Texo,
					    sourcefn_fname);
  } else {
    generate_source_function_sph_azi_sym(atm_tabular, Texo,
					 sourcefn_fname);
  }
}

void observation_fit::set_use_CO2_absorption(const bool use_CO2_absorption/* = true*/) {
  H_cross_section_options.no_CO2_absorption = !use_CO2_absorption;
  atm_tabular.no_CO2_absorption = !use_CO2_absorption;
}
void observation_fit::set_use_temp_dependent_sH(const bool use_temp_dependent_sH/* = true*/, const Real constant_temp_sH/* = -1*/) {
  H_cross_section_options.temp_dependent_sH = use_temp_dependent_sH;
  atm_tabular.temp_dependent_sH = use_temp_dependent_sH;
  
  H_cross_section_options.constant_temp_sH = constant_temp_sH;
  atm_tabular.constant_temp_sH = constant_temp_sH;
}

void observation_fit::set_sza_method_uniform() {
  hydrogen_RT.grid.szamethod = hydrogen_RT.grid.szamethod_uniform;
}
void observation_fit::set_sza_method_uniform_cos() {
  hydrogen_RT.grid.szamethod = hydrogen_RT.grid.szamethod_uniform_cos;
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
  std::cout << "Simulating O 102.6 nm brightness.\n";
  std::cout << "nOexo = " << nOexo << "; Texo = " << Texo << ".\n";

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
