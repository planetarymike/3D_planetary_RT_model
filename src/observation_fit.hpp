//observation_fit.hpp -- routines to fit an atmosphere observation

#ifndef __OBSERVATION_FIT_H
#define __OBSERVATION_FIT_H

#include <string> 
#include "Real.hpp"
#include "observation.hpp"
#include "atm/temperature.hpp"
#include "atm/chamberlain_exosphere.hpp"
#include "atm/chamb_diff_1d.hpp"
#include "atm/chamb_diff_1d_asymmetric.hpp"
#include "atm/chamb_diff_temp_asymmetric.hpp"
#include "atm/tabular_1d.hpp"
#include "RT_grid.hpp"
#include "grid/grid_plane_parallel.hpp"
#include "grid/grid_spherical_azimuthally_symmetric.hpp"
#include "emission/singlet_CFR.hpp"
#include "emission/H_lyman_multiplet.hpp"
#include "emission/H_lyman_multiplet_test.hpp"
#include "emission/O_1026.hpp"

class observation_fit {
protected:
  hydrogen_density_parameters H_thermosphere;
  deuterium_density_parameters D_thermosphere;
  oxygen_density_parameters O_thermosphere;
  
  krasnopolsky_temperature temp;
  const Real default_CO2_exobase_density = 2e8;//cm-3
  Real CO2_exobase_density;//cm-3
  //chamb_diff_1d atm;//make sure to use the same exobase alt as in Tconv
  //chamb_diff_1d_asymmetric atm_asym;//make sure to use the same quantities as in atm
  H_cross_sections H_cross_section_options;
  tabular_1d atm_tabular;

  //  // for comparison with Pratik
  // static const int n_radial_boundaries = 90;
  // static const int n_sza_boundaries = 32;/*20 for 10 deg increments with szamethod_uniform*/
  // static const int n_rays_theta = 24;
  // static const int n_rays_phi = 16;

  //standard resolution case
  static const int n_radial_boundaries = 40;
  static const int n_sza_boundaries = 20;/*20 for 10 deg increments with szamethod_uniform*/
  static const int n_rays_theta = 7;
  static const int n_rays_phi = 12;
  typedef plane_parallel_grid<n_radial_boundaries,
			      n_rays_theta> plane_parallel_grid_type;
  plane_parallel_grid_type grid_pp;

  typedef spherical_azimuthally_symmetric_grid<n_radial_boundaries,
					       n_sza_boundaries,
					       n_rays_theta,
					       n_rays_phi> grid_type;
  grid_type grid;

  // define hydrogen emissions
  static const int n_hydrogen_emissions = 2;

  typedef singlet_CFR<plane_parallel_grid_type::n_voxels> hydrogen_emission_type_pp;
  hydrogen_emission_type_pp lyman_alpha_pp, lyman_beta_pp;
  hydrogen_emission_type_pp *hydrogen_emissions_pp[n_hydrogen_emissions] = {&lyman_alpha_pp, &lyman_beta_pp};

  typedef RT_grid<hydrogen_emission_type_pp, n_hydrogen_emissions, plane_parallel_grid_type> H_RT_type_pp;
  H_RT_type_pp hydrogen_RT_pp;

  typedef singlet_CFR<grid_type::n_voxels> hydrogen_emission_type;
  hydrogen_emission_type lyman_alpha, lyman_beta;
  hydrogen_emission_type *hydrogen_emissions[n_hydrogen_emissions] = {&lyman_alpha, &lyman_beta};

  typedef RT_grid<hydrogen_emission_type, n_hydrogen_emissions, grid_type> H_RT_type;
  H_RT_type hydrogen_RT;

  observation<hydrogen_emission_type, n_hydrogen_emissions> hydrogen_obs;



  // Interplanetary H Lyman alpha
  string iph_sfn_fname; // quemerais IPH source function filename
  bool sim_iph;



  // define deuterium emissions
  hydrogen_emission_type_pp D_lyman_alpha_pp, D_lyman_beta_pp;
  hydrogen_emission_type_pp *deuterium_emissions_pp[n_hydrogen_emissions] = {&D_lyman_alpha_pp, &D_lyman_beta_pp};
  H_RT_type_pp deuterium_RT_pp;

  hydrogen_emission_type D_lyman_alpha, D_lyman_beta;
  hydrogen_emission_type *deuterium_emissions[n_hydrogen_emissions] = {&D_lyman_alpha, &D_lyman_beta};
  H_RT_type deuterium_RT;

  observation<hydrogen_emission_type, n_hydrogen_emissions> deuterium_obs;



  // H Lyman alpha multiplet models
  typedef H_lyman_multiplet<grid_type::n_voxels> ly_multiplet_type;
  ly_multiplet_type ly_multiplet;
  ly_multiplet_type *ly_multiplet_emissions[1] = {&ly_multiplet};

  typedef RT_grid<ly_multiplet_type, 1, grid_type> H_RT_ly_multiplet_type;
  H_RT_ly_multiplet_type ly_multiplet_RT;

  observation<ly_multiplet_type, 1> ly_multiplet_obs;


  // H Lyman alpha singlet model in multiplet framework
  typedef H_lyman_singlet<grid_type::n_voxels> ly_singlet_type;
  ly_singlet_type ly_singlet;
  ly_singlet_type *ly_singlet_emissions[1] = {&ly_singlet};

  typedef RT_grid<ly_singlet_type, 1, grid_type> H_RT_ly_singlet_type;
  H_RT_ly_singlet_type ly_singlet_RT;

  observation<ly_singlet_type, 1> ly_singlet_obs;
  

  // define oxygen emissions
  static const int n_oxygen_emissions = 1;

  typedef O_1026_emission<grid_type::n_voxels> oxygen_emission_type;
  oxygen_emission_type oxygen_1026;
  oxygen_emission_type *oxygen_emissions[n_oxygen_emissions] = {&oxygen_1026};

  RT_grid<oxygen_emission_type, n_oxygen_emissions, grid_type> oxygen_RT;

  observation<oxygen_emission_type, n_oxygen_emissions> oxygen_obs;

  // variables used to test sensitivity to density and temperature in each voxel
  bool tweak_H_density;
  vector<int> tweak_H_density_voxel_numbers;
  Real tweak_H_density_factor;

  bool tweak_H_temp;
  vector<int> tweak_H_temp_voxel_numbers;
  Real tweak_H_temp_factor;

public:
  observation_fit(const string iph_sfn_fnamee);

  Temp_converter Tconv;//also takes exobase alt argument

  void add_observation(const std::vector<vector<Real>> &MSO_locations,
		       const std::vector<vector<Real>> &MSO_directions);

  void get_unextincted_iph();
  
  void set_g_factor(vector<Real> &g);

  void simulate_iph(const bool sim_iphh);
  void add_observation_ra_dec(const std::vector<Real> &mars_ecliptic_coords,
			      const std::vector<Real> &RAA,
			      const std::vector<Real> &Decc);
    

  void generate_source_function(const Real &nHexo, const Real &Texo,
				const string atmosphere_fname = "",
				const string sourcefn_fname = "",
				const bool plane_parallel=false,
				const bool deuterium=false);
  void generate_source_function_effv(const Real &nHexo, const Real &effv_exo,
				     const string atmosphere_fname = "",
				     const string sourcefn_fname = "",
				     const bool plane_parallel=false,
				     const bool deuterium=false);
  void generate_source_function_lc(const Real &nHexo, const Real &lc_exo,
				   const string atmosphere_fname = "",
				   const string sourcefn_fname = "",
				   bool plane_parallel=false,
				   const bool deuterium=false);

  void generate_source_function_variable_thermosphere(const Real &nHexo,
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
						      const string atmosphere_fname = "",
						      const string sourcefn_fname = "",
						      const bool plane_parallel=false,
						      const bool deuterium=false);

  template <typename A, typename RT, typename E>
  void generate_source_function_plane_parallel(A &atmm, const Real &Texo,
					       RT &RT_obj,
					       E &lya_obj,
					       E &lyb_obj,
					       const string sourcefn_fname = "")
  {
    bool atmm_spherical = atmm.spherical;
    atmm.spherical = false;
    
    RT_obj.grid.setup_voxels(atmm);
    RT_obj.grid.setup_rays();
    
    //update the emission density values
    lya_obj.define("H Lyman alpha",
		   1.0,
		   Texo, atmm.sH_lya(Texo),
		   atmm,
		   &A::n_species_voxel_avg,   &A::Temp_voxel_avg,
		   &A::n_absorber_voxel_avg,  &A::sCO2_lya,
		   RT_obj.grid.voxels);
    lyb_obj.define("H Lyman beta",
		   lyman_beta_branching_ratio,
		   Texo, atmm.sH_lyb(Texo),
		   atmm,
		   &A::n_species_voxel_avg,   &A::Temp_voxel_avg,
		   &A::n_absorber_voxel_avg,  &A::sCO2_lyb,
		   RT_obj.grid.voxels);

    if (tweak_H_density) {
      lya_obj.tweak_species_density(tweak_H_density_voxel_numbers, tweak_H_density_factor);
      lyb_obj.tweak_species_density(tweak_H_density_voxel_numbers, tweak_H_density_factor);
    }
    if (tweak_H_temp) {
      lya_obj.tweak_species_temp(tweak_H_temp_voxel_numbers, tweak_H_temp_factor);
      lyb_obj.tweak_species_temp(tweak_H_temp_voxel_numbers, tweak_H_temp_factor);
    }
    
    atmm.spherical = atmm_spherical;  
    
    //compute source function on the GPU if compiled with NVCC
#ifdef __CUDACC__
    RT_obj.generate_S_gpu();
#else
    RT_obj.generate_S();
#endif
    
    if (sourcefn_fname!="")
      RT_obj.save_S(sourcefn_fname);
  }
  
  template <typename A, typename RT, typename E>
  void generate_source_function_sph_azi_sym(A &atmm, const Real &Texo,
					    RT &RT_obj,
					    E &lya_obj,
					    E &lyb_obj,
					    const string sourcefn_fname = "")
  {
    bool change_spherical = false;
    if (atmm.spherical != true) {
      change_spherical = true;
      atmm.spherical = true;
    }
    
    RT_obj.grid.setup_voxels(atmm);
    RT_obj.grid.setup_rays();
    
    //update the emission density values
    lya_obj.define("H Lyman alpha",
		   1.0,
		   Texo, atmm.sH_lya(Texo),
		   atmm,
		   &A::n_species_voxel_avg,   &A::Temp_voxel_avg,
		   &A::n_absorber_voxel_avg,  &A::sCO2_lya,
		   RT_obj.grid.voxels);
    lyb_obj.define("H Lyman beta",
		   lyman_beta_branching_ratio,
		   Texo, atmm.sH_lyb(Texo),
		   atmm,
		   &A::n_species_voxel_avg,   &A::Temp_voxel_avg,
		   &A::n_absorber_voxel_avg,  &A::sCO2_lyb,
		   RT_obj.grid.voxels);

    if (tweak_H_density) {
      lya_obj.tweak_species_density(tweak_H_density_voxel_numbers, tweak_H_density_factor);
      lyb_obj.tweak_species_density(tweak_H_density_voxel_numbers, tweak_H_density_factor);
    }
    if (tweak_H_temp) {
      lya_obj.tweak_species_temp(tweak_H_temp_voxel_numbers, tweak_H_temp_factor);
      lyb_obj.tweak_species_temp(tweak_H_temp_voxel_numbers, tweak_H_temp_factor);
    }
    
    if (change_spherical)
      atmm.spherical = false;    
    
    //compute source function on the GPU if compiled with NVCC
#ifdef __CUDACC__
    RT_obj.generate_S_gpu();
#else
    RT_obj.generate_S();
#endif
    
    if (sourcefn_fname!="")
      RT_obj.save_S(sourcefn_fname);
  }


  void generate_source_function_nH_asym(const Real &nHexo, const Real &Texo,
					const Real &asym,
					const string sourcefn_fname = "",
					const bool deuterium=false);
  void generate_source_function_temp_asym(const Real &nHavg,
					  const Real &Tnoon, const Real &Tmidnight,
					  const string sourcefn_fname = "",
					  const bool deuterium=false);

  void generate_source_function_temp_asym(const Real &nHavg,
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
					  const string sourcefn_fname = "",
					  const bool deuterium=false);

  
  void generate_source_function_tabular_atmosphere(const Real rmin, const Real rexo, const Real rmax,
 						   const std::vector<doubReal> &alt_nH, const std::vector<doubReal> &log_nH,
						   const std::vector<doubReal> &alt_nCO2, const std::vector<doubReal> &log_nCO2,
						   const std::vector<doubReal> &alt_temp, const std::vector<doubReal> &temp,
						   const bool compute_exosphere = false,
						   const bool plane_parallel = false,
						   const bool deuterium=false,
						   const string sourcefn_fname="");

  void set_use_CO2_absorption(const bool use_CO2_absorption = true);
  void set_use_temp_dependent_sH(const bool use_temp_dependent_sH = true, const Real constant_temp_sH = -1);

  void set_sza_method_uniform();
  void set_sza_method_uniform_cos();
  
  void reset_H_lya_xsec_coef(const Real xsec_coef = lyman_alpha_line_center_cross_section_coef);
  void reset_H_lyb_xsec_coef(const Real xsec_coef = lyman_beta_line_center_cross_section_coef);
  void reset_CO2_lya_xsec(const Real xsec = CO2_lyman_alpha_absorption_cross_section);
  void reset_CO2_lyb_xsec(const Real xsec = CO2_lyman_beta_absorption_cross_section);

  Real get_CO2_exobase_density();
  void reset_CO2_exobase_density();
  void set_CO2_exobase_density(const double nCO2);

  void save_influence_matrix(const string fname);
  void save_influence_matrix_O_1026(const string fname);

  void set_H_density_tweak(const bool tweak_H_densityy = false);
  void set_H_density_tweak_values(const vector<int> voxels_to_tweak, const Real tweak_factor);
  void set_H_temp_tweak(const bool tweak_H_tempp = false);
  void set_H_temp_tweak_values(const vector<int> voxels_to_tweak, const Real tweak_factor);
  
  std::vector<std::vector<Real>> brightness();
  std::vector<std::vector<Real>> species_col_dens();
  std::vector<std::vector<Real>> tau_species_final();
  std::vector<std::vector<Real>> tau_absorber_final();
  std::vector<std::vector<Real>> iph_brightness_observed();
  std::vector<std::vector<Real>> iph_brightness_unextincted();

  std::vector<std::vector<Real>> D_brightness();
  std::vector<std::vector<Real>> D_col_dens();
  std::vector<std::vector<Real>> tau_D_final();


  void O_1026_generate_source_function(const Real &nOexo,
  				       const Real &Texo,
  				       const Real &solar_brightness_lyman_beta, // 1.69e-3 is a good number for solar minimum
  				       const string atmosphere_fname = "",
  				       const string sourcefn_fname = "");
  std::vector<std::vector<Real>> O_1026_brightness();

  void lyman_multiplet_generate_source_function(const Real &nHexo,
						const Real &Texo,
						// const Real &solar_brightness_lyman_alpha,
						// const Real &solar_brightness_lyman_beta,
						const string atmosphere_fname = "",
						const string sourcefn_fname = "");
  std::vector<std::vector<Real>> lyman_multiplet_brightness();

  void lyman_singlet_generate_source_function(const Real &nHexo,
					      const Real &Texo,
					      // const Real &solar_brightness_lyman_alpha,
					      // const Real &solar_brightness_lyman_beta,
					      const string atmosphere_fname = "",
					      const string sourcefn_fname = "");
  std::vector<std::vector<Real>> lyman_singlet_brightness();
};

#endif
