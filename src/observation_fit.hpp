//observation_fit.hpp -- routines to fit an atmosphere observation

#ifndef __OBSERVATION_FIT_H
#define __OBSERVATION_FIT_H

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

class observation_fit {
protected:
  krasnopolsky_temperature temp;
  const Real CO2_exobase_density = 2e8;//cm-3
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
  static const int n_rays_theta = 6;
  static const int n_rays_phi = 12;
  typedef plane_parallel_grid<n_radial_boundaries,
			      n_rays_theta> plane_parallel_grid_type;
  plane_parallel_grid_type grid_pp;

  typedef spherical_azimuthally_symmetric_grid<n_radial_boundaries,
					       n_sza_boundaries,
					       n_rays_theta,
					       n_rays_phi> grid_type;
  grid_type grid;

  // define emissions
  static const int n_emissions = 2;

  typedef singlet_CFR<plane_parallel_grid_type::n_voxels> emission_type_pp;
  emission_type_pp lyman_alpha_pp, lyman_beta_pp;
  emission_type_pp *emissions_pp[n_emissions] = {&lyman_alpha_pp, &lyman_beta_pp};

  RT_grid<emission_type_pp, n_emissions, plane_parallel_grid_type> RT_pp;

  typedef singlet_CFR<grid_type::n_voxels> emission_type;
  emission_type lyman_alpha, lyman_beta;
  emission_type *emissions[n_emissions] = {&lyman_alpha, &lyman_beta};

  RT_grid<emission_type, n_emissions, grid_type> RT;

  observation<emission_type, n_emissions> obs;

  bool sim_iph;

public:
  observation_fit();

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
				bool plane_parallel=false);
  void generate_source_function_effv(const Real &nHexo, const Real &effv_exo,
				     const string atmosphere_fname = "",
				     const string sourcefn_fname = "",
				     bool plane_parallel=false);
  void generate_source_function_lc(const Real &nHexo, const Real &lc_exo,
				   const string atmosphere_fname = "",
				   const string sourcefn_fname = "",
				   bool plane_parallel=false);

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
						      bool plane_parallel=false);




  template <typename A>
  void generate_source_function_plane_parallel(A &atm, const Real &Texo,
					       const string sourcefn_fname = "");
  template <typename A>
  void generate_source_function_sph_azi_sym(A &atm, const Real &Texo,
					    const string sourcefn_fname = "");
  

  void generate_source_function_nH_asym(const Real &nHexo, const Real &Texo,
					const Real &asym,
					const string sourcefn_fname = "");
  void generate_source_function_temp_asym(const Real &nHavg,
					  const Real &Tnoon, const Real &Tmidnight,
					  const string sourcefn_fname = "");

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
					  const string sourcefn_fname = "");

  
  void generate_source_function_tabular_atmosphere(const Real rmin, const Real rexo, const Real rmax,
 						   const std::vector<double> &alt_nH, const std::vector<double> &log_nH,
						   const std::vector<double> &alt_nCO2, const std::vector<double> &log_nCO2,
						   const std::vector<double> &alt_temp, const std::vector<double> &temp,
						   const bool compute_exosphere = false,
						   const bool plane_parallel = false,
						   const string sourcefn_fname="");

  void set_use_CO2_absorption(const bool use_CO2_absorption = true);
  void set_use_temp_dependent_sH(const bool use_temp_dependent_sH = true, const Real constant_temp_sH = -1);

  void set_sza_method_uniform();
  void set_sza_method_uniform_cos();
  
  void reset_H_lya_xsec_coef(const Real xsec_coef = lyman_alpha_line_center_cross_section_coef);
  void reset_H_lyb_xsec_coef(const Real xsec_coef = lyman_beta_line_center_cross_section_coef);
  void reset_CO2_lya_xsec(const Real xsec = CO2_lyman_alpha_absorption_cross_section);
  void reset_CO2_lyb_xsec(const Real xsec = CO2_lyman_beta_absorption_cross_section);
  
  std::vector<std::vector<Real>> brightness();
  std::vector<std::vector<Real>> tau_species_final();
  std::vector<std::vector<Real>> tau_absorber_final();
  std::vector<std::vector<Real>> iph_brightness_observed();
  std::vector<std::vector<Real>> iph_brightness_unextincted();

  std::vector<Real> likelihood_and_derivatives(const Real &nHexo, const Real &Texo);
  void logl();
  void logl_gpu();
};

#endif
