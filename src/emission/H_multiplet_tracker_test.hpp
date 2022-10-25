// H_multiplet_tracker_test.hpp --- Hydrogen line test tracker for checking consistency of multiplet and singlet codes

#ifndef __H_MULTIPLET_TRACKER_TEST_H
#define __H_MULTIPLET_TRACKER_TEST_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "voxel_vector.hpp"
#include "constants.hpp"

namespace H_lyman_alpha_singlet_constants_detail {
  static constexpr int n_lines = 1; // Lyman alpha and Lyman beta
  static constexpr int n_multiplets = 1;
  static constexpr int n_lower = 1; 
  static constexpr int n_upper = 1;

  static constexpr Real co2_xsec = CO2_lyman_alpha_absorption_cross_section;

  // DECLARE_STATIC_ARRAY comes from cu 
  DECLARE_STATIC_ARRAY_HPP(int, n_lines, multiplet_index   , {0})
  DECLARE_STATIC_ARRAY_HPP(int, n_lines, lower_level_index , {0})
  DECLARE_STATIC_ARRAY_HPP(int, n_lines, upper_level_index , {0})

  // line data from NIST ASD
  
  //  rest wavelength  
  DECLARE_STATIC_ARRAY_HPP(Real, n_lines, line_wavelength,          {lyman_alpha_lambda*1e7}) // nm
  //  Einstein A  
  DECLARE_STATIC_ARRAY_HPP(Real, n_lines, line_A,                   {6.2648e8 /* s^-1 */}) // = one of hyperfine A's (they are nearly identical)
                                                                                           // technically should use f-value weighted average
  //  line f-value, used to compute absorption cross section  
  DECLARE_STATIC_ARRAY_HPP(Real, n_lines, line_f,                   {lyman_alpha_f}) // = sum of hyperfine f's

  //  line absorption cross section, from sigma_tot = pi*e^2/mc * f
  DECLARE_STATIC_ARRAY_HPP(Real, n_lines, line_sigma_total,         {line_f_coeff*line_f_array[0]}) // cm2 Hz
  
  //  sum of Einstein A's from upper state to all lower states (including all branching)
  DECLARE_STATIC_ARRAY_HPP(Real, n_lines, upper_state_decay_rate,   {6.2648e8 /* s^-1 */})
}

namespace H_lyman_beta_singlet_constants_detail {
  static constexpr int n_lines = 1; // Lyman alpha and Lyman beta
  static constexpr int n_multiplets = 1;
  static constexpr int n_lower = 1; 
  static constexpr int n_upper = 1;

  static constexpr Real co2_xsec = CO2_lyman_beta_absorption_cross_section;

  // DECLARE_STATIC_ARRAY comes from cu 
  DECLARE_STATIC_ARRAY_HPP(int, n_lines, multiplet_index   , {0})
  DECLARE_STATIC_ARRAY_HPP(int, n_lines, lower_level_index , {0})
  DECLARE_STATIC_ARRAY_HPP(int, n_lines, upper_level_index , {0})

  // line data from NIST ASD

  // note: although Lyman beta pumps the 2s state of hydrogen via
  // Balmer alpha emission, 2s decays by correlated two-photon
  // emission with an effective A of ~ 8.2/s
  // (https://ui.adsabs.harvard.edu/abs/1984A%26A...138..495N/abstract),
  // yielding a g-factor of ~6e-8 ph/s/molecule for solar moderate
  // conditions. This is ~1e5 times smaller than the Lyman alpha
  // g-factor, so for a column emitting 2 kR of Lyman alpha only 0.02 R
  // of two-photon emisison would be produced. This emission is spread
  // out over ~80-100 nm starting at Lyman alpha and continuing to
  // >200 nm, so the specific intensity at each wavelength would be
  // ~2e-4 R / nm, i.e. not detectable.
  
  //  rest wavelength  
  DECLARE_STATIC_ARRAY_HPP(Real, n_lines, line_wavelength,          {lyman_beta_lambda*1e7}) // nm
  //  Einstein A  
  DECLARE_STATIC_ARRAY_HPP(Real, n_lines, line_A,                   {1.6725e8 /* s^-1 */}) // = one of hyperfine A's (they are nearly identical)
                                                                                           // technically should use f-value weighted average
  //  line f-value, used to compute absorption cross section  
  DECLARE_STATIC_ARRAY_HPP(Real, n_lines, line_f,                   {lyman_beta_f}) // = sum of hyperfine f's

  //  line absorption cross section, from sigma_tot = pi*e^2/mc * f
  DECLARE_STATIC_ARRAY_HPP(Real, n_lines, line_sigma_total,         {line_f_coeff*line_f_array[0]}) // cm2 Hz
  
  //  sum of Einstein A's from upper state to all lower states (including all branching)
  DECLARE_STATIC_ARRAY_HPP(Real, n_lines, upper_state_decay_rate,   {1.6725e8 /* s^-1 */ 
								     +
								     // Balmer alpha
								     2.2449e7 }) 
}

template <bool is_influence, int N_VOXELS>
struct H_lyman_alpha_singlet_test_tracker {
  static const int n_lines      = H_lyman_alpha_singlet_constants_detail::n_lines;
  static const int n_multiplets = H_lyman_alpha_singlet_constants_detail::n_multiplets;
  static const int n_lower      = H_lyman_alpha_singlet_constants_detail::n_lower;
  static const int n_upper      = H_lyman_alpha_singlet_constants_detail::n_upper;

  // import the static array data as member functions that can be called
  CUDA_STATIC_ARRAY_MEMBER(H_lyman_alpha_singlet_constants_detail,  int, n_lines, multiplet_index   )
  CUDA_STATIC_ARRAY_MEMBER(H_lyman_alpha_singlet_constants_detail,  int, n_lines, lower_level_index )
  CUDA_STATIC_ARRAY_MEMBER(H_lyman_alpha_singlet_constants_detail,  int, n_lines, upper_level_index )

  CUDA_STATIC_ARRAY_MEMBER(H_lyman_alpha_singlet_constants_detail, Real, n_lines, line_wavelength               )
  CUDA_STATIC_ARRAY_MEMBER(H_lyman_alpha_singlet_constants_detail, Real, n_lines, line_A                        )
  CUDA_STATIC_ARRAY_MEMBER(H_lyman_alpha_singlet_constants_detail, Real, n_lines, line_f                        )
  CUDA_STATIC_ARRAY_MEMBER(H_lyman_alpha_singlet_constants_detail, Real, n_lines, line_sigma_total              )
  CUDA_STATIC_ARRAY_MEMBER(H_lyman_alpha_singlet_constants_detail, Real, n_lines, upper_state_decay_rate        )
  
  //  CO2 absorption cross section
  static constexpr Real co2_xsec = H_lyman_alpha_singlet_constants_detail::co2_xsec; //cm^2 
  
  // set absolute wavelength scale
  //   pick a single normalized wavelength scale for all lines and multiplets
  static constexpr Real doppler_width_reference_T          = 200; // K, should be within a factor of ~2 of expected atmospheric temp
  static constexpr Real doppler_width_reference_lambda     = H_lyman_alpha_singlet_constants_detail::line_wavelength_array[0]; // nm, smallest value of lambda
  static constexpr Real doppler_width_reference_velocity   = constexpr_sqrt(2*kB*doppler_width_reference_T/mH); // cm/s, velocity dispersion
  static constexpr Real doppler_width_wavelength_reference = (doppler_width_reference_lambda
							      * doppler_width_reference_velocity
							      / clight); // nm, doppler width in wavelength, absolute wavelength scale
  static constexpr Real doppler_width_frequency_reference  = (1.0/(doppler_width_reference_lambda*1e-7)
							      * doppler_width_reference_velocity); // Hz, frequency Doppler width
  static constexpr Real normalization = (one_over_sqrt_pi / doppler_width_frequency_reference); // Hz^-1, lineshape normalization
  
  // Tracked Variables

  // species column density
  Real species_col_dens[n_lower];
  
  // optical depths at line center, counting each line individually
  // even though they overlap
  Real tau_species_final[n_lines];
  Real tau_absorber_final;
  Real max_tau_species;

  // influence functions, one for each upper state
  Real holstein_T_final[n_lines]; // aka Holstein T function, transmission probability
  Real holstein_T_int[n_lines]; // integral of transmission probability across current voxel
  Real holstein_G_int[n_upper][n_upper]; /* coupling coefficient:
                                              integral of differential
                                              tansmission probability
                                              across current voxel */

  Real brightness[n_lines];

  // if we're computing influence coefficients, carry a voxel_arrray to track these
  typename std::conditional<is_influence,
			    voxel_array<N_VOXELS, n_upper>,
			    Real>::type influence[n_upper];
  
  // keep track of origin temperature and density for computing influence coefficients
  Real species_T_at_origin;
  Real species_density_at_origin[n_lower];
  
  // we carry one array of transmission probabilities per multiplet
  // each wavelength array is centered on the mean wavelength of the multiplet emission.
  // (line shapes need to incorporate the offset from this mean)
  static constexpr int n_lambda = 41; //number of wavelength bins
  static constexpr Real lambda_max = 4.0; //max wavelength from line center (dimensionless, units of wavelength Doppler width at Tref)
  static constexpr Real delta_lambda = 2*lambda_max/(n_lambda-1); // dimensionless, fraction of wavelength Doppler width at Tref

  // transfer probability as a function of wavelength, one for each multiplet and each wavelength
  Real transfer_probability_lambda_initial[n_multiplets][n_lambda];


  // functions to deal with line shapes
  CUDA_CALLABLE_MEMBER
  static Real lambda(const int &i_lambda) {
    return (-lambda_max + i_lambda*delta_lambda);
  }
  CUDA_CALLABLE_MEMBER
  static Real weight(__attribute__((unused)) const int &i_lambda) {
    return delta_lambda*doppler_width_frequency_reference; // Hz
  }
  CUDA_CALLABLE_MEMBER
  static Real line_shape_function(__attribute__((unused)) const int &i_line,
				  const int &i_lambda,
				  const Real &T) {
    Real lambda2 = lambda(i_lambda);
    lambda2 *= lambda2;
    lambda2 = lambda2*doppler_width_reference_T/T;
    lambda2 = exp(-lambda2);
    return lambda2;
  }
  CUDA_CALLABLE_MEMBER
  static Real line_shape_normalization(const Real &T) { 
    return normalization*sqrt(doppler_width_reference_T/T); // Hz-1
  }
  CUDA_CALLABLE_MEMBER
  static Real line_shape_function_normalized(const int &i_line, const int &i_lambda, const Real &T) {
    return line_shape_normalization(T)*line_shape_function(i_line,i_lambda,T);
  }  

  // RT functions required for interaction with rest of code
  // TODO: replace copy/paste here with something more elegant
  CUDA_CALLABLE_MEMBER
  void init() {
    max_tau_species = 0.0;
  }
  
  CUDA_CALLABLE_MEMBER
  void reset(const Real T_at_origin, const Real (&density_at_origin)[n_lower]) {
    species_T_at_origin=T_at_origin;

    for (int i_line=0;i_line<n_lines;i_line++) {
      tau_species_final[i_line] = 0.0;
      tau_absorber_final = 0.0;
      brightness[i_line] = 0.0;
    }
    // max_tau_species not reset because we want to track this across
    // all lines of sight

    for (int i_multiplet=0; i_multiplet<n_multiplets; i_multiplet++)
      for (int i_lambda = 0; i_lambda<n_lambda; i_lambda++)
	transfer_probability_lambda_initial[i_multiplet][i_lambda] = 1.0;

    for (int i_lower = 0; i_lower<n_lower; i_lower++) {
      species_density_at_origin[i_lower] = density_at_origin[i_lower];
      species_col_dens[i_lower] = 0.0;
    }
    
    for (int i_upper = 0; i_upper<n_upper; i_upper++)
      influence[i_upper] = 0.0;
  }
  
  CUDA_CALLABLE_MEMBER
  void update_end() {
    check_max_tau();
  }
  
  CUDA_CALLABLE_MEMBER
  void exits_bottom() {
    for (int i_line=0;i_line<n_lines;i_line++)
      tau_absorber_final = -1.0;
    // called when the ray exits the bottom of the atmosphere in
    // brightness calculations
  }

private:
  CUDA_CALLABLE_MEMBER
  void check_max_tau() {
    for (int i_line=0;i_line<n_lines;i_line++)
      if (tau_species_final[i_line] > max_tau_species)
	max_tau_species = tau_species_final[i_line];
  }
};


#endif
