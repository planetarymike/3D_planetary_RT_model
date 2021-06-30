//los_tracker.h -- track optical depth and brightness along a line of sight

#ifndef __O_1026_TRACKER_H
#define __O_1026_TRACKER_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "voxel_vector.hpp"
#include "constants.hpp"

template <bool is_influence, int N_VOXELS>
struct O_1026_tracker {
  // this is a tracker for the O 102.6 nm multiplet emission.

  // This multiplet has three upper states, three lower states, and
  // six lines. Each upper state is coupled to the others by the close
  // spacing of the lines, which are within 0.1-0.2 Doppler widths
  // and therefore coupled on the lower levels.

  // Each upper state can also radiatively decay to a different lower
  // state, so that each line has a branching ratio less than unity
  
  // We must keep track of each interacting set of lines seperately in
  // absolute wavelength space to compute optical depths and
  // transition probabilities.

  static const int n_lines = 6; // we model a total of 6 lines
  static const int n_multiplets = 3;
  static const int n_lower = 3; // there are three lower and three upper states
  static const int n_upper = 3;

  // upper, lower, and multiplet states of each line
  static constexpr int multiplet_identity[n_lines] = {1, 2, 2, 3, 3, 3};
  static constexpr int multiplet_index[n_lines]    = {0, 1, 1, 2, 2, 2};
  static constexpr int lower_level_J[n_lines]      = {0, 1, 1, 2, 2, 2};
  static constexpr int lower_level_index[n_lines]  = {0, 1, 1, 2, 2, 2};
  static constexpr int upper_level_J[n_lines]      = {1, 1, 2, 1, 2, 3};
  static constexpr int upper_level_index[n_lines]  = {0, 0, 1, 0, 1, 2};


  // line data

  //  rest wavelength
  static constexpr Real line_wavelength[n_lines]            = {// singlet
							       102.81571 /*nm*/,
							       // doublet
							       102.74313 /*nm*/, 102.74305 /*nm*/,
							       // triplet
							       102.57633 /*nm*/, 102.57626 /*nm*/, 102.57616 /*nm*/};
  //  offset from centroid of multiplet
  static constexpr Real line_wavelength_offset[n_lines]     = {// singlet
							       0.0 /*nm*/, 
							       // doublet
							       4e-5 /*nm*/, -4e-5 /*nm*/,
							       // triplet
							       8e-5 /*nm*/, 1e-5 /*nm*/, -9e-5 /*nm*/};

  //  Einstein A
  static constexpr Real line_A[n_lines]                     = { // singlet
							       4.22e7 /*s^-1*/,
							       // doublet
							       3.17e7 /*s^-1*/, 5.71e7 /*s^-1*/,
							       // triplet
							       2.11e6 /*s^-1*/, 1.91e7 /*s^-1*/, 7.66e7 /*s^-1*/};
  //  line f-value, used to compute absorption cross section
  static constexpr Real line_f[n_lines]                     = {// singlet
							       2.01e-2, // unitless
							       // doublet
							       5.02e-3, 1.51e-2,
							       // triplet
							       2.00e-4, 3.01e-3, 1.69e-2};

  //  line absorption cross section, from sigma_tot = pi*e^2/mc * f
  static constexpr Real line_sigma_total[n_lines]           = {line_f_coeff*line_f[0], // cm2 Hz
							       line_f_coeff*line_f[1],
							       line_f_coeff*line_f[2],
							       line_f_coeff*line_f[3],
							       line_f_coeff*line_f[4],
							       line_f_coeff*line_f[5]};

  //  sum of Einstein A's from upper state to all lower states (includes branching to ~1129nm)
  static constexpr Real upper_state_decay_rate[n_upper]     = {// J = 1
							       //   102.6 nm branch
							       2.11e6 /* to J=2*/ + 3.17e7 /* to J=1*/ + 4.22e7 /* to J=0*/
							       //   1129 nm branch
							       + 1.29e7 /* to J=1*/ + 8.6e5 /* to J=2*/ + 1.72e7 /* to J=0*/,

							       // J = 2
							       //   102.6 nm branch
							       1.91e7 /* to J=2*/ + 5.71e7 /* to J=1 */
							       //   1129 nm branch
							       + 2.32e7 /* to J=1*/ + 7.74e6 /* to J=2 */,
							       
							       // J = 3
							       //   102.6 nm branch
							       7.66e7 /* to J=2 */
							       //   1129 nm branch
							       + 3.09e7};

  //  energy of the lower states
  //    note: for atomic O, the J=2 state is the ground state and the
  //          lower J levels increase in energy
  static constexpr Real lower_state_energy[n_lower]              = {/* J = 0 */ REAL(0.0281416)*erg_per_eV, // erg
								    /* J = 1 */ REAL(0.0196224)*erg_per_eV,
								    /* J = 2 */ REAL(0.0      )*erg_per_eV};
  
  static constexpr int lower_state_statistical_weight[n_lower] = { /* J = 0 */ 1, 
								   /* J = 1 */ 3,
								   /* J = 2 */ 5};
  
  
  //  branching ratios for the line, A_ul / sum_L(AuL), includes branching to ~1129nm
  static constexpr Real line_branching_ratio[n_lines]       = { // singlet
							       0.394, 
							       // doublet
							       0.296, 0.533,
							       // triplet
							       0.020, 0.178, 0.713};

  //  CO2 cross section at 102.6 nm
  static constexpr Real co2_xsec = 3.53e-17; //cm^2 


  // set absolute wavelength scale
  //   pick a single normalized wavelength scale for all lines and multiplets
  //     (this slightly simplifies things, in reality there are
  //      factors of ~10^-4 difference between the different lines,
  //      resulting from different rest wavelengths).
  static constexpr Real doppler_width_reference_T          = 200; // K, should be within a factor of ~2 of expected atmospheric temp
  static constexpr Real doppler_width_reference_lambda     = line_wavelength[5]; // nm, smallest value of lambda
  static constexpr Real doppler_width_reference_velocity   = constexpr_sqrt(2*kB*doppler_width_reference_T/(16*mH)); // cm/s, velocity dispersion
  static constexpr Real doppler_width_wavelength_reference = (doppler_width_reference_lambda
							      * doppler_width_reference_velocity
							      / clight); // nm, doppler width in wavelength, absolute wavelength scale
  static constexpr Real doppler_width_frequency_reference  = (1.0/(doppler_width_reference_lambda*1e-7)
							      * doppler_width_reference_velocity); // Hz, frequency Doppler width
  static constexpr Real normalization = (one_over_sqrt_pi / doppler_width_frequency_reference); // Hz^-1, lineshape normalization

  // normalized wavelength offsets from line center
  static constexpr Real line_wavelength_offset_normalized[n_lines] = {line_wavelength_offset[0] / doppler_width_wavelength_reference,
								      line_wavelength_offset[1] / doppler_width_wavelength_reference,
								      line_wavelength_offset[2] / doppler_width_wavelength_reference,
								      line_wavelength_offset[3] / doppler_width_wavelength_reference,
								      line_wavelength_offset[4] / doppler_width_wavelength_reference,
								      line_wavelength_offset[5] / doppler_width_wavelength_reference};

  // Tracked Variables

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
			    double>::type influence[n_upper];
  
  // keep track of origin temperature and density for computing influence coefficients
  Real species_T_at_origin;
  Real species_density_at_origin[n_lower];
  
  // we carry one array of transmission probabilities per multiplet
  // each wavelength array is centered on the mean wavelength of the multiplet emission.
  // (line shapes need to incorporate the offset from this mean 
  static constexpr int n_lambda = 41; //number of wavelength bins
  static constexpr Real lambda_max = 8.0; //max wavelength from line center (dimensionless, units of wavelength Doppler width at Tref)
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
  static Real line_shape_function(const int &i_line,
				  const int &i_lambda,
				  const Real &T) {
    Real lambda2 = (lambda(i_lambda)-line_wavelength_offset_normalized[i_line]);
    lambda2 *= lambda2;
    return exp(-lambda2*doppler_width_reference_T/T);
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

    for (int i_lower = 0; i_lower<n_lower; i_lower++) 
      species_density_at_origin[i_lower] = density_at_origin[i_lower];
    
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
