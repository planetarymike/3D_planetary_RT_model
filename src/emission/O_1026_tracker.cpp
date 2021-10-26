// O_1026_tracker.cpp --- Oxygen 102.6 nm emission tracker
#include "O_1026_tracker.hpp"

// OK, so the only reason this CPP file exists is to provide a
// compilable version of the CUDA __constant__ arrays needed for the
// GPU version of the O 102.6 nm tracker. The below is copied directly
// from the .hpp file, with a different macro to set the values of the
// __constant__ arrays here to populate the extern directive in the
// hpp file.


namespace O_1026_constants_detail {
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, multiplet_identity, {1, 2, 2, 3, 3, 3})
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, multiplet_index   , {0, 1, 1, 2, 2, 2})
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, lower_level_J     , {0, 1, 1, 2, 2, 2})
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, lower_level_index , {0, 1, 1, 2, 2, 2})
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, upper_level_J     , {1, 1, 2, 1, 2, 3})
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, upper_level_index , {0, 0, 1, 0, 1, 2})

  // line data from W. L. Wiese, J. R. Fuhr, and T. M. Deters, AIP Press, Melville, NY, 532 pp. (1996)
  //                ^^^^^^^^^^^ This paper did a survey of calculations up to that point.
  //                            For these transitions something better than about 5% accuracy
  //                            is expected based on comparison between calculations.
  //                Tayal 2009 "Oscillator strengths for allowed transitions in neutral oxygen"
  //                            reports very similar f-values for these transitions
  //                Cashman+2017 https://ui.adsabs.harvard.edu/abs/2017ApJS..230....8C/abstract
  //                            this paper and the O I study they cite have values perhaps
  //                            5% lower than those used here.
  // the values used here are the same as those used in Meier+1987
  
  //  rest wavelength  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_wavelength,          {// singlet
								     102.81571 /*nm*/,
								     // doublet
								     102.74313 /*nm*/, 102.74305 /*nm*/,
								     // triplet
								     102.57633 /*nm*/, 102.57626 /*nm*/, 102.57616 /*nm*/})
  //  offset from centroid of multiplet  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_wavelength_offset,   {// singlet
								     0.0 /*nm*/, 
								     // doublet
								     4e-5 /*nm*/, -4e-5 /*nm*/,
								     // triplet
								     8e-5 /*nm*/, 1e-5 /*nm*/, -9e-5 /*nm*/})
  //  Einstein A  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_A,                   { // singlet
								     4.22e7 /*s^-1*/,
								     // doublet
								     3.17e7 /*s^-1*/, 5.71e7 /*s^-1*/,
								     // triplet
								     2.11e6 /*s^-1*/, 1.91e7 /*s^-1*/, 7.66e7 /*s^-1*/})
  //  line f-value, used to compute absorption cross section  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_f,                   {// singlet
								     2.01e-2, // unitless
								     // doublet
								     5.02e-3, 1.51e-2,
								     // triplet
								     2.00e-4, 3.01e-3, 1.69e-2})

  //  line absorption cross section, from sigma_tot = pi*e^2/mc * f
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_sigma_total,         {line_f_coeff*line_f_array[0], // cm2 Hz
								     line_f_coeff*line_f_array[1],
								     line_f_coeff*line_f_array[2],
								     line_f_coeff*line_f_array[3],
								     line_f_coeff*line_f_array[4],
								     line_f_coeff*line_f_array[5]})
  
  //  sum of Einstein A's from upper state to all lower states (includes branching to ~1129nm)
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, upper_state_decay_rate,   {// J = 1
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
								     + 3.09e7})
  //  energy of the lower states
  //    note: for atomic O, the J=2 state is the ground state and the
  //          lower J levels increase in energy
  DECLARE_STATIC_ARRAY_CPP(Real, n_lower, lower_state_energy,             {/* J = 0 */ REAL(0.0281416)*erg_per_eV, // erg
									   /* J = 1 */ REAL(0.0196224)*erg_per_eV,
									   /* J = 2 */ REAL(0.0      )*erg_per_eV})
  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lower, lower_state_statistical_weight, { /* J = 0 */ 1, 
									    /* J = 1 */ 3,
									    /* J = 2 */ 5})
}
