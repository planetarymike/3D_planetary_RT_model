// H_multiplet_tracker_test.cpp --- Hydrogen line test tracker for checking consistency of multiplet and singlet codes
#include "H_multiplet_tracker_test.hpp"

// OK, so the only reason this CPP file exists is to provide a
// compilable version of the CUDA __constant__ arrays needed for the
// GPU version of the O 102.6 nm tracker. The below is copied directly
// from the .hpp file, with a different macro to set the values of the
// __constant__ arrays here to populate the extern directive in the
// hpp file.

namespace H_lyman_alpha_singlet_constants_detail {
   // DECLARE_STATIC_ARRAY comes from cuda_compatability.hpp
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, multiplet_index   , {0})
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, lower_level_index , {0})
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, upper_level_index , {0})

  // line data from NIST ASD
  
  //  rest wavelength  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_wavelength,          {lyman_alpha_lambda})
  //  Einstein A  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_A,                   {6.2648e8 /* s^-1 */}) // = one of hyperfine A's (they are nearly identical)

  //  line f-value, used to compute absorption cross section  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_f,                   {0.41641}) // = sum of hyperfine f's

  //  line absorption cross section, from sigma_tot = pi*e^2/mc * f
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_sigma_total,         {line_f_coeff*line_f_array[0]}) // cm2 Hz
  
  //  sum of Einstein A's from upper state to all lower states (including all branching)
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, upper_state_decay_rate,   {6.2648e8 /* s^-1 */})
}

namespace H_lyman_beta_singlet_constants_detail {
   // DECLARE_STATIC_ARRAY comes from cuda_compatability.hpp
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, multiplet_index   , {0})
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, lower_level_index , {0})
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, upper_level_index , {0})

  // line data from NIST ASD
  
  //  rest wavelength  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_wavelength,          {lyman_beta_lambda})
  //  Einstein A  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_A,                   {1.6725e8 /* s^-1 */}) // = one of hyperfine A's (they are nearly identical)

  //  line f-value, used to compute absorption cross section  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_f,                   {0.079142}) // = sum of hyperfine f's

  //  line absorption cross section, from sigma_tot = pi*e^2/mc * f
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_sigma_total,         {line_f_coeff*line_f_array[0]}) // cm2 Hz
  
  //  sum of Einstein A's from upper state to all lower states (including all branching)
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, upper_state_decay_rate,   {1.6725e8 /* s^-1 */ 
								     +
								     // Balmer alpha
								     2.2449e7 })
}
