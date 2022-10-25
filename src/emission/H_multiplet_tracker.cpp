// H_multiplet_tracker.cpp --- Lyman multiplet tracker including fine structure doublet
#include "H_multiplet_tracker.hpp"

// OK, so the only reason this CPP file exists is to provide a
// compilable version of the CUDA __constant__ arrays needed for the
// GPU version of the O 102.6 nm tracker. The below is copied directly
// from the .hpp file, with a different macro to set the values of the
// __constant__ arrays here to populate the extern directive in the
// hpp file.

namespace H_lyman_multiplet_constants_detail {
  // DECLARE_STATIC_ARRAY comes from cuda_compatability.hpp
  // DECLARE_STATIC_ARRAY comes from cu 
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, multiplet_index   , {0, 0, 1, 1})
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, lower_level_index , {0, 0, 0, 0})
  DECLARE_STATIC_ARRAY_CPP(int, n_lines, upper_level_index , {0, 1, 2, 3})

  // line data from NIST ASD
  
  //  rest wavelength  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_wavelength,          {121.5668237310 /* nm */,
								     121.5673644608,
								     102.572182505,
								     102.572296565})

  //  offset from centroid of multiplet
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_wavelength_offset,   {-2.70365e-4 /* nm */,
								      2.70365e-4,
								     -5.703e-5,
								      5.703e-5})

  //  Einstein A  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_A,                   {6.2648e8, /* s^-1 */
								     6.2649e8,
								     1.6725e8,
								     1.6725e8 }) 
  //  line f-value, used to compute absorption cross section  
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_f,                   {0.2776,
								     0.13881,
								     5.2761e-2,
								     2.6381e-2}) 

  //  line absorption cross section, from sigma_tot = pi*e^2/mc * f
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, line_sigma_total,         {line_f_coeff*line_f_array[0],  // cm2 Hz
								     line_f_coeff*line_f_array[1],
								     line_f_coeff*line_f_array[2],
								     line_f_coeff*line_f_array[3]}) 
  
  //  sum of Einstein A's from upper state to all lower states (including all branching)
  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, upper_state_decay_rate,   {6.2648e8, /* s^-1 */
								     6.2649e8,
								     1.6725e8 + 2.2449e7,
								     1.6725e8 + 2.2449e7})

  DECLARE_STATIC_ARRAY_CPP(Real, n_lines, absorber_xsec,            {6.3e-20, // cm2
								     6.3e-20,
								     3.53e-17,
								     3.52e-17})
}
