//los_tracker.h -- track optical depth and brightness along a line of sight

#ifndef __LOS_TRACKER_H
#define __LOS_TRACKER_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "emission.hpp"
#include "voxel_vector.hpp"

struct los_tracker {
  // base class for radiative transfer influence calculations.
  
  // populated by update_start
  // final optical depths
  Real tau_species_final;
  Real tau_absorber_final;
  Real max_tau_species;

  // influence functions
  Real holstein_T_final; // aka Holstein T function, transmission probability
  Real holstein_T_int; // integral of transmission probability across current voxel
  Real holstein_G_int; // coupling coefficient: integral of
		       //   differential tansmission probability across
		       //   current voxel



  CUDA_CALLABLE_MEMBER
  void init() {
    max_tau_species = 0;
  }

  CUDA_CALLABLE_MEMBER
  void reset() {
    tau_species_final  = 0.0;
    tau_absorber_final = 0.0;
    //max_tau_species not reset because we want to track this across
    //all lines of sight
  }
  
  CUDA_CALLABLE_MEMBER
  void update_end() {
    check_max_tau();
  }

  CUDA_CALLABLE_MEMBER
  void exits_bottom() {
    // called when the ray exits the bottom of the atmosphere in
    // brightness calculations
    tau_absorber_final = -1.0;
  }

protected:  
  CUDA_CALLABLE_MEMBER
  void check_max_tau() {
    if (tau_species_final > max_tau_species)
      max_tau_species = tau_species_final;
  }
};

struct brightness_tracker : los_tracker {
  // container for emission brightness
  Real brightness;

  CUDA_CALLABLE_MEMBER
  void reset() {
    los_tracker::reset();
    brightness=0;
  }
};

template <int N_VOXELS>
struct influence_tracker : los_tracker {
  // container for influence coefficients
protected:
  typedef voxel_array<N_VOXELS, 1> vv;

public:
  vv influence;
  
  CUDA_CALLABLE_MEMBER
  void init() {
    los_tracker::init();
    reset_influence();
  }

  CUDA_CALLABLE_MEMBER
  void reset_influence() {
    for (int j_pt = 0; j_pt < N_VOXELS; j_pt++)
      influence(j_pt) = 0.0;
  }

  CUDA_CALLABLE_MEMBER
  void reset() {
    los_tracker::reset();
    reset_influence();
  }
};

template <bool influence, int N_VOXELS>
struct H_lyman_series_tracker : std::conditional<influence, influence_tracker<N_VOXELS>, brightness_tracker>::type {
  //holstein function tracker with absorption

  typedef typename std::conditional<influence, influence_tracker<N_VOXELS>, brightness_tracker>::type parent;
  static constexpr int n_lambda = 10; //number of wavelength bins
  static constexpr Real lambda_max = 4.0; //max wavelength from line center
  static constexpr Real delta_lambda = lambda_max/(n_lambda-1);
  
  Real species_T_ratio_at_origin;
  Real transfer_probability_lambda_initial[n_lambda];//exp(-(tau_species_initial+tau_absorber_initial)) at each lambda

  // it is faster to compute lambda and line shapes than store these
  // (GPUs have more compute bandwidth than kernel memory)
  CUDA_CALLABLE_MEMBER
  static Real lambda(const int &i_lambda) {
    return i_lambda*delta_lambda;
  }
  CUDA_CALLABLE_MEMBER
  static Real weight(const int &i_lambda) {
    if (i_lambda==0 || i_lambda==n_lambda-1)
      return REAL(0.5)*delta_lambda; //start and end are half-bins
    else
      return delta_lambda;
  }
  CUDA_CALLABLE_MEMBER
  static Real line_shape_function(const Real &lambda2, const Real &T_ratio) {
    return exp(-lambda2*T_ratio);
  }
  CUDA_CALLABLE_MEMBER
  static Real line_shape_function_normalized(const Real &lambda2, const Real &T_ratio) {
    return std::sqrt(T_ratio)*line_shape_function(lambda2,T_ratio);
  }  

  using parent::reset;
  CUDA_CALLABLE_MEMBER
  void reset(const Real &T_ratio) {
    reset();
    species_T_ratio_at_origin=T_ratio;

    for (int i_lambda = 0; i_lambda<n_lambda; i_lambda++)
      transfer_probability_lambda_initial[i_lambda] = 1.0;
  }

};


template <bool influence, int N_VOXELS>
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

  static const int n_lines = 6;

  // whether the line is part of the singlet, doublet, or triplet
  static constexpr int multiplet_identity[n_lines] = {1, 2, 2, 3, 3, 3};
  static constexpr int lower_level[n_lines]        = {0, 1, 1, 2, 2, 2};
  static constexpr int upper_level[n_lines]        = {1, 1, 2, 1, 2, 3};

  static constexpr Real line_wavelengths[n_lines]     = {// singlet
							 102.81571 /*nm*/,
							 // doublet
							 102.74313 /*nm*/, 102.74305 /*nm*/,
							 // triplet
							 102.57633 /*nm*/, 102.57626 /*nm*/, 102.57616 /*nm*/};
  
  static constexpr Real line_A[n_lines]               = {// singlet
							 4.22e7 /*s^-1*/,
							 // doublet
							 3.17e7 /*s^-1*/, 5.71e7 /*s^-1*/,
							 // triplet
							 2.11e6 /*s^-1*/, 1.91e7 /*s^-1*/, 7.66e7 /*s^-1*/};
  
  static constexpr Real line_f[n_lines]               = {// singlet
							 2.01e-2,
							 // doublet
							 5.02e-3, 1.51e-2,
							 // triplet
							 2.00e-4, 3.01e-3, 1.69e-2};
  
  static constexpr Real line_branching_ratio[n_lines] = {// singlet
							 0.394,
							 // doublet
							 0.296, 0.533,
							 // triplet
							 0.020, 0.178, 0.713};
  
  // final optical depths for each line
  Real tau_species_final[n_lines];
  Real tau_absorber_final[n_lines];
  Real max_tau_species[n_lines];

  // influence functions
  Real holstein_T_final[n_lines]; // aka Holstein T function, used to
				  // calculate single scattering

  Real holstein_T_int[n_lines];   // integral of T across current
				  // voxel, used to calculate
				  // brightness

  Real holstein_G_int[n_lines][n_lines];   // integral of differential
					   // tansmission, used to
					   // compute internal
					   // scattering
  //   ^^^line blending means this needs to have two dimensions to
  //      account for emission in one line and absorption by another.

  


};


#endif
