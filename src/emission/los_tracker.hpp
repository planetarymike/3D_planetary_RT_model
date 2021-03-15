//los_tracker.h -- track optical depth and brightness along a line of sight

#ifndef __LOS_TRACKER_H
#define __LOS_TRACKER_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "emission.hpp"

struct los_tracker {
  // base class for radiative transfer influence calculations.
  
  // populated by update_start_implicit
  // final optical depths
  Real tau_species_final;
  Real tau_absorber_final;
  Real max_tau_species;

  // influence functions
  Real holstein_T_final; // aka Holstein T function
  Real holstein_T_int;   // integral of T across current voxel
  Real holstein_G_int;   // integral of differential tansmission
                         //   probability across current voxel

  CUDA_CALLABLE_MEMBER
  void init() {
    max_tau_species = 0;
  }

  CUDA_CALLABLE_MEMBER
  void reset() {
    tau_species_final=0;
    tau_absorber_final=0;
    holstein_T_final=1.0;
    //max_tau_species not reset because we want to track this across
    //all lines of sight
  }
  
  CUDA_CALLABLE_MEMBER
  void update_end() {
    check_max_tau();
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

  using los_tracker::init;

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
  static const int n_voxels = N_VOXELS;

public:
  Real influence[n_voxels];
  
  CUDA_CALLABLE_MEMBER
  void init() {
    los_tracker::init();
    reset_influence();
  }

  CUDA_CALLABLE_MEMBER
  void reset_influence() {
    for (int j_pt = 0; j_pt < n_voxels; j_pt++)
      influence[j_pt] = 0.0;
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

#endif
