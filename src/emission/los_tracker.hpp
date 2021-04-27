//los_tracker.h -- track optical depth and brightness along a line of sight

#ifndef __LOS_TRACKER_H
#define __LOS_TRACKER_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "voxel_vector.hpp"
#include "constants.hpp"

struct los_tracker {
  // base class for radiative transfer influence calculations.
  
  // populated by update_start
  // final optical depths
  Real tau_species_final;
  Real tau_absorber_final;
  Real max_tau_species;

  // influence functions
  Real holstein_T_final; // aka Holstein T function, transmission probability
  //  Real holstein_T_initial;
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
    holstein_T_final = 1.0;
    // holstein_T_initial = 1.0;
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

//#include "voigt.hpp"
//CUDA_CALLABLE_MEMBER
//double Voigt(double xx, double sigma, double lg, int r);
// // Voigt needs n_lambda = 60, lambda_max = 12.0

template <bool influence, int N_VOXELS>
struct singlet_CFR_tracker : std::conditional<influence, influence_tracker<N_VOXELS>, brightness_tracker>::type {
  //holstein function tracker with absorption

  typedef typename std::conditional<influence, influence_tracker<N_VOXELS>, brightness_tracker>::type parent;
  static constexpr int n_lambda = 20; //number of wavelength bins (I usually use 20)
  static constexpr Real lambda_max = 4.0; //max wavelength from line center (units of Doppler width) (4.0 is usually enough)
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
      return delta_lambda; //start and end are half-bins
    else
      return REAL(2.0)*delta_lambda;
  }
  CUDA_CALLABLE_MEMBER
  static Real line_shape_function(const int &i_lambda, const Real &T_ratio) {
    // Real lambda0, lambda1;

    // if (i_lambda==0) {
    //   lambda0 = REAL(0.0);
    //   lambda1 = REAL(0.5)*delta_lambda;    
    // } else if (i_lambda==n_lambda-1) {
    //   lambda0 = lambda_max - REAL(0.5)*delta_lambda;
    //   lambda1 = lambda_max;
    // } else {
    //   lambda0 = lambda(i_lambda) - REAL(0.5)*delta_lambda;
    //   lambda1 = lambda(i_lambda) + REAL(0.5)*delta_lambda;
    // }
      
    // return Real(0.5)*(exp(-lambda0*lambda0*T_ratio)
    // 		      + exp(-lambda1*lambda1*T_ratio));

    
    Real lambda2 = lambda(i_lambda);
    lambda2 *= lambda2;
    return exp(-lambda2*T_ratio);

    // return Voigt(lambda(i_lambda),1.0,3.3e-3,3);
  }
  CUDA_CALLABLE_MEMBER
  static Real line_shape_normalization(const Real &T_ratio) {
    return one_over_sqrt_pi*std::sqrt(T_ratio);
    // return REAL(0.3989423);
  }
  CUDA_CALLABLE_MEMBER
  static Real line_shape_function_normalized(const int &i_lambda, const Real &T_ratio) {
    return line_shape_normalization(T_ratio)*line_shape_function(i_lambda,T_ratio);
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
