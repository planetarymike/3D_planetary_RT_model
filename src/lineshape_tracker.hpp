//lineshape_tracker.hpp --- holstein integrals computed JIT as lines of sight are traversed
#ifndef __LINESHAPE_TRACKER_H
#define __LINESHAPE_TRACKER_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "emission.hpp"
#include "boundaries.hpp"
//#include "shizgal_quadrature_pts.hpp"

//#define N_SHIZGAL_LAMBDA 16
//#define SHIZGAL_LAMBDAS shizgal_lambdas_16
//#define SHIZGAL_WEIGHTS shizgal_weights_16

struct lineshape_tracker {
  //holstein function tracker with absorption

protected:
  static const int n_lambda=12;// = N_SHIZGAL_LAMBDA;
  static constexpr Real lambda_max=5.0;
  //Real lambda[n_lambda];// assigned in constructor
  //Real weight[n_lambda];// ^
  //Real lambda2[n_lambda];//square of lambda points
  //Real weightfn[n_lambda];//weighting function for this integration rule ( exp(-l**2) for shizgal )

  Real lineshape_at_origin[n_lambda];//line shape at integral origin
  //Real lineshape[n_lambda];//lineshape in current voxel

  Real tau_species_lambda_initial[n_lambda];//tau_species at each lambda
  //Real tau_species_lambda_final[n_lambda];

  Real transfer_probability_lambda_initial[n_lambda];//exp(-(tau_species_initial+tau_absorber_initial)) at each lambda
  //Real transfer_probability_lambda_voxel[n_lambda];//above for just the optical depth across the voxel
  //Real transfer_probability_lambda_final[n_lambda];//above for just the final optical depth in this cell

  CUDA_CALLABLE_MEMBER
  Real lambda(int i_lambda);
  CUDA_CALLABLE_MEMBER
  Real weight(int i_lambda);
  
public:
  Real tau_species_initial;//line center species optical depth
  Real tau_species_final;

  Real tau_absorber_initial;
  Real tau_absorber_final;

  Real max_tau_species;

  Real holstein_T_initial;//holT at tau_initial for current voxel
  Real holstein_T_final;
  Real holstein_T_int;//integral across current voxel
  Real holstein_G_int;//integral across current voxel
  
  // CUDA_CALLABLE_MEMBER
  // lineshape_tracker();
  // CUDA_CALLABLE_MEMBER
  // ~lineshape_tracker();
  // CUDA_CALLABLE_MEMBER
  // lineshape_tracker(const lineshape_tracker &copy);
  // CUDA_CALLABLE_MEMBER
  // lineshape_tracker& operator=(const lineshape_tracker &rhs);

  CUDA_CALLABLE_MEMBER
  void init();
  
  CUDA_CALLABLE_MEMBER
  void reset(const Real &T, const Real &T_ref);

  CUDA_CALLABLE_MEMBER
  void check_max_tau();
  
  template <typename S, typename E>
  CUDA_CALLABLE_MEMBER
  void update_start(const S &stepper,
		    const E &emission) {
    Real T             = emission.species_T(stepper.current_voxel);
    Real T_ref         = emission.species_T_ref;
    Real dtau_species  = emission.dtau_species(stepper.current_voxel);
    Real dtau_absorber = emission.dtau_absorber(stepper.current_voxel);
    Real abs           = emission.abs(stepper.current_voxel);
    Real pathlength    = stepper.pathlength;

    update_start(T, T_ref,
		 dtau_species,
		 dtau_absorber,
		 abs,
		 pathlength);
  }
  

  CUDA_CALLABLE_MEMBER
  void update_start(const Real &T, const Real &T_ref,
		    const Real &dtau_species,
		    const Real &dtau_absorber,
		    const Real &abs,
		    const Real &pathlength);

  CUDA_CALLABLE_MEMBER
  void update_end();
};

#endif
