//los_tracker.h -- track optical depth and brightness along a line of sight

#ifndef __LOS_TRACKER_H
#define __LOS_TRACKER_H

#include "cuda_compatibility.h"
#include "emission.h"
#include "boundaries.h"

template <int N_EMISS>
struct tau_tracker {
  double tau_species_initial[N_EMISS];
  double tau_species_final[N_EMISS];
  double tau_absorber_initial[N_EMISS];
  double tau_absorber_final[N_EMISS];
  double max_tau_species;

  CUDA_CALLABLE_MEMBER
  tau_tracker() {
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++) {
      tau_species_initial[i_emission]=0;
      tau_species_final[i_emission]=0;
      tau_absorber_initial[i_emission]=0;
      tau_absorber_final[i_emission]=0;
    }
    max_tau_species = 0;
  }
  CUDA_CALLABLE_MEMBER
  ~tau_tracker() { }
  CUDA_CALLABLE_MEMBER
  tau_tracker(const tau_tracker<N_EMISS> &copy) {
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++) {
      tau_species_initial[i_emission]=copy.tau_species_initial[i_emission];
      tau_species_final[i_emission]=copy.tau_species_final[i_emission];
      tau_absorber_initial[i_emission]=copy.tau_absorber_initial[i_emission];
      tau_absorber_final[i_emission]=copy.tau_absorber_final[i_emission];
    }
    max_tau_species=copy.max_tau_species;
  }
  CUDA_CALLABLE_MEMBER
  tau_tracker& operator=(const tau_tracker<N_EMISS> &rhs) {
    if(this == &rhs) return *this;
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++) {
      tau_species_initial[i_emission]=rhs.tau_species_initial[i_emission];
      tau_species_final[i_emission]=rhs.tau_species_final[i_emission];
      tau_absorber_initial[i_emission]=rhs.tau_absorber_initial[i_emission];
      tau_absorber_final[i_emission]=rhs.tau_absorber_final[i_emission];
    }
    max_tau_species=rhs.max_tau_species;
    return *this;
  }

  CUDA_CALLABLE_MEMBER
  void reset() {
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++) {
      tau_species_initial[i_emission]=0;
      tau_species_final[i_emission]=0;
      tau_absorber_initial[i_emission]=0;
      tau_absorber_final[i_emission]=0;
    }
  }
  
  CUDA_CALLABLE_MEMBER
  void check_max_tau() {
    for (int i_emission=0; i_emission < N_EMISS; i_emission++) {
      if (tau_species_final[i_emission] > max_tau_species)
	max_tau_species = tau_species_final[i_emission];
    }
  }

  template <typename S, typename E>
  CUDA_CALLABLE_MEMBER
  void update_start(S& stepper,
		    E (&emissions)[N_EMISS]) {
    for (int i_emission=0; i_emission < N_EMISS; i_emission++) {
      tau_species_final[i_emission] = ( tau_species_initial[i_emission]
					+ (emissions[i_emission].dtau_species(stepper.current_voxel)
					   * stepper.pathlength));
      tau_absorber_final[i_emission] = ( tau_absorber_initial[i_emission]
					 + (emissions[i_emission].dtau_absorber(stepper.current_voxel)
					    * stepper.pathlength)); 
    }
  }
  
  CUDA_CALLABLE_MEMBER
  void update_end() {
    for (int i_emission=0; i_emission < N_EMISS; i_emission++) {
      tau_species_initial[i_emission]=tau_species_final[i_emission];
      tau_absorber_initial[i_emission]=tau_absorber_final[i_emission];
    }
    check_max_tau();
  }
};

template <int N_EMISS>
struct brightness_tracker : tau_tracker<N_EMISS> {
  double brightness[N_EMISS];

  CUDA_CALLABLE_MEMBER
  brightness_tracker() : tau_tracker<N_EMISS>() {
    for (int i_emission=0; i_emission < N_EMISS; i_emission++)
       brightness[i_emission]=0;
  }
  CUDA_CALLABLE_MEMBER
  ~brightness_tracker() { }
  CUDA_CALLABLE_MEMBER
  brightness_tracker(const brightness_tracker<N_EMISS> &copy)
    : tau_tracker<N_EMISS>(copy) 
  {
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++)
      brightness[i_emission]=copy.brightness[i_emission];
  }
  CUDA_CALLABLE_MEMBER
  brightness_tracker& operator=(const brightness_tracker<N_EMISS> &rhs) {
    if(this == &rhs) return *this;
    tau_tracker<N_EMISS>::operator=(rhs);
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++)
      brightness[i_emission]=rhs.brightness[i_emission];
    return *this;
  }


  CUDA_CALLABLE_MEMBER
  void reset() {
    tau_tracker<N_EMISS>::reset();
    for (int i_emission=0; i_emission < N_EMISS; i_emission++)
      brightness[i_emission]=0;
  }
};

//interpolation support
template <int N>
class interpolated_values {
public:
  double dtau_species_interp[N];
  double dtau_absorber_interp[N];
  double sourcefn_interp[N];
};

#endif
