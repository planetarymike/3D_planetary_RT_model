//los_tracker.h -- track optical depth and brightness along a line of sight

#ifndef __LOS_TRACKER_H
#define __LOS_TRACKER_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "emission.hpp"
#include "boundaries.hpp"
#include "lineshape_tracker.hpp"
#include <cmath>

template <int N_EMISS>
struct tau_tracker {
  lineshape_tracker line[N_EMISS];
  Real max_tau_species;

  // CUDA_CALLABLE_MEMBER
  // tau_tracker() { }
  // CUDA_CALLABLE_MEMBER
  // ~tau_tracker() { }
  // CUDA_CALLABLE_MEMBER
  // tau_tracker(const tau_tracker<N_EMISS> &copy) {
  //   max_tau_species = copy.max_tau_species;
  //   for (int i_emission = 0; i_emission<N_EMISS; i_emission++)
  //     line[i_emission] = copy.line[i_emission];
  // }
  // CUDA_CALLABLE_MEMBER
  // tau_tracker& operator=(const tau_tracker<N_EMISS> &rhs) {
  //   if(this == &rhs) return *this;
  //   max_tau_species = rhs.max_tau_species;
  //   for (int i_emission = 0; i_emission<N_EMISS; i_emission++)
  //     line[i_emission] = rhs.line[i_emission];
  //   return *this;
  // }

  CUDA_CALLABLE_MEMBER
  void init() {
    //this is not in the constructor so that this object can be
    //declared in CUDA shared memory, which requires an empty
    //constructor
    max_tau_species = 0.0;
    for (int i_emiss=0;i_emiss<N_EMISS;i_emiss++)
      line[i_emiss].init();
  }

  template <typename E>
  CUDA_CALLABLE_MEMBER
  void reset(const E (&emissions)[N_EMISS], const int i_voxel) {
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++)
      line[i_emission].reset(emissions[i_emission].species_T[i_voxel],
			     emissions[i_emission].species_T_ref);
  }
  template <typename E>
  CUDA_CALLABLE_MEMBER
  void reset(const E (&emissions)[N_EMISS]) {
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++)
      line[i_emission].reset(emissions[i_emission].species_T_ref,
			     emissions[i_emission].species_T_ref);
  }
  
  template <typename S, typename E>
  CUDA_CALLABLE_MEMBER
  void update_start(S& stepper,
		    E (&emissions)[N_EMISS]) {
    for (int i_emission=0; i_emission < N_EMISS; i_emission++) 
      line[i_emission].update_start(stepper, emissions[i_emission]);
  }

  CUDA_CALLABLE_MEMBER
  void update_end() {
    for (int i_emission=0; i_emission < N_EMISS; i_emission++) {
      line[i_emission].update_end();
      if (line[i_emission].max_tau_species > max_tau_species)
	max_tau_species = line[i_emission].max_tau_species;
    }
  }
};

template <int N_EMISS, int N_VOXELS>
struct influence_tracker : tau_tracker<N_EMISS> {
  Real influence[N_EMISS][N_VOXELS];

  // CUDA_CALLABLE_MEMBER
  // influence_tracker() { }
  // CUDA_CALLABLE_MEMBER
  // ~influence_tracker() { }
  // CUDA_CALLABLE_MEMBER
  // influence_tracker(const influence_tracker<N_EMISS,N_VOXELS> &copy)
  //   : tau_tracker<N_EMISS>(copy) 
  // {
  //   for (int i_emission=0; i_emission < N_EMISS; i_emission++)
  //     for (int j_pt = 0; j_pt < N_VOXELS; j_pt++)
  // 	influence[i_emission][j_pt]=copy.influence[i_emission][j_pt];
  // }
  // CUDA_CALLABLE_MEMBER
  // influence_tracker& operator=(const influence_tracker<N_EMISS,N_VOXELS> &rhs) {
  //   if(this == &rhs) return *this;
  //   tau_tracker<N_EMISS>::operator=(rhs);
  //   for (int i_emission=0; i_emission < N_EMISS; i_emission++)
  //     for (int j_pt = 0; j_pt < N_VOXELS; j_pt++)
  // 	influence[i_emission][j_pt]=rhs.influence[i_emission][j_pt];
  //   return *this;
  // }

  CUDA_CALLABLE_MEMBER
  void init() {
    //this is not in the constructor so that this object can be
    //declared in CUDA shared memory, which requires an empty
    //constructor
    tau_tracker<N_EMISS>::init();
    for (int i_emission=0; i_emission < N_EMISS; i_emission++)
      for (int j_pt = 0; j_pt < N_VOXELS; j_pt++)
	influence[i_emission][j_pt] = 0.0;
  }

  template <typename E>
  CUDA_CALLABLE_MEMBER
  void reset(const E (&emissions)[N_EMISS],const int i_voxel) {
    tau_tracker<N_EMISS>::reset(emissions, i_voxel);
    for (int i_emission=0; i_emission < N_EMISS; i_emission++)
      for (int j_pt = 0; j_pt < N_VOXELS; j_pt++)
	influence[i_emission][j_pt] = 0.0;
  }

  CUDA_CALLABLE_MEMBER
  influence_tracker<N_EMISS,N_VOXELS> operator+(const influence_tracker<N_EMISS,N_VOXELS> &other) const {
    influence_tracker<N_EMISS,N_VOXELS> ret;
    if (tau_tracker<N_EMISS>::max_tau_species > other.max_tau_species) {
      ret.max_tau_species = tau_tracker<N_EMISS>::max_tau_species;
    } else {
      ret.max_tau_species = other.max_tau_species;
    }
    
    for (int i_emission=0; i_emission < N_EMISS; i_emission++)
      for (int j_pt = 0; j_pt < N_VOXELS; j_pt++)
	ret.influence[i_emission][j_pt]=influence[i_emission][j_pt]+other.influence[i_emission][j_pt];
  
    return ret;
  }

  CUDA_CALLABLE_MEMBER
  influence_tracker<N_EMISS,N_VOXELS> & operator+=(const influence_tracker<N_EMISS,N_VOXELS> &other) {
    if (tau_tracker<N_EMISS>::max_tau_species < other.max_tau_species)
      tau_tracker<N_EMISS>::max_tau_species = other.max_tau_species;
    
    for (int i_emission=0; i_emission < N_EMISS; i_emission++)
      for (int j_pt = 0; j_pt < N_VOXELS; j_pt++)
	influence[i_emission][j_pt]+=other.influence[i_emission][j_pt];
  
    return *this;
  }
};

template <int N_EMISS>
struct brightness_tracker : tau_tracker<N_EMISS> {
  Real brightness[N_EMISS];

  // CUDA_CALLABLE_MEMBER
  // brightness_tracker() { }
  // CUDA_CALLABLE_MEMBER
  // ~brightness_tracker() { }
  // CUDA_CALLABLE_MEMBER
  // brightness_tracker(const brightness_tracker<N_EMISS> &copy)
  //   : tau_tracker<N_EMISS>(copy) 
  // {
  //   for (int i_emission = 0; i_emission<N_EMISS; i_emission++)
  //     brightness[i_emission]=copy.brightness[i_emission];
  // }
  // CUDA_CALLABLE_MEMBER
  // brightness_tracker& operator=(const brightness_tracker<N_EMISS> &rhs) {
  //   if(this == &rhs) return *this;
  //   tau_tracker<N_EMISS>::operator=(rhs);
  //   for (int i_emission = 0; i_emission<N_EMISS; i_emission++)
  //     brightness[i_emission]=rhs.brightness[i_emission];
  //   return *this;
  // }

  CUDA_CALLABLE_MEMBER
  void init() {
    //this is not in the constructor so that this object can be
    //declared in CUDA shared memory, which requires an empty
    //constructor
    tau_tracker<N_EMISS>::init();
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++)
      brightness[i_emission]=0;
  }

  template <typename E>
  CUDA_CALLABLE_MEMBER
  void reset(const E (&emissions)[N_EMISS]) {
    tau_tracker<N_EMISS>::reset(emissions);
    for (int i_emission=0; i_emission < N_EMISS; i_emission++)
      brightness[i_emission]=0;
  }

};

//interpolation support
template <int N>
class interpolated_values {
public:
  Real dtau_species_interp[N];
  Real species_T_interp[N];
  Real dtau_absorber_interp[N];
  Real abs_interp[N];
  Real sourcefn_interp[N];
};


template <int N_EMISS>
struct observed {
  //structure for observed brightnesses and uncertainties

  Real brightness[N_EMISS];
  Real brightness_unc[N_EMISS];

  CUDA_CALLABLE_MEMBER
  observed() {
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++) {
      brightness[i_emission]=0;
      brightness_unc[i_emission]=0;
    }
  }
  CUDA_CALLABLE_MEMBER
  ~observed() { }
  CUDA_CALLABLE_MEMBER
  observed(const observed<N_EMISS> &copy) {
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++) {
      brightness[i_emission]=copy.brightness[i_emission];
      brightness_unc[i_emission]=copy.brightness_unc[i_emission];
    }
  }
  CUDA_CALLABLE_MEMBER
  observed& operator=(const observed<N_EMISS> &rhs) {
    if(this == &rhs) return *this;
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++) {
      brightness[i_emission]=rhs.brightness[i_emission];
      brightness_unc[i_emission]=rhs.brightness_unc[i_emission];
    }
    return *this;
  }
};


#endif
