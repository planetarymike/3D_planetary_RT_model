//los_tracker.h -- track optical depth and brightness along a line of sight

#ifndef __LOS_TRACKER_H
#define __LOS_TRACKER_H

#include <vector>
#include "emission.h"
#include "boundaries.h"

using std::vector;

template <int N_EMISS>
struct tau_tracker {
  double tau_species_initial[N_EMISS];
  double tau_species_final[N_EMISS];
  double tau_absorber_initial[N_EMISS];
  double tau_absorber_final[N_EMISS];
  double max_tau_species;

  tau_tracker() {
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++) {
      tau_species_initial[i_emission]=0;
      tau_species_final[i_emission]=0;
      tau_absorber_initial[i_emission]=0;
      tau_absorber_final[i_emission]=0;
    }
    max_tau_species = 0;
  }

  void reset() {
    for (int i_emission = 0; i_emission<N_EMISS; i_emission++) {
      tau_species_initial[i_emission]=0;
      tau_species_final[i_emission]=0;
      tau_absorber_initial[i_emission]=0;
      tau_absorber_final[i_emission]=0;
    }
  }
  
  void check_max_tau() {
    for (int i_emission=0; i_emission < N_EMISS; i_emission++) {
      if (tau_species_final[i_emission] > max_tau_species)
	max_tau_species = tau_species_final[i_emission];
    }
  }

  template <typename S, typename E>
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

  brightness_tracker() : tau_tracker<N_EMISS>() {
    for (int i_emission=0; i_emission < N_EMISS; i_emission++)
       brightness[i_emission]=0;
  }

  void reset() {
    tau_tracker<N_EMISS>::reset();
    for (int i_emission=0; i_emission < N_EMISS; i_emission++)
      brightness[i_emission]=0;
  }
};


#endif
