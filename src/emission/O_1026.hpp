//O_1026.hpp --- routines to compute O 102.6 nm emission

// most of the code to compute multiplet emission is generic, and
// stored in multiplet_CFR_emission.hpp

// the routines to define the emission densities / temperatures are
// specific to O 102.6 and stored here, so are the solar excitation
// calculations


#ifndef __O_1026_h
#define __O_1026_h

#include <boost/type_traits/type_identity.hpp> //for type deduction in define
#include "emission_voxels.hpp"
#include "atmo_vec.hpp"
#include "multiplet_CFR_emission.hpp"
#include "O_1026_tracker.hpp"

template <int N_VOXELS>
struct O_1026_emission : multiplet_CFR_emission<N_VOXELS,
						/*emission_type = */ O_1026_emission<N_VOXELS>,
						/*los_tracker_type = */ O_1026_tracker> {
protected:
  typedef multiplet_CFR_emission<N_VOXELS,
				 /*emission_type = */ O_1026_emission<N_VOXELS>,
				 /*los_tracker_type = */ O_1026_tracker> parent;
  friend parent;

public:
  // wavelength info and line shape routines are stored in tracker
  // object b/c this compiles to code that is 30% faster, not sure why
  using parent::n_lines;
  using parent::n_multiplets;
  using parent::n_lambda;
  
  template <bool influence>
  using los = O_1026_tracker<influence, N_VOXELS>;
  using typename parent::brightness_tracker;
  using typename parent::influence_tracker;

protected:
  using parent::n_voxels;
  using parent::n_upper_elements;
  using parent::n_upper;
  using parent::n_lower;
  using parent::internal_name;
  using parent::internal_init;
  
  using parent::influence_matrix;
  using parent::tau_species_single_scattering;
  using parent::singlescat; 
  using parent::sourcefn;

  using parent::species_density; 
  using parent::species_density_pt;
  using parent::species_T; 
  using parent::species_T_pt;

  using parent::absorber_density; 
  using parent::absorber_density_pt; 
  using parent::tau_absorber_single_scattering; 

  //Radiative transfer parameters
  Real solar_lyman_beta_brightness_Hz; /* ph / cm2 / s / Hz <--!
					  solar lyman beta brightness
					  pumping the J=2 lower state */
  // For lyman beta, F_lyb ~= F_lya/66
  //                       ~~ 1/66*(2.5-6.0 x 10^12 ph / cm2 / s / nm)
  //                       ~~ 3.8-9.1 x 10^10 ph / cm2 / s/ nm
  // conversion between 1/nm and 1/Hz is lambda^2 / c = 3.51e-14 nm / Hz
  //
  // F_lyb ~= 1.3-3.2 x 10^-3 ph / cm2 / s / Hz  
  
public:
  void set_solar_brightness(const Real &g) {
    solar_lyman_beta_brightness_Hz = g;
  };
  
  Real get_solar_brightness() const {
    return solar_lyman_beta_brightness_Hz;
  };
  
  // compute the single scattering from the input tracker
  CUDA_CALLABLE_MEMBER
  void compute_single_scattering(const int &start_voxel, influence_tracker &tracker, bool sun_visible = true) {
    for (int i_line=0;i_line<n_lines;i_line++) {
      const int i_lower = tracker.lower_level_index(i_line);
      const int i_upper = tracker.upper_level_index(i_line);
      //const int i_multiplet = tracker.multiplet_index(i_line);

      if (!sun_visible) {
	// single scattering point is behind limb, populate the tracker with the appropriate values
	tracker.tau_species_final[i_line]  = -1.0;
	tracker.tau_absorber_final[i_line] = -1.0;
	tracker.holstein_T_final[i_line]   = 0.0;
      }

      tau_species_single_scattering(start_voxel, i_line) = tracker.tau_species_final[i_line];//line center optical depth
      assert(!std::isnan(tau_species_single_scattering(start_voxel, i_line))
	     && (tau_species_single_scattering(start_voxel, i_line) >= 0
		 || tau_species_single_scattering(start_voxel, i_line) == -1)
	     && "optical depth must be real and positive, or -1 if point is behind limb");
	
      tau_absorber_single_scattering(start_voxel, i_line) = tracker.tau_absorber_final[i_line];
      assert(!std::isnan(tau_absorber_single_scattering(start_voxel, i_line))
	     && (tau_absorber_single_scattering(start_voxel, i_line) >= 0
		 ||  tau_absorber_single_scattering(start_voxel, i_line) == -1)
	     && "optical depth must be real and positive, or -1 if point is behind limb");
	
      if (tracker.lower_level_J(i_line) == 2) {
	// assuming the only excitation of the multiplet is from the
	// J=2 state, whose lines overlap Lyman beta

	// one line connects the J=2 state to each upper state.
	
	// this puts some excitation into every upper state, so we
	// don't need to worry about non-initialized values
	Real solar_line_excitation = (solar_lyman_beta_brightness_Hz // ph / cm2 / s / Hz
				      * species_density(start_voxel, i_lower) // cm-3
				      * tracker.line_sigma_total(i_line) // cm2 Hz
				      / tracker.upper_state_decay_rate(i_upper) // 1/(ph/s)
				      ); // with solar flux included all units except density cancel, we are left with cm-3
      
	singlescat(start_voxel, i_upper) = solar_line_excitation*tracker.holstein_T_final[i_line]; // cm-3, single scattering upper state density
	assert(!std::isnan(singlescat(start_voxel, i_upper))
	       && singlescat(start_voxel, i_upper) >= 0
	       && "single scattering coefficient must be real and positive");
      }
    }
  }

public:

  // definition of class members needed for RT
  template<typename C>
  void define(string emission_name,
	      const C &atmosphere,
	      void (boost::type_identity<C>::type::*species_density_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
	      void (boost::type_identity<C>::type::*species_T_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
	      void (boost::type_identity<C>::type::*absorber_density_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
	      const atmo_voxel (&voxels)[N_VOXELS]) {

    strcpy(internal_name, emission_name.c_str());
    
    for (unsigned int i_voxel=0;i_voxel<N_VOXELS;i_voxel++) {

      (atmosphere.*species_T_function)(voxels[i_voxel],
				       species_T(i_voxel),
				       species_T_pt(i_voxel));
      assert(!std::isnan(species_T(i_voxel))
	     && species_T(i_voxel) >= 0
	     && "temperatures must be real and positive");
      assert(!std::isnan(species_T_pt(i_voxel))
	     && species_T_pt(i_voxel) >= 0
	     && "temperatures must be real and positive");
      
      (atmosphere.*absorber_density_function)(voxels[i_voxel],
					      absorber_density(i_voxel),
					      absorber_density_pt(i_voxel));
      assert(!std::isnan(absorber_density(i_voxel))
	     && absorber_density(i_voxel) >= 0
	     && "densities must be real and positive");
      assert(!std::isnan(absorber_density_pt(i_voxel))
	     && absorber_density_pt(i_voxel) >= 0
	     && "densities must be real and positive");

      // we assume the lower states are collisionally populated, with
      // no RT contributions (can check this after solution by
      // comparing the average radiatiave redistribution rates based
      // on the calculated upper state density with the expected
      // collisional rates)

      Real bulk_density;
      Real bulk_density_pt;
      (atmosphere.*species_density_function)(voxels[i_voxel],
					     bulk_density,
					     bulk_density_pt);
      assert(!std::isnan(bulk_density)
	     && bulk_density >= 0
	     && "densities must be real and positive");
      assert(!std::isnan(bulk_density_pt)
	     && bulk_density_pt >= 0
	     && "densities must be real and positive");

      // now compute the statistical weight of each J level
      Real J_state_fraction[n_lower];
      Real total_statistical_weight = 0;
      Real J_state_fraction_pt[n_lower];
      Real total_statistical_weight_pt = 0;
      brightness_tracker temp_tracker;
      for (int i_lower=0; i_lower<n_lower; i_lower++) {
	J_state_fraction[i_lower] = (temp_tracker.lower_state_statistical_weight(i_lower)
				     * exp(-temp_tracker.lower_state_energy(i_lower)/kB/species_T(i_voxel)));
	total_statistical_weight += J_state_fraction[i_lower];
				     
	J_state_fraction_pt[i_lower] = (temp_tracker.lower_state_statistical_weight(i_lower)
					* exp(-temp_tracker.lower_state_energy(i_lower)/kB/species_T_pt(i_voxel)));
	total_statistical_weight_pt += J_state_fraction_pt[i_lower];

      }
      
      // divide by the computed total weight to get the population
      // fraction, and assign this to the species density
      for (int i_lower=0; i_lower<n_lower; i_lower++) {
	J_state_fraction[i_lower] /= total_statistical_weight;
	species_density(i_voxel, i_lower) = J_state_fraction[i_lower] * bulk_density;

	J_state_fraction_pt[i_lower] /= total_statistical_weight_pt;
	species_density_pt(i_voxel, i_lower) = J_state_fraction_pt[i_lower] * bulk_density_pt;
      }
    }

    parent::reset_solution();

    internal_init=true;
  }

#ifdef __CUDACC__
  using parent::device_emission;
  using parent::vector_to_device;
  using parent::copy_trivial_member_to_device;

  void copy_trivial_members_to_device() {
    parent::copy_trivial_members_to_device();
    copy_trivial_member_to_device(solar_lyman_beta_brightness_Hz, device_emission->solar_lyman_beta_brightness_Hz);
  }

  using parent::copy_trivial_member_to_host;
  void copy_trivial_members_to_host() {
    parent::copy_trivial_members_to_host();
    copy_trivial_member_to_host(solar_lyman_beta_brightness_Hz, device_emission->solar_lyman_beta_brightness_Hz);
  }

  using parent::copy_to_device_influence;
  using parent::copy_to_device_brightness;
  using parent::copy_solved_to_host;
  using parent::device_clear;
#endif

};

#endif
