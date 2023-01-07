//H_lyman_multiplet_test.hpp --- check multiplet and singlet codes for
//consistency

// most of the code to compute multiplet emission is generic, and
// stored in multiplet_CFR_emission.hpp

// the routines to define the emission densities / temperatures are
// specific to H Lyman alpha and stored here, so are the solar excitation
// calculations


#ifndef __H_lyman_multiplet_test
#define __H_lyman_multiplet_test

#include <boost/type_traits/type_identity.hpp> //for type deduction in define
#include "emission_voxels.hpp"
#include "atmo_vec.hpp"
#include "multiplet_CFR_emission.hpp"
#include "H_multiplet_tracker_test.hpp"

template <int N_VOXELS>
struct H_lyman_singlet : multiplet_CFR_emission<N_VOXELS,
						  /*emission_type = */ H_lyman_singlet<N_VOXELS>,
						  /*los_tracker_type = */ H_lyman_singlet_tracker> {
protected:
  typedef multiplet_CFR_emission<N_VOXELS,
				 /*emission_type = */ H_lyman_singlet<N_VOXELS>,
				 /*los_tracker_type = */ H_lyman_singlet_tracker> parent;
  friend parent;

public:
  // wavelength info and line shape routines are stored in tracker
  // object b/c this compiles to code that is 30% faster, not sure why
  using parent::n_lines;
  using parent::n_multiplets;
  using parent::n_lambda;
  
  template <bool influence>
  using los = H_lyman_singlet_tracker<influence, N_VOXELS>;
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
  Real solar_brightness_lya_Hz; /* ph / cm2 / s / Hz <--!
			           solar lyman alpha / beta brightness
			           pumping the lower state */
  Real solar_brightness_lyb_Hz; /* ph / cm2 / s / Hz <--!
			           solar lyman alpha / beta brightness
			           pumping the lower state */

  // For lyman alpha, F_lya ~= 2.5-6.0 x 10^12 ph / cm2 / s / nm
  // For lyman beta,  F_lyb ~= F_lya/66
  //                        ~= 3.8-9.1 x 10^10 ph / cm2 / s/ nm

  // conversion between 1/nm and 1/Hz is lambda^2 / c = 4.93e-14 nm / Hz for lya
  //                                                  = 3.51e-14 nm / Hz for lyb
  
  // F_lya ~= 1.2-3.0 x 10^-1 ph / cm2 / s / Hz
  // F_lyb ~= 1.3-3.2 x 10^-3 ph / cm2 / s / Hz


  bool constant_temp_RT_internal; // whether to use a constant temperature for RT calculations
  Real constant_temp_RT;     // ^^^^^^^ Temp value to use

  bool CO2_absorption;
  
public:
  void set_solar_brightness(const Real &lya, const Real &lyb) {
    solar_brightness_lya_Hz = lya;
    solar_brightness_lyb_Hz = lyb;
  };
  
  Real get_solar_brightness_lya() const {
    return solar_brightness_lya_Hz;
  };

  Real get_solar_brightness_lyb() const {
    return solar_brightness_lyb_Hz;
  };

  void set_atmosphere_temp_RT() {
    constant_temp_RT_internal = false;
  }
  
  void set_constant_temp_RT(Real temp) {
    constant_temp_RT_internal = true;
    constant_temp_RT = temp;
  }

  void set_CO2_absorption_off() {
    CO2_absorption = false;
  }

  void set_CO2_absorption_on() {
    CO2_absorption = true;
  }

  H_lyman_singlet()
    : constant_temp_RT_internal(false),
      CO2_absorption(true)
  {}
  
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
	
      Real solar_brightness_Hz = i_line == 0 ? solar_brightness_lya_Hz : solar_brightness_lyb_Hz;
      Real solar_line_excitation = (solar_brightness_Hz // ph / cm2 / s / Hz
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

      if (constant_temp_RT_internal) {
	// use a constant temperature specified by the user in an earlier call of set_constant_temp_RT()
	species_T(i_voxel)    = constant_temp_RT;
	species_T_pt(i_voxel) = constant_temp_RT;
      } else {
	// get the temperature to use from the atmosphere object
	(atmosphere.*species_T_function)(voxels[i_voxel],
					 species_T(i_voxel),
					 species_T_pt(i_voxel));
      }
      assert(!std::isnan(species_T(i_voxel))
	     && species_T(i_voxel) >= 0
	     && "temperatures must be real and positive");
      assert(!std::isnan(species_T_pt(i_voxel))
	     && species_T_pt(i_voxel) >= 0
	     && "temperatures must be real and positive");

      if (CO2_absorption) {
	(atmosphere.*absorber_density_function)(voxels[i_voxel],
						absorber_density(i_voxel),
						absorber_density_pt(i_voxel));
      } else {
	absorber_density(i_voxel) = 0.0;
	absorber_density_pt(i_voxel) = 0.0;
      }
      assert(!std::isnan(absorber_density(i_voxel))
	     && absorber_density(i_voxel) >= 0
	     && "densities must be real and positive");
      assert(!std::isnan(absorber_density_pt(i_voxel))
	     && absorber_density_pt(i_voxel) >= 0
	     && "densities must be real and positive");

      (atmosphere.*species_density_function)(voxels[i_voxel],
					     species_density(i_voxel, 0),
					     species_density_pt(i_voxel, 0));
      assert(!std::isnan(species_density(i_voxel, 0))
	     && species_density(i_voxel, 0) >= 0
	     && "densities must be real and positive");
      assert(!std::isnan(species_density_pt(i_voxel, 0))
	     && species_density_pt(i_voxel, 0) >= 0
	     && "densities must be real and positive");

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
    copy_trivial_member_to_device(solar_brightness_lya_Hz, device_emission->solar_brightness_lya_Hz);
    copy_trivial_member_to_device(solar_brightness_lyb_Hz, device_emission->solar_brightness_lyb_Hz);
    copy_trivial_member_to_device(constant_temp_RT_internal, device_emission->constant_temp_RT_internal);
    copy_trivial_member_to_device(constant_temp_RT, device_emission->constant_temp_RT);
    copy_trivial_member_to_device(CO2_absorption, device_emission->CO2_absorption);
  }

  using parent::copy_trivial_member_to_host;
  void copy_trivial_members_to_host() {
    parent::copy_trivial_members_to_host();
    copy_trivial_member_to_host(solar_brightness_lya_Hz, device_emission->solar_brightness_lya_Hz);
    copy_trivial_member_to_host(solar_brightness_lyb_Hz, device_emission->solar_brightness_lyb_Hz);
    copy_trivial_member_to_host(constant_temp_RT_internal, device_emission->constant_temp_RT_internal);
    copy_trivial_member_to_host(constant_temp_RT, device_emission->constant_temp_RT);
    copy_trivial_member_to_host(CO2_absorption, device_emission->CO2_absorption);
  }

  using parent::copy_to_device_influence;
  using parent::copy_to_device_brightness;
  using parent::copy_solved_to_host;
  using parent::device_clear;
#endif

};

#endif
