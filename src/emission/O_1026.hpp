//O_1026.hpp --- routines to compute O 102.6 nm emission

#ifndef __O_1026_h
#define __O_1026_h

#include <boost/type_traits/type_identity.hpp> //for type deduction in define
#include "emission_voxels.hpp"
#include "atmo_vec.hpp"
#include "O_1026_tracker.hpp"

template <int N_VOXELS>
struct O_1026_emission : emission_voxels<N_VOXELS,
					 /*emission_type = */ O_1026_emission<N_VOXELS>,
					 /*los_tracker_type = */ O_1026_tracker> {
protected:
  typedef emission_voxels<N_VOXELS,
			  3,
			  O_1026_emission<N_VOXELS>,
			  O_1026_tracker> parent;
  friend parent;

  // wavelength info and line shape routines are stored in tracker
  // object b/c this compiles to code that is 30% faster, not sure why
  static const int n_lines      = O_1026_tracker<true, N_VOXELS>::n_lines; // = 6
  static const int n_multiplets = O_1026_tracker<true, N_VOXELS>::n_multiplets; // = 3
  static const int n_lambda     = O_1026_tracker<true, N_VOXELS>::n_lambda; // ~= 20

public:
  template <bool transmission, int NV>
  using los = O_1026_tracker<transmission, NV>;
  using typename parent::brightness_tracker;
  using typename parent::influence_tracker;

  
  ~O_1026_emission() {
#if defined(__CUDACC__) and not defined(__CUDA_ARCH__)
    device_clear();
#endif
  }

protected:
  using parent::n_elements;
  using parent::n_voxels;
  using parent::n_upper;
  using parent::n_lower;
  using parent::internal_name;
  using parent::internal_init;
  using parent::internal_solved;

  //Radiative transfer parameters
  Real solar_lyman_beta_brightness_Hz; /* ph / cm2 / s / Hz <--!
					  solar lyman beta brightness
					  pumping the J=2 lower state */
  
  //parent RT parameters
  using vv_1 = typename parent::vv_1;
  using vv_line = typename parent::vv_line;
  using vv_lower = typename parent::vv_lower;
  using vv_upper = typename parent::vv_upper;

  using parent::influence_matrix;
  using parent::tau_species_single_scattering;
  using parent::tau_absorber_single_scattering;
  using parent::singlescat; 
  using parent::sourcefn;

  vv_lower species_density; //average and point densities of species lower states on the grid
  vv_lower species_density_pt;
  vv_1 species_T; // temperature of the species on the grid
  vv_1 species_T_pt;

  vv_1 absorber_density; // average and point densities of absorber on the grid 
  vv_1 absorber_density_pt; 

  // main update routine
  template<bool influence>
  CUDA_CALLABLE_MEMBER
  void update_tracker_start(const Real &current_species_T,
			    const Real (&current_species_density)[n_lower],
			    const Real &current_absorber_density,
			    const Real &pathlength,
			    los<influence,n_elements> &tracker) const {

    // update the transfer probabilities across this voxel
    Real lineshape[n_lines][n_lambda];
    Real dtau_lambda_voxel[n_multiplets][n_lambda];
    Real transfer_probability_lambda_voxel[n_multiplets][n_lambda];
    Real transfer_probability_lambda_final[n_multiplets][n_lambda];

    // initialize
    for (int i_multiplet = 0; i_multiplet < n_multiplets; i_multiplet++) {
      for (int i_lambda = 0; i_lambda < n_lambda; i_lambda++) {
	Real dtau_lambda_voxel[n_multiplets][n_lambda] = 0.0;
	Real transfer_probability_lambda_voxel[n_multiplets][n_lambda] = 1.0;
	Real transfer_probability_lambda_final[n_multiplets][n_lambda] = 1.0;
      }
    }

    // compute optical depth and transfer probabilities
    for (int i_line=0; i_line<n_lines; i_line++) {
      const int i_lower = tracker.lower_level[i_line];
      const int i_upper = tracker.upper_level[i_line];
      const int i_multiplet = tracker.multiplet_identity[i_line];
      
      for (int i_lambda = 0; i_lambda < n_lambda; i_lambda++) {
	lineshape[i_line][i_lambda] = tracker.line_shape_function_normalized(i_line,
									     i_lambda,
									     current_species_T);
	assert(!std::isnan(lineshape[i_line][i_lambda]) && lineshape[i_line][i_lambda] > 0 &&
	       "lineshape must be real and positive");
	
	dtau_lambda_voxel[i_multiplet][i_lambda] += (current_species_density[i_lower]
						     *tracker.line_center_cx_coef[i_line]
						     *lineshape[i_line][i_lambda]
						     +
						     current_absorber_density
						     *tracker.co2_xsec);
	assert(!std::isnan(dtau_lambda_voxel[i_multiplet][i_lambda])
	     && dtau_lambda_voxel[i_multiplet][i_lambda]>=0
	     && "optical depths must be positive numbers");
	
	transfer_probability_lambda_voxel[i_multiplet][i_lambda] *= std::exp(-dtau_lambda_voxel[i_multiplet][i_lambda]*pathlength);
	assert(!std::isnan(transfer_probability_lambda_voxel[i_multiplet][i_lambda]) &&
	       transfer_probability_lambda_voxel[i_multiplet][i_lambda] >= 0 &&
	       transfer_probability_lambda_voxel[i_multiplet][i_lambda] <= 1 &&
	       "transfer probability is a probability.");
	
	transfer_probability_lambda_final[i_multiplet][i_lambda] = (tracker.transfer_probability_lambda_initial[i_multiplet][i_lambda]
								    * transfer_probability_lambda_voxel[i_multiplet][i_lambda]);
      }

      // get the line center optical depth due to each individual line
      tracker.tau_species_final[i_line] += (current_species_density[i_lower]
					    *tracker.line_center_cx_coef[0]
					    /tracker.line_shape_normalization(current_species_T))*pathlength;
      assert(!std::isnan(tracker.tau_species_final[i_line])
      	     && tracker.tau_species_final[i_line]>=0
      	     && "optical depths must be positive numbers");

      tau_absorber_voxel[i_line] += (current_absorber_density
				     * tracker.co2_xsec
				     * pathlength);
      assert(!std::isnan(tracker.tau_absorber_final[i_line])
	     && tracker.tau_absorber_final[i_line]>=0
	     && "optical depths must be positive numbers");
    }

    // ^^^
    // don't combine these loops--- we need to know the transfer
    // probability due to scattering / absorption of all lines before
    // we can compute influence coefficients.
    // vvv


    // reset tracker quantities
    for (int i_line=0; i_line<n_lines; i_line++) {
      tracker.holstein_T_int[i_line] = 0;
      tracker.holstein_T_final[i_line] = 0;
    }
    for (int i_upper=0; i_upper<n_upper; i_upper++)
      for (int j_upper=0; j_upper<n_upper; j_upper++)
	tracker.holstein_G_int[i_upper][j_upper] = 0;

    Real holstein_T_int_coef[n_lines][n_lambda];
    Real lineshape_at_origin[n_lines][n_lambda];
    
    // now we can compute influence coefficients
    for (int i_line = 0; i_line < n_lines; i_line++) {
      const int i_lower = tracker.lower_level[i_line];
      const int i_upper = tracker.upper_level[i_line];
      const int i_multiplet = tracker.multiplet_identity[i_line];

      for (int i_lambda = 0; i_lambda < n_lambda; i_lambda++) {
	// emission for holstein T int happens at current voxel, use the
	// normalization there
	Real holcoef = tracker.weight(i_lambda);

	// holstein_T_int represents a frequency averaged emission in this voxel observed at the start location
	
	// therefore, we use the lineshape in this voxel to compute holstein T int
	holstein_T_int_coef[i_line][i_lambda] = (holcoef
						 * lineshape[i_line][i_lambda]
						 * tracker.transfer_probability_lambda_initial[i_multiplet][i_lambda]
						 * (REAL(1.0) - transfer_probability_lambda_voxel[i_multiplet][i_lambda])
						 / (dtau_lambda_voxel[i_multiplet][i_lambda]));
	tracker.holstein_T_int[i_line] += holstein_T_int_coef[i_line][i_lambda];
	assert(!std::isnan(tracker.holstein_T_int[i_line]) && tracker.holstein_T_int[i_line] >= 0 &&
	       (tracker.holstein_T_int[i_line] <= tau_species_voxel[i_line] ||
		std::abs(tracker.holstein_T_int[i_line] - tau_species_voxel[i_line]) < EPS)
	       //  ^^ this allows for small rounding errors
	       && "holstein integral must be between 0 and Delta tau b/c 0<=HolT<=1");

	if (i_lambda == 0 || i_lambda == n_lambda-1) {
	      //check that the last element is not contributing too much to the integral
	  assert(!((holstein_T_int_coef[i_line][i_lambda] > 1e-6)
		   && (holstein_T_int_coef[i_line][i_lambda]/tracker.holstein_T_int[i_line] > 1e-2))
		 && "wings of line contribute too much to transmission. Increase LAMBDA_MAX in singlet_CFR_tracker.");
	}
	  
	
	if (influence) {
	  // for single scattering calculation, we want the frequency
	  // integrated absorption probability in the start voxel.

	  // we need the lineshape at the start voxel
	  lineshape_at_origin[i_line][i_lambda] = tracker.line_shape_function_normalized(i_line,
											 i_lambda,
											 tracker.species_T_at_origin);
	  assert(!std::isnan(lineshape_at_origin[i_line][i_lambda])
		 && lineshape_at_origin[i_line][i_lambda]>0
		 && "lineshape must be real and positive");

	  // holstein T final represents a frequency-averaged absorption and
	  // uses lineshape_at_origin
	  tracker.holstein_T_final[i_line] += (holcoef
					       * lineshape_at_origin[i_line][i_lambda]
					       * transfer_probability_lambda_final[i_multiplet][i_lambda]);
	  assert(!std::isnan(tracker.holstein_T_final[i_line])
		 && tracker.holstein_T_final[i_line] >= 0
		 && tracker.holstein_T_final[i_line] <= 1
		 && "holstein function represents a probability");

	}
      }
    }

    // ^^^
    // don't combine these loops--- we need to know all lineshapes and holstein_T_int_coefs
    // before we can compute line-to-line transfer.
    // vvv

    if (influence) {
      // now that we've computed all the line shapes we need, we can
      // compute the mixing between levels due to each line
      // emission/absorption
      for (int i_line_origin = 0; i_line_origin < n_lines; i_line_origin++) {
	const int i_lower_origin = tracker.lower_level[i_line];
	const int i_upper_origin = tracker.upper_level[i_line];
	const int i_multiplet_origin = tracker.multiplet_identity[i_line];

	for (int j_line_current = 0; j_line_current < n_lines; j_line_current++) {
	  const int j_lower_current = tracker.lower_level[j_line];
	  const int j_upper_current = tracker.upper_level[j_line];
	  const int j_multiplet_current = tracker.multiplet_identity[j_line];

	  if (i_multiplet_origin == j_multiplet_current) {
	    // the emission in the distant voxel can be captured into any of the upper states

            for (int i_lambda = 0; i_lambda < n_lambda; i_lambda++) {
              // the influence coefficient represents emission in the
              // current voxel being absorbed in the start voxel,
              // possibly into a different upper state
              tracker.holstein_G_int[i_upper_origin][j_upper_current] += (tracker.line_sigma_total[i_line_origin]          // cm2 Hz
									  * species_density_at_origin[i_lower_origin]      // cm-3
									  * tracker.line_A[j_line_current]                 // ph/s
									  / tracker.upper_state_decay_rate[i_upper_origin] // 1/(ph/s)
									  * lineshape_at_origin[i_line_origin][i_lambda]   // Hz^-1
									  * holstein_T_int_coef[j_line_current][i_lambda]  // cm
									  ); // unitless influence coefficient
              assert(!std::isnan(tracker.holstein_G_int[i_upper_origin][j_upper_current])
		     && tracker.holstein_G_int[i_upper_origin][j_upper_current] >= 0
		     && tracker.holstein_G_int[i_upper_origin][j_upper_current] <= 1
		     && "holstein G integral represents a probability");
              // note: this coefficient should be added to the influence matrix at
	      //       (row, col) = (start_voxel_state, current_voxel_state),
	      //       because it represents emission from the current voxel
	      //       and absorption into the start voxel.
            }
          }
        }
      }
    }

    for (int i_multiplet = 0; i_multiplet < n_multiplets; i_multiplet++) {
      for (int i_line = 0; i_line < n_lines; i_line++) {
        // update the initial values
        tracker.transfer_probability_lambda_initial[i_multiplet][i_lambda] = transfer_probability_lambda_final[i_multiplet][i_lambda];
      }
    }

    for (int i_line = 0; i_line < n_lines; i_line++) {
      // if holstein T is larger than physically possible due to rounding
      // errors, reduce it to the physical limit
      if (tracker.holstein_T_int[i_line] > tau_species_voxel[i_line])
        tracker.holstein_T_int[i_line] = tau_species_voxel[i_line]; /*we don't need to worry
								      about rounding errors
								      here because we checked
								      earlier*/
    }
  }

  CUDA_CALLABLE_MEMBER
  void update_tracker_brightness(const Real (&sourcefn_temp)[n_upper], brightness_tracker &tracker) const {
    // called by parent method that specifies interp or nointerp

    for (int i_line=0;i_line<n_lines;i_line++) {
      const int i_upper = tracker.upper_level[i_line];
      
      //bishop formulation
      tracker.brightness[i_line] += (sourcefn_temp[i_upper] // cm-3
				     * tracker.line_A[i_line] // ph/s-1
				     * tracker.holstein_T_int[i_line] // cm
				     / REAL(1e9)) // converts to kR, 10^9 ph/cm2/s, see C&H pg 280-282
      assert(!isnan(tracker.brightness[i_line]) && "brightness must be a real number");
      assert(tracker.brightness[i_line]>=0 && "brightness must be positive");
    }
  }
  
public:

  void set_solar_brightness(const Real &g) {
    solar_lyman_beta_brightness_Hz = g;
  };

  Real get_solar_brightness() const {
    return solar_lyman_beta_brightness_Hz;
  };
  
  //overloads of RT methods
  template<bool influence>
  CUDA_CALLABLE_MEMBER
  void reset_tracker(const int &start_voxel,
		     los<influence,n_elements> &tracker) const {
    Real density_at_origin[n_lower];
    if (!influence)
      // start voxel values don't matter for brightness calculations
      tracker.reset(0.0, density_at_origin);
    else {
      for (int i_lower=0;i_lower<n_lower;i_lower++)
	density_at_origin[i_lower] = species_density(start_voxel, i_lower);
      tracker.reset(species_T[start_voxel], density_at_origin);
    }
  }
  
  template<bool influence>
  CUDA_CALLABLE_MEMBER
  void update_tracker_start(const int &current_voxel,
			    const Real & pathlength,
			    los<influence,n_elements> &tracker) const {
    Real current_species_density[n_lower];
    for (int i_lower=0;i_lower<n_lower;i_lower++)
      current_species_density[i_lower] = species_density(current_voxel, i_lower);

    update_tracker_start(species_T(current_voxel),
			 current_species_density,
			 absorber_density(current_voxel),
			 pathlength,
			 tracker);
  }

  template<bool influence>  
  CUDA_CALLABLE_MEMBER
  void update_tracker_start_interp(const int &n_interp_points,
				   const int *indices,
				   const Real *weights,
				   const Real &pathlength,
				   los<influence,n_elements> &tracker) const {


    Real species_T_interp[1];
    parent::interp_voxel_vector(n_interp_points, indices, weights, species_T_pt, species_T_interp);

    Real species_density_interp[n_lower];
    parent::interp_voxel_vector(n_interp_points, indices, weights, species_density_pt, species_density_interp);

    Real absorber_density_interp[1];
    parent::interp_voxel_vector(n_interp_points, indices, weights, absorber_density_pt, absorber_density_interp);

    
    update_tracker_start(species_T_interp[0],
			 species_density_interp,
			 absorber_density_interp[0],
			 pathlength,
			 tracker);

  }

  using parent::update_tracker_end;

  // update the influence tracker with the contribution from this voxel
  CUDA_CALLABLE_MEMBER
  void update_tracker_influence(const int &current_voxel,
				const Real &pathlength,
				const Real &domega,
				influence_tracker &tracker) const {
    //update influence functions for this voxel
    update_tracker_start(current_voxel, pathlength, tracker);

    for (int i_upper_origin=0; i_upper_origin<n_upper; i_upper_origin++) {
      for (int j_upper_current=0; j_upper_current<n_upper; j_upper_current++) {
	
	//see Bishop1999 for derivation of this formula
	Real coef = domega;
	coef *= tracker.holstein_G_int[i_upper_origin][j_upper_current];
	
	assert(!isnan(coef) && "influence coefficients must be real numbers");
	assert(0<=coef && coef<=1 && "influence coefficients represent transition probabilities");
	
	tracker.influence[i_upper_origin](current_voxel, j_upper_current) += coef;
      }
    }
    update_tracker_end(tracker);
  }

  // compute the single scattering from the input tracker
  CUDA_CALLABLE_MEMBER
  void compute_single_scattering(const int &start_voxel, influence_tracker &tracker, bool sun_visible = true) {
    for (int i_line=0;i_line<n_lines;i_line++) {
      const int i_lower = tracker.lower_level[i_line];
      const int i_upper = tracker.upper_level[i_line];
      const int i_multiplet = tracker.multiplet_identity[i_line];

      if (!sun_visible) {
	// single scattering point is behind limb, populate the tracker with the appropriate values
	tracker.tau_species_final[i_line]  = -1.0;
	tracker.tau_absorber_final[i_line] = -1.0;
	tracker.holstein_T_final[i_line]   = 0.0;
      }

      tau_species_single_scattering(start_voxel, i_line) = tracker.tau_species_final[i_line];//line center optical depth
      assert(!isnan(tau_species_single_scattering(start_voxel, i_line))
	     && (tau_species_single_scattering(start_voxel, i_line) >= 0
		 || tau_species_single_scattering(start_voxel, i_line) == -1)
	     && "optical depth must be real and positive, or -1 if point is behind limb");
	
      tau_absorber_single_scattering(start_voxel, i_line) = tracker.tau_absorber_final[i_line];
      assert(!isnan(tau_absorber_single_scattering(start_voxel, i_line))
	     && (tau_absorber_single_scattering(start_voxel, i_line) >= 0
		 ||  tau_absorber_single_scattering(start_voxel, i_line) == -1)
	     && "optical depth must be real and positive, or -1 if point is behind limb");
	
      if (i_lower == 2) {
	// assuming the only excitation of the multiplet is from the
	// J=2 state, whose lines overlap Lyman beta

	// one line connects the J=2 state to each upper state.
	
	// this puts some excitation into every upper state, so we
	// don't need to worry about non-initialized values
	Real solar_line_excitation = (solar_lyman_beta_brightness_Hz // ph / cm2 / s / Hz
				      * species_density(i_voxel, i_lower) // cm-3
				      * tracker.line_sigma_total[i_line] // cm2 Hz
				      / tracker.upper_state_decay_rate[i_upper] // 1/(ph/s)
				      ); // all units except density cancel, we are left with cm-3
      
	singlescat(start_voxel, i_upper) = solar_line_excitation*tracker.holstein_T_final[i_line]; // cm-3, single scattering upper state density
	assert(!isnan(singlescat(start_voxel, i_upper))
	       && singlescat(start_voxel, i_upper) >= 0
	       && "single scattering coefficient must be real and positive");
      }
    }
  }

  using parent::accumulate_influence;

  void pre_solve() { }
  void pre_solve_gpu(); //defined below
protected:
#ifdef __CUDACC__
  template <int NV>
  friend __global__ void O_1026_prepare_for_solution(O_1026_emission<NV> *emission);
#endif
public:

  using parent::solve;
  using parent::update_tracker_brightness_interp;
  using parent::update_tracker_brightness_nointerp;
  
  // definition of class members needed for RT
  template<typename C>
  void define(const C &atmosphere,
	      void (boost::type_identity<C>::type::*species_density_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
	      void (boost::type_identity<C>::type::*species_T_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
	      void (boost::type_identity<C>::type::*absorber_density_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
	      const atmo_voxel (&voxels)[N_VOXELS]) {

    strcpy(internal_name, emission_name.c_str());
    branching_ratio     = emission_branching_ratio;
    species_T_ref       = species_T_reff;
    species_sigma_T_ref = species_sigma_T_reff;
    
    for (unsigned int i_voxel=0;i_voxel<N_VOXELS;i_voxel++) {

      (atmosphere.*species_T_function)(voxels[i_voxel],
				       species_T(i_voxel),
				       species_T_pt(i_voxel));
      assert(!isnan(species_T(i_voxel))
	     && species_T(i_voxel) >= 0
	     && "temperatures must be real and positive");
      assert(!isnan(species_T_pt(i_voxel))
	     && species_T_pt(i_voxel) >= 0
	     && "temperatures must be real and positive");
      
      (atmosphere.*absorber_density_function)(voxels[i_voxel],
					      absorber_density(i_voxel),
					      absorber_density_pt(i_voxel));
      assert(!isnan(absorber_density(i_voxel))
	     && absorber_density(i_voxel) >= 0
	     && "densities must be real and positive");
      assert(!isnan(absorber_density_pt(i_voxel))
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
					     bulk_density_pt));
      assert(!isnan(bulk_density)
	     && bulk_density >= 0
	     && "densities must be real and positive");
      assert(!isnan(bulk_density_pt)
	     && bulk_density_pt >= 0
	     && "densities must be real and positive");

      // now compute the statistical weight of each J level
      Real J_state_fraction[n_lower];
      Real total_statistical_weight = 0;
      Real J_state_fraction_pt[n_lower];
      Real total_statistical_weight_pt = 0;
      for (int i_lower=0; i_lower<n_lower; i_lower++) {
	J_state_fraction[i_lower] = (brightness_tracker::lower_state_statistical_weight[i_lower]
				     * exp(-brightness_tracker::lower_state_energy[i_lower]/kB/species_T(i_voxel)));
	total_statistical_weight += J_state_fraction[i_lower];
				     
	J_state_fraction_pt[i_lower] = (brightness_tracker::lower_state_statistical_weight[i_lower]
				     * exp(-brightness_tracker::lower_state_energy[i_lower]/kB/species_T_pt(i_voxel));
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

  template <int N_STATES>
  void save_voxel_state(std::ostream &file, VectorX (*function)(VectorX, int), const int i,
			voxel_vector<n_voxels, N_STATES> vector, string message) {
    VectorX voxel_quantity;
    voxel_quantity.resize(n_voxels);
    for (int i_state=0; i_state<N_STATES; i_state++) {
      for (int i_voxel=0; i_voxel<n_voxels; i_voxel++)
	voxel_quantity[i_voxel] = vector(i_voxel, i_state);

      file << "      " << message << " " << i_state << ": " << "\n"
	   << function(voxel_quantity, i).transpose() << "\n";
    }
  }
  
  
  void save(std::ostream &file, VectorX (*function)(VectorX, int), const int i) const {

    file << "    Species density [cm-3]: " << "\n";
    save_voxel_state(file, function, i,
		     species_density, "lower state", n_lower);
      
    file << "    Species single scattering tau: " << "\n";
    save_voxel_state(file, function, i,
		     tau_species_single_scattering, "line", n_lines);
    
    file << "    Absorber density [cm-3]: " 
	 <<      function(absorber_density.eigen(), i).transpose() << "\n"
      
	 << "    Absorber single scattering tau: " 
	 <<	 function(tau_absorber_single_scattering.eigen(), i).transpose() << "\n";
      
    file << "    Species single scattering source function S0: " 
    save_voxel_state(file, function, i,
		     singlescat, "upper state", n_upper);
      
    file << "    Source function: " 
    save_voxel_state(file, function, i,
		     sourcefn, "upper state", n_upper);
  }

  void save_brightness(std::ostream &file, const gpu_vector<brightness_tracker> &los_brightness) const {
    // save a list of brightness values from brightness trackers to file

    MatrixX brightness_write_out;
    brightness_write_out.resize(los_brightness.size(), n_lines);
    
    for (int i_obs=0;i_obs<los_brightness.size();i_obs++)
      for (int i_line=0;i_line<n_lines;i_line++)
	brightness_write_out(i_obs, i_line) = los_brightness[i].brightness[i_line];

    file << internal_name
	 << " brightness [kR]: \n";
    for (int i_line=0;i_line<n_lines;i_line++)
      file << "   line " << i_line << ": ";
      for (int i_obs=0;i_obs<los_brightness.size();i_obs++)
	file << brightness_write_out(i_obs,i_line);
      file << "\n";
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

  void copy_to_device_influence() {
    copy_trivial_members_to_device();
    parent::copy_to_device_influence();
  
    //copy the read-only atmosphere arrays
    bool transfer = true;
    vector_to_device(device_emission->species_density, species_density, transfer);
    vector_to_device(device_emission->species_T, species_T, transfer);
    vector_to_device(device_emission->absorber_density, absorber_density, transfer);
  }

  void copy_to_device_brightness()
  {
    copy_trivial_members_to_device();
    parent::copy_to_device_brightness();

    //free some of the influnce stuff if we used it
    species_density.free_d_vec();
    species_T.free_d_vec();
    absorber_density.free_d_vec();

    //copy the stuff we need to do the calculation on the device
    bool transfer = true;
    vector_to_device(device_emission->species_density_pt, species_density_pt, transfer);
    vector_to_device(device_emission->species_T_pt, species_T_pt, transfer);
    vector_to_device(device_emission->absorber_density_pt, absorber_density_pt, transfer);
  }
  
  void device_clear() {
    species_density.free_d_vec();
    species_T.free_d_vec();
    absorber_density.free_d_vec();

    species_density_pt.free_d_vec();
    species_T_pt.free_d_vec();
    absorber_density_pt.free_d_vec();

    parent::device_clear();
  }
#endif

};

#ifdef __CUDACC__

template<int N_VOXELS>
__global__
void O_1026_prepare_for_solution(O_1026_emission<N_VOXELS> *emission)
//TODO: move this to parent class?
{
  //each block prepares one row of the influence matrix
  int i_vox = blockIdx.x;
  
  //convert to kernel from influence (kernel = identity - branching_ratio*influence)
  for (int j_vox = threadIdx.x; j_vox < emission->n_elements; j_vox+=blockDim.x)
    emission->influence_matrix(i_vox,j_vox) *= -1.0;

  if (threadIdx.x == 0) {
    //add the identity matrix
    emission->influence_matrix(i_vox,i_vox) += 1;
    
    //now copy singlescat to sourcefn, preparing for in-place solution
    emission->sourcefn(i_vox) = emission->singlescat(i_vox);
  }
}


template <int N_VOXELS>
void O_1026_emission<N_VOXELS>::pre_solve_gpu() {
  const int n_threads = 32;
  O_1026_prepare_for_solution<<<n_elements, n_threads>>>(device_emission);
  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );
}

#endif


#endif
