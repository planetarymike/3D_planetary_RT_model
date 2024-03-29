//multiplet_CFR.hpp --- generic routines to compute multiplet CFR emission

#ifndef __multiplet_CFR_h
#define __multiplet_CFR_h

#include <boost/type_traits/type_identity.hpp> //for type deduction in define
#include "emission_voxels.hpp"
#include "atmo_vec.hpp"

template <int N_VOXELS, typename derived_emission, template<bool,int> class los_tracker_type>
struct multiplet_CFR_emission : emission_voxels<N_VOXELS,
						/*emission_type = */ derived_emission,
						/*los_tracker_type = */ los_tracker_type> {
protected:
  typedef emission_voxels<N_VOXELS,
			  /*emission_type = */ derived_emission,
			  /*los_tracker_type = */ los_tracker_type> parent;
  friend parent;

public:
  // wavelength info and line shape routines are stored in tracker
  // object b/c this compiles to code that is 30% faster, not sure why
  static const int n_lines      = los_tracker_type<true, N_VOXELS>::n_lines; // = 6
  static const int n_multiplets = los_tracker_type<true, N_VOXELS>::n_multiplets; // = 3
  static const int n_lambda     = los_tracker_type<true, N_VOXELS>::n_lambda; // ~= 21

  template <bool influence>
  using los = los_tracker_type<influence, N_VOXELS>;
  using typename parent::brightness_tracker;
  using typename parent::influence_tracker;

  
  ~multiplet_CFR_emission() {
#if defined(__CUDACC__) and not defined(__CUDA_ARCH__)
    device_clear();
#endif
  }

protected:
  using parent::n_voxels;
  using parent::n_upper_elements;
  using parent::n_upper;
  using parent::n_lower;
  using parent::internal_name;
  using parent::internal_init;
  using parent::internal_solved;

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
			    los<influence> &tracker) const {

    // update the transfer probabilities across this voxel in these
    // working arrays, then transfer into tracker
    Real lineshape[n_lines][n_lambda];
    Real tau_lambda_voxel[n_multiplets][n_lambda];
    Real tau_species_voxel[n_lines];
    Real tau_absorber_voxel[n_lines];
    Real transfer_probability_lambda_voxel[n_multiplets][n_lambda];
    Real transfer_probability_lambda_final[n_multiplets][n_lambda];

    // update column density across voxel
    for (int i_lower = 0; i_lower < n_lower; i_lower++) {
      tracker.species_col_dens[i_lower] += current_species_density[i_lower]*pathlength;
      assert_positive(tracker.species_col_dens[i_lower]);
    }

    // initialize
    for (int i_multiplet=0;i_multiplet<n_multiplets;i_multiplet++)
      for (int i_lambda = 0; i_lambda < n_lambda; i_lambda++)
	tau_lambda_voxel[i_multiplet][i_lambda] = 0.0;

    // compute optical depth contribution from each line
    for (int i_line=0; i_line<n_lines; i_line++) {
      const int i_lower = tracker.lower_level_index(i_line);
      //const int i_upper = tracker.upper_level_index(i_line);
      const int i_multiplet = tracker.multiplet_index(i_line);
      
      for (int i_lambda = 0; i_lambda < n_lambda; i_lambda++) {
	lineshape[i_line][i_lambda] = tracker.line_shape_function_normalized(i_line,
									     i_lambda,
									     current_species_T);
	assert_positive(lineshape[i_line][i_lambda]);
	
	tau_lambda_voxel[i_multiplet][i_lambda] += ((current_species_density[i_lower]  // cm-3
						     *tracker.line_sigma_total(i_line) // cm2 Hz
						     *lineshape[i_line][i_lambda]      // Hz-1
						     +
						     current_absorber_density          // cm-3
						     * tracker.absorber_xsec(i_line))   // cm2
						    *pathlength);                      // cm
	assert_positive(tau_lambda_voxel[i_multiplet][i_lambda]);
      }

      // get the line center optical depth due to each individual line
      tau_species_voxel[i_line] = ((current_species_density[i_lower]                               // cm-3
				    *tracker.line_sigma_total(i_line)                              // cm2 Hz
				    *tracker.line_shape_normalization(i_line, current_species_T))  // Hz-1
				   *pathlength);                                                    // cm
      tracker.tau_species_final[i_line] += tau_species_voxel[i_line];
      assert_positive(tracker.tau_species_final[i_line]);

      tau_absorber_voxel[i_line] = (current_absorber_density         // cm-3
				    * tracker.absorber_xsec(i_line)  // cm2
				    * pathlength);                   // cm
      tracker.tau_absorber_final[i_line] += tau_absorber_voxel[i_line];
      assert_positive(tracker.tau_absorber_final[i_line]);
    }

    // ^^^
    // don't combine these loops--- we need to know the optical depths
    // before we can compute the transfer probability
    // vvv

    for (int i_multiplet=0;i_multiplet<n_multiplets;i_multiplet++) {
      for (int i_lambda = 0; i_lambda < n_lambda; i_lambda++) {
    	// because the loop above updated all of the optical depths we can now compute transfer probabilities
	transfer_probability_lambda_voxel[i_multiplet][i_lambda] = std::exp(-tau_lambda_voxel[i_multiplet][i_lambda]);
	assert_probability(transfer_probability_lambda_voxel[i_multiplet][i_lambda]);
	
	transfer_probability_lambda_final[i_multiplet][i_lambda] = (tracker.transfer_probability_lambda_initial[i_multiplet][i_lambda]
								    * transfer_probability_lambda_voxel[i_multiplet][i_lambda]);
      }
    }
    
    // ^^^
    // don't combine these loops--- we need to know the transfer probabilities
    // before we can compute the influence functions
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
      //const int i_lower = tracker.lower_level_index(i_line);
      //const int i_upper = tracker.upper_level_index(i_line);
      const int i_multiplet = tracker.multiplet_index(i_line);

      for (int i_lambda = 0; i_lambda < n_lambda; i_lambda++) {
	// emission for holstein T int happens at current voxel, use the
	// normalization there
	Real holcoef = tracker.weight(i_line, i_lambda); // Hz

	// holstein_T_int represents an effective path length for the
	// optically thick emission

	// integrate emission across the cell at this wavelength
	if (tau_lambda_voxel[i_multiplet][i_lambda] < 1e-3)
	  // this solves a floating point issue.
	  //   when tau << 1, (1-exp(-tau))/tau ~= 1.0 - tau/2
	  //                                   ^^^^ for tau < 1e-3, this expansion is accurate to >6 decimal digits
	  //   could also use (1-exp(-tau))/tau = 1 - tau/2! + tau^2/3! - tau^3/4! + ...
	  holstein_T_int_coef[i_line][i_lambda] = (REAL(1.0)
						   - (REAL(0.5)*tau_lambda_voxel[i_multiplet][i_lambda]));
						   // *(REAL(1.0) - tau_lambda_voxel[i_multiplet][i_lambda]/REAL(3.0)));
	    
	else
	  holstein_T_int_coef[i_line][i_lambda] = ((REAL(1.0) - transfer_probability_lambda_voxel[i_multiplet][i_lambda])
						   / (tau_lambda_voxel[i_multiplet][i_lambda]));

	// the function (1-exp(-tau))/tau is bounded by 0 and 1. Only
	// floating point errors can put it outside this range.
	assert_probability(holstein_T_int_coef[i_line][i_lambda]); // voxel contribution must be between 0 and 1

	// we use the lineshape in this voxel to compute
	// holstein_T_int because it represents emission in the
	// current voxel
	holstein_T_int_coef[i_line][i_lambda] *= (holcoef                                                               // Hz
						  * lineshape[i_line][i_lambda]                                         // Hz-1
						  * tracker.transfer_probability_lambda_initial[i_multiplet][i_lambda]  // unitless
						  * pathlength);                                                        // cm
	tracker.holstein_T_int[i_line] += holstein_T_int_coef[i_line][i_lambda];
	assert_leq(tracker.holstein_T_int[i_line], pathlength); // holstein integral must be between 0 and pathlength b/c 0<=HolT<=1
	
	if (influence) {
	  // for single scattering calculation, we want the frequency
	  // integrated absorption probability in the start voxel.

	  // we need the lineshape at the start voxel
	  lineshape_at_origin[i_line][i_lambda] = tracker.line_shape_function_normalized(i_line,
											 i_lambda,
											 tracker.species_T_at_origin);
	  assert_positive(lineshape_at_origin[i_line][i_lambda]);

	  // holstein T final represents a frequency-averaged
	  // absorption probability in the origin cell and uses
	  // lineshape_at_origin
	  tracker.holstein_T_final[i_line] += (holcoef                                                       // Hz
					       * lineshape_at_origin[i_line][i_lambda]                       // Hz-1
					       * transfer_probability_lambda_final[i_multiplet][i_lambda]);  // unitless
	  assert_probability(tracker.holstein_T_final[i_line]);

	}
      }
      // check that the first and last element are not contributing too much to the holstein T integral
      assert_small_contribution(holstein_T_int_coef[i_line][0         ], tracker.holstein_T_int[i_line]);
      assert_small_contribution(holstein_T_int_coef[i_line][n_lambda-1], tracker.holstein_T_int[i_line]);
      // If either of the above fail, wings of line contribute too much to transmission.
      // Increase LAMBDA_MAX in singlet_CFR_tracker.
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
	const int i_lower_origin = tracker.lower_level_index(i_line_origin);
	const int i_upper_origin = tracker.upper_level_index(i_line_origin);
	const int i_multiplet_origin = tracker.multiplet_index(i_line_origin);

	for (int j_line_current = 0; j_line_current < n_lines; j_line_current++) {
	  //const int j_lower_current = tracker.lower_level_index(j_line_current);
	  const int j_upper_current = tracker.upper_level_index(j_line_current);
	  const int j_multiplet_current = tracker.multiplet_index(j_line_current);
	  
	  if (i_multiplet_origin == j_multiplet_current) {
	    // lines potentially overlap

            for (int i_lambda = 0; i_lambda < n_lambda; i_lambda++) {
              // the holstein G influence coefficient tracks emission
              // in the current voxel being absorbed in the start
              // voxel, possibly into a different upper state.
              tracker.holstein_G_int[i_upper_origin][j_upper_current] += (tracker.line_sigma_total(i_line_origin)              // cm2 Hz
									  * tracker.species_density_at_origin[i_lower_origin]  // cm-3
									  * tracker.line_A(j_line_current)                     // ph/s
									  / tracker.upper_state_decay_rate(i_upper_origin)     // 1/(ph/s)
									  * lineshape_at_origin[i_line_origin][i_lambda]       // Hz^-1
									  * holstein_T_int_coef[j_line_current][i_lambda]      // cm
									  ); // unitless influence coefficient
	      // holstein G integral ALMOST represents a probability --- but it needs to be multipled by differential solid angle.
	      // this is done in update_tracker_influence, after which we can check that it's between 0 and 1.

              // note: this coefficient should be added to the influence matrix at
	      //       (row, col) = (start_voxel_state, current_voxel_state),
	      //       because it represents emission from the current voxel
	      //       and absorption into the start voxel.
            }
          }
        }
      }
    }

    // calculation is done, do some bookkeeping so the tracker is
    // ready to be updated in the next voxel
    for (int i_multiplet = 0; i_multiplet < n_multiplets; i_multiplet++) {
      for (int i_lambda = 0; i_lambda < n_lambda; i_lambda++) {
        // update the initial values
        tracker.transfer_probability_lambda_initial[i_multiplet][i_lambda] = transfer_probability_lambda_final[i_multiplet][i_lambda];
      }
    }

    for (int i_line = 0; i_line < n_lines; i_line++) {
      // if holstein T int is larger than physically possible due to rounding
      // errors, reduce it to the physical limit
      if (tracker.holstein_T_int[i_line] > pathlength)
        tracker.holstein_T_int[i_line] = pathlength; /*we don't need to worry
						       about rounding errors
						       here because we checked
						       earlier*/
    }
  }

  CUDA_CALLABLE_MEMBER
  void update_tracker_brightness(const Real (&sourcefn_temp)[n_upper], brightness_tracker &tracker) const {
    // called by parent method that specifies interp or nointerp

    for (int i_line=0;i_line<n_lines;i_line++) {
      const int i_upper = tracker.upper_level_index(i_line);
      
      tracker.brightness[i_line] += (sourcefn_temp[i_upper] // cm-3
				     * tracker.line_A(i_line) // ph/s
				     * tracker.holstein_T_int[i_line] // cm
				     / REAL(1e9)); // converts to kR, 10^9 ph/cm2/s, see C&H pg 280-282

      assert_positive(tracker.brightness[i_line]);
    }
  }
  
public:

  //overloads of RT methods
  template<bool influence>
  CUDA_CALLABLE_MEMBER
  void reset_tracker(const int &start_voxel,
		     los<influence> &tracker) const {
    Real density_at_origin[n_lower];
    if (!influence)
      // start voxel values don't matter for brightness calculations
      tracker.reset(0.0, density_at_origin);
    else {
      for (int i_lower=0;i_lower<n_lower;i_lower++)
	density_at_origin[i_lower] = species_density(start_voxel, i_lower);
      tracker.reset(species_T(start_voxel), density_at_origin);
    }
  }
  
  template<bool influence>
  CUDA_CALLABLE_MEMBER
  void update_tracker_start(const int &current_voxel,
			    const Real & pathlength,
			    los<influence> &tracker) const {
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
				   los<influence> &tracker) const {


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

	// influence coefficients represent transition probabilities	
	assert_probability(coef);
	
	tracker.influence[i_upper_origin](current_voxel, j_upper_current) += coef;
      }
    }
    update_tracker_end(tracker);
  }

  using parent::accumulate_influence;

  void pre_solve() { }
  void pre_solve_gpu(); //defined below
protected:
#ifdef __CUDACC__
  template <int NV, typename d_e, template<bool,int> class l_t_t> // need different names here to avoid override
  friend __global__ void multiplet_CFR_prepare_for_solution(multiplet_CFR_emission<NV, d_e, l_t_t> *emission);
#endif
public:

  using parent::solve;
  using parent::update_tracker_brightness_interp;
  using parent::update_tracker_brightness_nointerp;
  
  template <int N_STATES>
  void save_voxel_state(std::ostream &file, VectorX (*function)(VectorX, int), const int i,
			const voxel_vector<n_voxels, N_STATES> &vector, const string message) const {
    VectorX voxel_quantity;
    voxel_quantity.resize(n_voxels);
    for (int i_state=0; i_state<N_STATES; i_state++) {
      for (int i_voxel=0; i_voxel<n_voxels; i_voxel++)
	voxel_quantity[i_voxel] = vector(i_voxel, i_state);

      file << "      " << message << " " << i_state << ": "
	   << function(voxel_quantity, i).transpose() << "\n";
    }
  }
  
  
  void save(std::ostream &file, VectorX (*function)(VectorX, int), const int i) const {

    file << "    Species density [cm-3]: \n";
    save_voxel_state(file, function, i,
		     species_density, "lower state");

    file << "    Temperature [K]: " 
	 <<      function(species_T.eigen(), i).transpose() << "\n";
      
    file << "    Species single scattering tau: \n";
    save_voxel_state(file, function, i,
		     tau_species_single_scattering, "line");
    
    file << "    Absorber density [cm-3]: " 
	 <<      function(absorber_density.eigen(), i).transpose() << "\n"
      
	 << "    Absorber single scattering tau: \n";
    save_voxel_state(file, function, i,
		     tau_absorber_single_scattering, "line");
      
    file << "    Species single scattering source function S0: \n";
    save_voxel_state(file, function, i,
		     singlescat, "upper state");
      
    file << "    Source function: \n";
    save_voxel_state(file, function, i,
		     sourcefn, "upper state");
  }

  void save_brightness(std::ostream &file, const gpu_vector<brightness_tracker> &los_brightness) const {
    // save a list of brightness values from brightness trackers to file
    file << internal_name
	 << " brightness [kR]: \n";
    for (int i_line=0;i_line<n_lines;i_line++) {
      file << "   line " << i_line << ": ";
      for (int i_obs=0;i_obs<los_brightness.size();i_obs++)
	file << los_brightness[i_obs].brightness[i_line] << " ";
      file << "\n";
    }
  }

#ifdef __CUDACC__
  using parent::device_emission;
  using parent::vector_to_device;
  using parent::copy_trivial_member_to_device;
  using parent::copy_trivial_members_to_device;
  using parent::copy_trivial_members_to_host;

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

  void copy_solved_to_host() {
    parent::copy_solved_to_host();
    
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
template <int N_VOXELS, typename derived_emission, template<bool,int> class los_tracker_type>
__global__
void multiplet_CFR_prepare_for_solution(multiplet_CFR_emission<N_VOXELS, derived_emission, los_tracker_type> *emission)
//TODO: move this to parent class?
{
  //each block prepares one row of the influence matrix
  int i_el = blockIdx.x;
  
  //convert to kernel from influence (kernel = identity - branching_ratio*influence)
  for (int j_el = threadIdx.x; j_el < emission->influence_matrix.n_elements; j_el+=blockDim.x)
    emission->influence_matrix(i_el,j_el) *= REAL(-1.0);

  if (threadIdx.x == 0) {
    //add the identity matrix
    emission->influence_matrix(i_el,i_el) += REAL(1.0);
    
    //now copy singlescat to sourcefn, preparing for in-place solution
    emission->sourcefn(i_el) = emission->singlescat(i_el);
  }
}


template <int N_VOXELS, typename derived_emission, template<bool,int> class los_tracker_type>
void multiplet_CFR_emission<N_VOXELS, derived_emission, los_tracker_type>::pre_solve_gpu() {
  const int n_threads = 32;
  multiplet_CFR_prepare_for_solution<<<influence_matrix.n_elements, n_threads>>>(device_emission);
  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );
}

#endif


#endif
