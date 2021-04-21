//singlet_CFR.hpp --- routines to compute singlet emissions with complete frequency redistribution

#ifndef __singlet_CFR_h
#define __singlet_CFR_h

#include <boost/type_traits/type_identity.hpp> //for type deduction in define
#include "emission_voxels.hpp"
#include "atmo_vec.hpp"
#include "los_tracker.hpp"

template <int N_VOXELS>
struct singlet_CFR : emission_voxels<N_VOXELS,
					/*N_STATES_PER_VOXEL = */ 1, // this class tracks only singlet emissions
					/*emission_type = */ singlet_CFR<N_VOXELS>,
					/*los_tracker_type = */ singlet_CFR_tracker> {
protected:
  typedef emission_voxels<N_VOXELS,
			  1,
			  singlet_CFR<N_VOXELS>,
			  singlet_CFR_tracker> parent;
  friend parent;

  // wavelength info and line shape routines are stored in tracker
  // object b/c this compiles to code that is 30% faster
  static const int n_lambda = singlet_CFR_tracker<true, N_VOXELS>::n_lambda;

public:
  template <bool transmission, int NV>
  using los = singlet_CFR_tracker<transmission,NV>;
  using typename parent::brightness_tracker;
  using typename parent::influence_tracker;

  using parent::n_elements;
  using parent::n_voxels;
  using parent::n_states;
  
  ~singlet_CFR() {
#if defined(__CUDACC__) and not defined(__CUDA_ARCH__)
    device_clear();
#endif
  }

protected:
  using parent::internal_name;
  using parent::internal_init;
  using parent::internal_solved;

  //Radiative transfer parameters
  Real branching_ratio;

  Real species_T_ref;//reference temp, K (anything near mean works OK as long
                     //as variations are less than a factor of 5ish)
  Real species_sigma_T_ref;//species cross section at this temperature

  Real emission_g_factor; // g-factor for this emission
  
  //parent RT parameters
  using vv = typename parent::vv;
  using vm = typename parent::vm;
  using parent::influence_matrix;
  using parent::tau_species_single_scattering;
  using parent::tau_absorber_single_scattering;
  using parent::singlescat; 
  using parent::sourcefn;

  vv species_density; //average and point densities of species on the grid
  vv species_density_pt;
  vv species_T_ratio; //average and point T_ref/T of species on the grid
  vv species_T_ratio_pt;
  vv dtau_species; // species density * species_sigma_T_ref * sqrt(species_T_ref/species_T) 
  vv dtau_species_pt;

  vv absorber_density; //same as above for absorber
  vv absorber_density_pt; 
  vv absorber_sigma; //pure absorption does not need T info; can be implicitly T-dependent
  vv absorber_sigma_pt;
  vv dtau_absorber; // absorber_density * absorber_sigma for each voxel independently
  vv dtau_absorber_pt;

  vv abs; //ratio (dtau_absorber / dtau_species)
  vv abs_pt;

  // main update routine
  template<bool influence>
  CUDA_CALLABLE_MEMBER
  void update_tracker_start(const Real &current_species_T_ratio,
			    const Real &current_dtau_species,
			    const Real &current_dtau_absorber,
			    const Real &current_abs,
			    const Real &pathlength,
			    los<influence,n_elements> &tracker) const {
    Real tau_species_voxel = current_dtau_species * pathlength;
    tracker.tau_species_final += tau_species_voxel;
    assert(!std::isnan(tracker.tau_species_final)
	   && tracker.tau_species_final>=0
	   && "optical depths must be real numbers");

    tracker.tau_absorber_final += current_dtau_absorber * pathlength; 
    assert(!std::isnan(tracker.tau_absorber_final)
	   && tracker.tau_absorber_final>=0
	   && "optical depths must be real numbers");
    
    tracker.holstein_T_int = 0;
#ifndef NDEBUG
    Real holstein_T_initial = tracker.holstein_T_final;
    Real test_holstein_T_int = 0;
#endif
    tracker.holstein_T_final = 0;
    tracker.holstein_G_int = 0;

    Real holTcoef;
    Real holstein_T_int_coef;

    for (int i_lambda = 0; i_lambda < n_lambda; i_lambda++) {
      Real lineshape = tracker.line_shape_function(i_lambda, current_species_T_ratio);
      assert(!std::isnan(lineshape) && lineshape > 0 &&
	     "lineshape must be real and positive");

      Real transfer_probability_lambda_voxel = std::exp(
							-(current_dtau_absorber
							  + (current_dtau_species
							     * lineshape))
							* pathlength
							);
      assert(!std::isnan(transfer_probability_lambda_voxel) &&
	     transfer_probability_lambda_voxel >= 0 &&
	     transfer_probability_lambda_voxel <= 1 &&
	     "transfer probability is a probability.");

      Real transfer_probability_lambda_final = (tracker.transfer_probability_lambda_initial[i_lambda]
						* transfer_probability_lambda_voxel);

      // emission for holstein T int happens at current voxel, use the
      // normalization there
      holTcoef = tracker.weight(i_lambda);

      // holstein_T_int represents a frequency averaged emission in
      // this voxel observed at the start location

      // therefore, we use the lineshape in this voxel to compute holstein T int
      holstein_T_int_coef = (holTcoef
			     * lineshape
			     * tracker.transfer_probability_lambda_initial[i_lambda]
			     * (REAL(1.0) - transfer_probability_lambda_voxel)
			     / (current_abs + lineshape));
      tracker.holstein_T_int += holstein_T_int_coef;
      assert(!std::isnan(tracker.holstein_T_int) && tracker.holstein_T_int >= 0 &&
	     (tracker.holstein_T_int*tracker.line_shape_normalization(current_species_T_ratio) <= tau_species_voxel ||
	      std::abs(tracker.holstein_T_int*tracker.line_shape_normalization(current_species_T_ratio) - tau_species_voxel) < EPS)
	     //  ^^ this allows for small rounding errors
	     && "holstein integral must be between 0 and Delta tau b/c 0<=HolT<=1");

      if (influence) {
	// for single scattering calculation, we want the frequency
	// integrated absorption probability in the start voxel.

	// we need the lineshape at the start voxel
	Real lineshape_at_origin = tracker.line_shape_function(i_lambda,
							       tracker.species_T_ratio_at_origin);
	assert(!std::isnan(lineshape_at_origin) && lineshape_at_origin>0 &&
	       "lineshape must be real and positive");

	// we also need to remove the normalization used before and
	// replace it with the appropriate normalization in this voxel
	Real renormalize_to_origin = tracker.line_shape_normalization(tracker.species_T_ratio_at_origin);

	// holstein T final represents a frequency-averaged absorption and
	// uses lineshape_at_origin
	tracker.holstein_T_final += (holTcoef
				     * renormalize_to_origin // this is needed so T=1 when tau = 0
				     * lineshape_at_origin
				     * transfer_probability_lambda_final);
	assert(!std::isnan(tracker.holstein_T_final)
	       && tracker.holstein_T_final >= 0
	       && tracker.holstein_T_final <= 1
	       && "holstein function represents a probability");

	// finally, the influence coefficient represents emission at
	// the start voxel followed by absorption in the current
	// voxel. Because the emission happens in the start voxel we
	// need to use the start voxel normalization as well as the
	// lineshape in the current voxel.

	// note: even though this is the influence coefficient for
	// emission in the start voxel and absorption in the current
	// voxel, this coefficient should be added to the influence
	// matrix at (row, col) = (start_voxel, current_voxel),
	// because Anderson&Hord1977 say so (it's a result of the
	// expression of the source function in terms of piecewise
	// constant basis functions --- upon substitution back into
	// the original integral equation this inverts the expected
	// relationship of absorber / emitter)
	tracker.holstein_G_int += (holstein_T_int_coef
				   * renormalize_to_origin
				   * lineshape_at_origin);
	assert(!std::isnan(tracker.holstein_G_int)
	       && tracker.holstein_G_int >= 0
	       && tracker.holstein_G_int <= 1
	       && "holstein G integral represents a probability");


#ifndef NDEBUG
	// a test value computing holstein T integral this voxel so we
	// can check that integral(G) = T0 - T1 - a*integral(T)
	
	// We use lineshape_at_origin on top because it's on top in
	// holstein T, which is normalized to the origin, and
	// lineshape on the bottom because the integral is over the
	// current voxel
	test_holstein_T_int += (holTcoef
				* lineshape_at_origin
				* renormalize_to_origin
				* tracker.transfer_probability_lambda_initial[i_lambda]
				* (REAL(1.0) - transfer_probability_lambda_voxel)
				/ (current_abs + lineshape));
#endif
      }
      
      // update the initial values
      tracker.transfer_probability_lambda_initial[i_lambda] = transfer_probability_lambda_final;
    }

    //check that the last element is not contributing too much to the integral
    assert(!((holstein_T_int_coef > 1e-6)
	     && (holstein_T_int_coef/tracker.holstein_T_int > 1e-2))
	   && "wings of line contribute too much to transmission. Increase LAMBDA_MAX in singlet_CFR_tracker.");

    //if holstein T is larger than physically possible due to rounding errors, reduce it to the physical limit
    if (tracker.holstein_T_int > tau_species_voxel)
      tracker.holstein_T_int = tau_species_voxel; /*we don't need to worry
						     about rounding errors
						     here because we checked
						     earlier*/

#ifndef NDEBUG
    if (influence) {
      // check that the integral(G) = (T0-T1) - abs*integral(T) 
      Real holTdiff = (tracker.holstein_G_int
		       + (test_holstein_T_int*current_abs)
		       - (holstein_T_initial-tracker.holstein_T_final));
      assert(std::abs(holTdiff) < EPS
	     && "Integral of holstein G must equal delta holstein T");
    }
#endif
  }

  CUDA_CALLABLE_MEMBER
  void update_tracker_brightness(const Real (&sourcefn_temp)[1], brightness_tracker &tracker) const {
    //bishop formulation
    tracker.brightness += (sourcefn_temp[0] // unitless
			   * emission_g_factor / species_sigma_T_ref * one_over_sqrt_pi// ph / cm2
			   /* ^^^
			      g-factor = (solar flux)*(total cross section) = (solar flux)*(line center sigma @ Tref)*(sqrt(pi)*doppler width @ Tref)
			      therefore (g-factor)/(sqrt(pi)*line center sigma @ Tref) = (solar flux)*(doppler width at Tref), which is what we need.
			   */
			   * branching_ratio // unitless
			   * tracker.holstein_T_int // unitless, carries information about path length in this voxel
			   / REAL(1e9)); // converts to kR, 10^9 ph/cm2/s, see C&H pg 280-282
    assert(!isnan(tracker.brightness) && "brightness must be a real number");
    assert(tracker.brightness>=0 && "brightness must be positive");
  }
  

  
public:

  void set_emission_g_factor(const Real &g) {
    emission_g_factor = g;
  };

  Real get_emission_g_factor() const {
    return emission_g_factor;
  };
  
  //overloads of RT methods
  template<bool influence>
  CUDA_CALLABLE_MEMBER
  void reset_tracker(const int &start_voxel,
		     los<influence,n_elements> &tracker) const {
    if (influence)
      tracker.reset(species_T_ratio(start_voxel));
    else
      tracker.reset(0.0); // start voxel T_ratio doesn't matter for brightness calculations
  }
  
  template<bool influence>
  CUDA_CALLABLE_MEMBER
  void update_tracker_start(const int &current_voxel,
			    const Real & pathlength,
			    los<influence,n_elements> &tracker) const {
    update_tracker_start(species_T_ratio(current_voxel),
			 dtau_species(current_voxel),
			 dtau_absorber(current_voxel),
			 abs(current_voxel),
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

    Real swap_array[1];
    
    parent::interp_voxel_vector(n_interp_points, indices, weights, species_T_ratio_pt, swap_array);
    Real species_T_ratio_interp = swap_array[0];

    parent::interp_voxel_vector(n_interp_points, indices, weights, dtau_species_pt, swap_array);
    Real dtau_species_interp = swap_array[0];

    parent::interp_voxel_vector(n_interp_points, indices, weights, dtau_absorber_pt, swap_array);
    Real dtau_absorber_interp = swap_array[0];
    
    parent::interp_voxel_vector(n_interp_points, indices, weights, abs_pt, swap_array);
    Real abs_interp = swap_array[0];

    update_tracker_start(species_T_ratio_interp,
			 dtau_species_interp,
			 dtau_absorber_interp,
			 abs_interp,
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
    
    //see Bishop1999 for derivation of this formula
    Real coef = domega;
    coef *= tracker.holstein_G_int;
    
    assert(!isnan(coef) && "influence coefficients must be real numbers");
    assert(0<=coef && coef<=1 && "influence coefficients represent transition probabilities");
    // #ifndef __CUDA_ARCH__ // shared memory
    tracker.influence(current_voxel) += coef;
    // #else
    //     atomicAdd(&tracker.influence(current_voxel), coef);
    // #endif
    update_tracker_end(tracker);
  }

  // compute the single scattering from the input tracker
  CUDA_CALLABLE_MEMBER
  void compute_single_scattering(const int &start_voxel, influence_tracker &tracker, bool sun_visible = true) {
    if (!sun_visible) {
      // single scattering point is behind limb, populate the tracker with the appropriate values
      tracker.tau_species_final  = REAL(-1.0);
      tracker.tau_absorber_final = REAL(-1.0);
      tracker.holstein_T_final   = REAL(0.0);
    }
    
    tau_species_single_scattering(start_voxel) = tracker.tau_species_final;//line center optical depth
    assert(!isnan(tau_species_single_scattering(start_voxel))
	   && (tau_species_single_scattering(start_voxel) >= 0
	       || tau_species_single_scattering(start_voxel) == -1)
	   && "optical depth must be real and positive, or -1 if point is behind limb");
    
    tau_absorber_single_scattering(start_voxel) = tracker.tau_absorber_final;
    assert(!isnan(tau_absorber_single_scattering(start_voxel))
	   && (tau_absorber_single_scattering(start_voxel) >= 0
	       ||  tau_absorber_single_scattering(start_voxel) == -1)
	   && "optical depth must be real and positive, or -1 if point is behind limb");
    
    singlescat(start_voxel) = tracker.holstein_T_final;
    assert(!isnan(singlescat(start_voxel))
	   && singlescat(start_voxel) >= 0
	   && "single scattering coefficient must be real and positive");
  }

  using parent::accumulate_influence;

  void pre_solve() {
    influence_matrix.eigen() *= branching_ratio;
  }
  void pre_solve_gpu(); //defined below
protected:
#ifdef __CUDACC__
  template <int NV>
  friend __global__ void singlet_CFR_prepare_for_solution(singlet_CFR<NV> *emission);
#endif
public:

  
  using parent::solve;
  using parent::update_tracker_brightness_interp;
  using parent::update_tracker_brightness_nointerp;
  
  // definition of class members needed for RT
  template<typename C>
  void define(const string &emission_name,
	      const Real &emission_branching_ratio,
	      const Real &species_T_reff, const Real &species_sigma_T_reff,
	      const C &atmosphere,
	      void (boost::type_identity<C>::type::*species_density_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
	      void (boost::type_identity<C>::type::*species_T_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
	      void (boost::type_identity<C>::type::*absorber_density_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
	      Real (boost::type_identity<C>::type::*absorber_sigma_function)(const Real &T) const,
	      const atmo_voxel (&voxels)[N_VOXELS]) {

    strcpy(internal_name, emission_name.c_str());
    branching_ratio     = emission_branching_ratio;
    species_T_ref       = species_T_reff;
    species_sigma_T_ref = species_sigma_T_reff;
    
    for (unsigned int i_voxel=0;i_voxel<N_VOXELS;i_voxel++) {
      (atmosphere.*species_density_function)(voxels[i_voxel],
					     species_density(i_voxel),
					     species_density_pt(i_voxel));
      assert(!isnan(species_density(i_voxel))
	     && species_density(i_voxel) >= 0
	     && "densities must be real and positive");
      assert(!isnan(species_density_pt(i_voxel))
	     && species_density_pt(i_voxel) >= 0
	     && "densities must be real and positive");
      
      Real species_T;
      Real species_T_pt;
      (atmosphere.*species_T_function)(voxels[i_voxel],
				       species_T,
				       species_T_pt);
      species_T_ratio(i_voxel) = species_T_ref/species_T;
      species_T_ratio_pt(i_voxel) = species_T_ref/species_T_pt;
      assert(!isnan(species_T_ratio(i_voxel))
	     && species_T_ratio(i_voxel) >= 0
	     && "temperatures must be real and positive");
      assert(!isnan(species_T_ratio_pt(i_voxel))
	     && species_T_ratio_pt(i_voxel) >= 0
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
      

      absorber_sigma(i_voxel) = (atmosphere.*absorber_sigma_function)(species_T);
      absorber_sigma_pt(i_voxel) = (atmosphere.*absorber_sigma_function)(species_T_pt);
      assert(!isnan(absorber_sigma(i_voxel))
	     && absorber_sigma(i_voxel) >= 0
	     && "cross sections must be real and positive");
      assert(!isnan(absorber_sigma_pt(i_voxel))
	     && absorber_sigma_pt(i_voxel) >= 0
	     && "cross sections must be real and positive");
    }
    
    //define differential optical depths by coefficientwise multiplication
    dtau_species = species_density.eigen().array() * species_sigma_T_ref * species_T_ratio.eigen().array().sqrt();
    dtau_species_pt = species_density_pt.eigen().array() * species_sigma_T_ref * species_T_ratio_pt.eigen().array().sqrt();
    dtau_absorber = absorber_density.eigen().array() * absorber_sigma.eigen().array();
    dtau_absorber_pt = absorber_density_pt.eigen().array() * absorber_sigma_pt.eigen().array();
    abs = dtau_absorber.eigen().array() / dtau_species.eigen().array();
    abs_pt = dtau_absorber_pt.eigen().array() / dtau_species_pt.eigen().array();
    
    parent::reset_solution();

    internal_init=true;
  }

  void save(std::ostream &file, VectorX (*function)(VectorX, int), const int i) const {
    file << "    Species density [cm-3]: "
	 <<      function(species_density.eigen(), i).transpose() << "\n"
      
	 << "    Species single scattering tau: " 
	 <<	 function(tau_species_single_scattering.eigen(), i).transpose() << "\n"
      
	 << "    Species cross section [cm2]: " 
	 <<      function(species_sigma_T_ref*species_T_ratio.eigen().array().sqrt(), i).transpose() << "\n"
      
	 << "    Absorber density [cm-3]: " 
	 <<      function(absorber_density.eigen(), i).transpose() << "\n"
      
	 << "    Absorber single scattering tau: " 
	 <<	 function(tau_absorber_single_scattering.eigen(), i).transpose() << "\n"
      
	 << "    Absorber cross section [cm2]: " 
	 <<	 function(absorber_sigma.eigen(), i).transpose() << "\n"
      
	 << "    Species single scattering source function S0: " 
	 <<	 function(singlescat.eigen(), i).transpose() << "\n"	      
      
	 << "    Source function: " 
	 <<      function(sourcefn.eigen(), i).transpose() << "\n\n";
  }

  void save_brightness(std::ostream &file, const gpu_vector<brightness_tracker> &los_brightness) const {
    // save a list of brightness trackers to file

    VectorX brightness_write_out;
    brightness_write_out.resize(los_brightness.size());
    
    for (int i=0;i<los_brightness.size();i++)
      brightness_write_out[i] = los_brightness[i].brightness;

    file << internal_name
	 << " brightness [kR]: "
	 << brightness_write_out.transpose()
	 << "\n";
  }

#ifdef __CUDACC__
  using parent::device_emission;
  using parent::vector_to_device;
  using parent::copy_trivial_member_to_device;

  void copy_trivial_members_to_device() {
    parent::copy_trivial_members_to_device();
    copy_trivial_member_to_device(branching_ratio, device_emission->branching_ratio);
    copy_trivial_member_to_device(species_T_ref, device_emission->species_T_ref);
    copy_trivial_member_to_device(species_sigma_T_ref, device_emission->species_sigma_T_ref);
    copy_trivial_member_to_device(emission_g_factor, device_emission->emission_g_factor);
  }

  using parent::copy_trivial_member_to_host;
  void copy_trivial_members_to_host() {
    parent::copy_trivial_members_to_host();
    copy_trivial_member_to_host(branching_ratio, device_emission->branching_ratio);
    copy_trivial_member_to_host(species_T_ref, device_emission->species_T_ref);
    copy_trivial_member_to_host(species_sigma_T_ref, device_emission->species_sigma_T_ref);
    copy_trivial_member_to_host(emission_g_factor, device_emission->emission_g_factor);
  }
  
  void copy_to_device_influence() {
    copy_trivial_members_to_device();
    parent::copy_to_device_influence();
  
    //copy the read-only atmosphere arrays
    bool transfer = true;
    vector_to_device(device_emission->species_T_ratio, species_T_ratio, transfer);
    vector_to_device(device_emission->dtau_species, dtau_species, transfer);
    vector_to_device(device_emission->dtau_absorber, dtau_absorber, transfer);
    vector_to_device(device_emission->abs, abs, transfer);
  }

  void copy_to_device_brightness()
  {
    copy_trivial_members_to_device();
    parent::copy_to_device_brightness();

    //free some of the influnce stuff if we used it
    species_T_ratio.free_d_vec();
    dtau_species.free_d_vec();
    dtau_absorber.free_d_vec();
    abs.free_d_vec();

    //copy the stuff we need to do the calculation on the device
    bool transfer = true;
    vector_to_device(device_emission->dtau_species_pt, dtau_species_pt, transfer);
    vector_to_device(device_emission->species_T_ratio_pt, species_T_ratio_pt, transfer);  
    vector_to_device(device_emission->dtau_absorber_pt, dtau_absorber_pt, transfer);
    vector_to_device(device_emission->abs_pt, abs_pt, transfer);
  }
  
  void device_clear() {
    species_density.free_d_vec();
    species_density.free_d_vec();
    species_density_pt.free_d_vec();
    species_T_ratio.free_d_vec();
    species_T_ratio_pt.free_d_vec();
    dtau_species.free_d_vec();
    dtau_species_pt.free_d_vec();

    absorber_density.free_d_vec();
    absorber_density_pt.free_d_vec();
    absorber_sigma.free_d_vec();
    absorber_sigma_pt.free_d_vec();
    dtau_absorber.free_d_vec();
    dtau_absorber_pt.free_d_vec();

    abs.free_d_vec();
    abs_pt.free_d_vec();

    parent::device_clear();
  }
#endif

};

#ifdef __CUDACC__

template<int N_VOXELS>
__global__
void singlet_CFR_prepare_for_solution(singlet_CFR<N_VOXELS> *emission)
{
  //each block prepares one row of the influence matrix
  int i_vox = blockIdx.x;
  
  //convert to kernel from influence (kernel = identity - branching_ratio*influence)
  for (int j_vox = threadIdx.x; j_vox < N_VOXELS; j_vox+=blockDim.x)
    emission->influence_matrix(i_vox,j_vox) *= -emission->branching_ratio;

  if (threadIdx.x == 0) {
    //add the identity matrix
    emission->influence_matrix(i_vox,i_vox) += 1;
    
    //now copy singlescat to sourcefn, preparing for in-place solution
    emission->sourcefn(i_vox) = emission->singlescat(i_vox);
  }
}


template <int N_VOXELS>
void singlet_CFR<N_VOXELS>::pre_solve_gpu() {
  const int n_threads = 32;
  singlet_CFR_prepare_for_solution<<<n_voxels, n_threads>>>(device_emission);
  checkCudaErrors( cudaPeekAtLastError() );
  checkCudaErrors( cudaDeviceSynchronize() );
}

#endif

#endif
