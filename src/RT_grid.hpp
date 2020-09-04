 //RT_grid.h -- basic RT grid definitions
//             specifics for different geometries are in other files

#ifndef __RT_grid_H
#define __RT_grid_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include <iostream> // for file output and dialog
#include <cmath>    // for cos and sin
#include "constants.hpp" // basic parameter definitions
#include "my_clock.hpp"
#include "atmo_vec.hpp"
#include "boundaries.hpp"
#include "influence.hpp"
#include "emission.hpp"
#include "observation.hpp"
#include "los_tracker.hpp"
#include <Eigen/Dense> //matrix solvers
#include <Eigen/StdVector>
#include <string>
#include <cassert>
#include <boost/type_traits/type_identity.hpp> //for type deduction in define_emission

//structure to hold the atmosphere grid
template<int N_EMISSIONS, typename grid_type, typename influence_type>
struct RT_grid {
  static const int n_emissions = N_EMISSIONS;//number of emissions to evaluate at each point in the grid
  emission<grid_type::n_voxels> emissions[n_emissions];

  grid_type grid;//this stores all of the geometrical info

  //need to refactor so this lives inside each emission
  //the influence function
  influence_type transmission;

  //initialization parameters
  bool all_emissions_init;

  //GPU interface
  typedef RT_grid<N_EMISSIONS,grid_type,influence_type> RT_grid_type;
  RT_grid_type *d_RT=NULL; //pointer to the GPU partner of this object
  void RT_to_device();
  void RT_to_device_brightness();
  void RT_to_device_influence();
  void emissions_to_device_influence();
  void emissions_to_device_brightness();
  void emissions_influence_to_host();
  void emissions_solved_to_host();
  
  RT_grid() : all_emissions_init(false) { }
  
  RT_grid(const string (&emission_names)[n_emissions]) : RT_grid()
  {
    set_names(emission_names);
  }
  ~RT_grid() {
#ifdef __CUDACC__
    if(d_RT!=NULL)
      checkCudaErrors(cudaFree(d_RT));
#endif
  }

  void set_names(const string (&emission_names)[n_emissions]) {
    for (int i_emission=0; i_emission < n_emissions; i_emission++)
      emissions[i_emission].name = emission_names[i_emission];
  }

  template<typename C>
  void define_emission(const string &emission_name,
		       const Real &emission_branching_ratio,
		       const Real &species_T_ref, const Real &species_sigma_T_ref,
		       const C &atmosphere,
		       void (boost::type_identity<C>::type::*species_density_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
		       void (boost::type_identity<C>::type::*species_T_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
		       void (boost::type_identity<C>::type::*absorber_density_function)(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const,
		       Real (boost::type_identity<C>::type::*absorber_sigma_function)(const Real &T) const) {

    //find emission name (dumb search but n_emissions is small and these are called infrequently)
    int n;
    for (n=0;n<n_emissions;n++) {
      if (emissions[n].name == emission_name)
	break;
    }
    if (n == n_emissions) {
      assert(false && "can't find emission name in define_emission");
    } else {
      assert((0<=n&&n<n_emissions) && "attempt to set invalid emission in define_emitter");
      
      emissions[n].define(emission_branching_ratio,
			  species_T_ref, species_sigma_T_ref,
			  atmosphere,
			  species_density_function,
			  species_T_function,
			  absorber_density_function,
			  absorber_sigma_function,
			  grid.voxels);
      
      all_emissions_init = true;
      for (int i_emission=0;i_emission<n_emissions;i_emission++)
	if (!emissions[i_emission].init)
	  all_emissions_init = false;
    }
  }

  
  template <typename R>
  CUDA_CALLABLE_MEMBER
  void voxel_traverse(const atmo_vector &v,
		      void (RT_grid::*function)(boundary_intersection_stepper<grid_type::n_dimensions,
						                              grid_type::n_max_intersections>& , R& ),
		      R &retval)
  {
    assert(all_emissions_init && "!nitialize grid and influence function before calling voxel_traverse.");
  
    //get boundary intersections
    boundary_intersection_stepper<grid_type::n_dimensions,
				  grid_type::n_max_intersections> stepper;
    grid.ray_voxel_intersections(v, stepper);

    if (stepper.boundaries.size() == 0)
      return;
    
    stepper.origin();

    while(stepper.inside) {
      //call the desired function in each voxel
      (this->*function)(stepper, retval);
    
      stepper.next();
    }
  }

  CUDA_CALLABLE_MEMBER
  void influence_update(boundary_intersection_stepper<grid_type::n_dimensions,
			                              grid_type::n_max_intersections>& stepper,
			influence_tracker<n_emissions,grid_type::n_voxels>& temp_influence) {
    //update the influence matrix for each emission

    temp_influence.update_start(stepper,emissions);
      
    for (int i_emission=0; i_emission < n_emissions; i_emission++) {
      
      //see Bishop1999 for derivation of this formula
      Real coef = stepper.vec.ray.domega;

      coef *= temp_influence.line[i_emission].holstein_G_int;

      assert(!isnan(coef) && "influence coefficients must be real numbers");
      assert(0<=coef && coef<=1 && "influence coefficients represent transition probabilities");
      
      temp_influence.influence[i_emission][stepper.current_voxel] += coef;
    
    }

    temp_influence.update_end();

  }

  CUDA_CALLABLE_MEMBER
  void get_single_scattering_optical_depths(boundary_intersection_stepper<grid_type::n_dimensions,
					                                  grid_type::n_max_intersections>& stepper,
					    influence_tracker<n_emissions,grid_type::n_voxels>& temp_influence)
  {
    temp_influence.update_start(stepper,emissions);
    temp_influence.update_end();
  }
  
  CUDA_CALLABLE_MEMBER
  void get_single_scattering(const atmo_point &pt, influence_tracker<n_emissions,grid_type::n_voxels>& temp_influence) {

    if (pt.z<0&&pt.x*pt.x+pt.y*pt.y<grid.rmin*grid.rmin) {
      //if the point is behind the planet, no single scattering
      for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	emissions[i_emission].tau_species_single_scattering(pt.i_voxel) = -1;
	emissions[i_emission].tau_absorber_single_scattering(pt.i_voxel) = -1;
	emissions[i_emission].singlescat(pt.i_voxel)=0.0;
      }
    } else {
      atmo_vector vec = atmo_vector(pt, grid.sun_direction);
      voxel_traverse(vec,
		     &RT_grid::get_single_scattering_optical_depths,
		     temp_influence);
      
      for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	emissions[i_emission].tau_species_single_scattering(pt.i_voxel) = temp_influence.line[i_emission].tau_species_final;//line center optical depth
	assert(!isnan(emissions[i_emission].tau_species_single_scattering(pt.i_voxel))
	       && emissions[i_emission].tau_species_single_scattering(pt.i_voxel) >= 0
	       && "optical depth must be real and positive");

	emissions[i_emission].tau_absorber_single_scattering(pt.i_voxel) = temp_influence.line[i_emission].tau_absorber_final;
	assert(!isnan(emissions[i_emission].tau_absorber_single_scattering(pt.i_voxel))
	       && emissions[i_emission].tau_absorber_single_scattering(pt.i_voxel) >= 0
	       && "optical depth must be real and positive");

	
	emissions[i_emission].singlescat(pt.i_voxel) = temp_influence.line[i_emission].holstein_T_final;
	assert(!isnan(emissions[i_emission].singlescat(pt.i_voxel))
	       && emissions[i_emission].singlescat(pt.i_voxel) >= 0
	       && "single scattering coefficient must be real and positive");

      }
    }
  }

  void solve() {
    for (int i_emission=0;i_emission<n_emissions;i_emission++) {
      emissions[i_emission].influence_matrix *= emissions[i_emission].branching_ratio;

      MatrixX kernel = MatrixX::Identity(grid.n_voxels,grid.n_voxels);
      kernel -= emissions[i_emission].influence_matrix;

      emissions[i_emission].sourcefn=kernel.partialPivLu().solve(emissions[i_emission].singlescat);//partialPivLu has multithreading support

      // // iterative solution.
      // Real err = 1;
      // int it = 0;
      // VectorX sourcefn_old(grid.n_voxels);
      // sourcefn_old = emissions[i_emission].singlescat;
      // while (err > ABS && it < 500) {
      // 	emissions[i_emission].sourcefn = emissions[i_emission].singlescat + emissions[i_emission].influence_matrix * sourcefn_old;

      // 	err=((emissions[i_emission].sourcefn-sourcefn_old).array().abs()/sourcefn_old.array()).maxCoeff();
      // 	sourcefn_old = emissions[i_emission].sourcefn;
      // 	it++;
      // }
      // std::cout << "For " << emissions[i_emission].name << std::endl;
      // std::cout << "  Scattering up to order: " << it << " included.\n";
      // std::cout << "  Error at final order is: " << err << " .\n";
      
      for (int i=0;i<grid.n_voxels;i++)
	emissions[i_emission].log_sourcefn(i) = emissions[i_emission].sourcefn(i) == 0 ? -1e5 : log(emissions[i_emission].sourcefn(i));
      emissions[i_emission].solved=true;
    }
  }
  //gpu stuff is defined in RT_gpu.cu
  void transpose_influence_gpu();
  void solve_emission_gpu(emission<grid_type::n_voxels> & emiss);
  void solve_gpu();


  //generate source functions on the grid
  void generate_S() {
  
    //start timing
    my_clock clk;
    clk.start();

    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      emissions[i_emission].influence_matrix.setZero();
    
    atmo_vector vec;
    influence_tracker<n_emissions,grid_type::n_voxels> temp_influence;
    Real max_tau_species = 0;
  
#pragma omp parallel for firstprivate(vec,temp_influence) shared(max_tau_species,std::cout) default(none)
    for (int i_vox = 0; i_vox < grid.n_voxels; i_vox++) {
      Real omega = 0.0; // make sure sum(domega) = 4*pi
    
      //now integrate outward along the ray grid:
      for (int i_ray=0; i_ray < grid.n_rays; i_ray++) {
	vec = atmo_vector(grid.voxels[i_vox].pt, grid.rays[i_ray]);
	omega += vec.ray.domega;

	temp_influence.reset(emissions, i_vox);
	voxel_traverse(vec, &RT_grid::influence_update, temp_influence);
	for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	  double rowsum = 0.0;
	  for (int j_vox = 0; j_vox < grid.n_voxels; j_vox++) {
	    emissions[i_emission].influence_matrix(i_vox,j_vox) += temp_influence.influence[i_emission][j_vox];
	    rowsum += emissions[i_emission].influence_matrix(i_vox,j_vox);
	  }
	  assert(0.0 <= rowsum && rowsum <= 1.0 && "row represents scattering probability from this voxel");
	}
	
	if (temp_influence.max_tau_species > max_tau_species)
	  max_tau_species = temp_influence.max_tau_species;
      }

    assert(std::abs(omega - 1.0) < ABS && "omega must = 4*pi\n");
    
      //now compute the single scattering function:
      temp_influence.reset(emissions, i_vox);
      get_single_scattering(grid.voxels[i_vox].pt, temp_influence);
      if (temp_influence.max_tau_species > max_tau_species)
	max_tau_species = temp_influence.max_tau_species;
    }
  
    //solve for the source function
    solve();

    //std::cout << "max_tau_species = " << max_tau_species << std::endl;
  
    // print time elapsed
    clk.stop();
    clk.print_elapsed("source function generation takes ");
    std::cout << std::endl;
    
    return;
  }
  void generate_S_gpu();
  
  void save_influence(const string fname = "test/influence_matrix.dat") {
    std::ofstream file(fname);
    if (file.is_open())
      for (int i_emission=0;i_emission<n_emissions;i_emission++)
	file << "Here is the influence matrix for " << emissions[i_emission].name <<":\n" 
	     << emissions[i_emission].influence_matrix << "\n\n";
  }

  void save_S(const string fname) {
    grid.save_S(fname,emissions,n_emissions);
  }


  template <typename V>
  CUDA_CALLABLE_MEMBER
  Real interp_array(const int *indices, const Real *weights, const V &arr) const {
    Real retval=0;
    for (int i=0;i<grid.n_interp_points;i++)
      retval+=weights[i]*arr[indices[i]];
    return retval;
  }
  CUDA_CALLABLE_MEMBER
  void interp(const int &ivoxel, const atmo_point &pt, interpolated_values<n_emissions> &retval) const {

    int indices[grid_type::n_interp_points];
    Real weights[grid_type::n_interp_points];

    grid.interp_weights(ivoxel,pt,indices,weights);
    
    for (int i_emission=0;i_emission<n_emissions;i_emission++) {
      retval.dtau_species_interp[i_emission]  = interp_array(indices, weights,
							     emissions[i_emission].dtau_species_pt);
      retval.species_T_interp[i_emission]     = interp_array(indices, weights,
							     emissions[i_emission].species_T_pt);
      retval.dtau_absorber_interp[i_emission] = interp_array(indices, weights,
							     emissions[i_emission].dtau_absorber_pt);
      retval.abs_interp[i_emission]           = interp_array(indices, weights,
							     emissions[i_emission].abs_pt);
      
      retval.sourcefn_interp[i_emission]      = interp_array(indices, weights,
							     emissions[i_emission].sourcefn);
    }
  }
  
  //interpolated brightness routine
  CUDA_CALLABLE_MEMBER
  void brightness(const atmo_vector &vec, const Real (&g)[n_emissions],
		  brightness_tracker<n_emissions> &los,
		  const int n_subsamples=20) const {
    assert(n_subsamples!=1 && "choose either 0 or n>1 voxel subsamples.");
    
    boundary_intersection_stepper<grid_type::n_dimensions,
				  grid_type::n_max_intersections> stepper;
    grid.ray_voxel_intersections(vec, stepper);

    los.reset(emissions);

    atmo_point pt;
    interpolated_values<n_emissions> interp_vals;
    
    //brightness is zero if we do not intersect the grid
    if (stepper.boundaries.size() == 0)
      return;

    int n_subsamples_distance = n_subsamples;
    if (n_subsamples==0)
      n_subsamples_distance=2;
    
    for(unsigned int i_bound=1;i_bound<stepper.boundaries.size();i_bound++) {
      //add option to step in increments of optical depth here?
      Real d_start=stepper.boundaries[i_bound-1].distance;
      Real d_step=(stepper.boundaries[i_bound].distance-d_start)/(n_subsamples_distance-1);

      //account for small rounding errors in boundary crossings
      const Real eps = ABS;
      d_start += eps/2.*d_step;
      d_step *= 1-eps;

      int current_voxel = stepper.boundaries[i_bound-1].entering;

      
      for (int i_step=1;i_step<n_subsamples_distance;i_step++) {

	pt = vec.extend(d_start+i_step*d_step);
	
	if (n_subsamples!=0) {
	  interp(current_voxel,pt,interp_vals);
	}
	for (int i_emission=0;i_emission<n_emissions;i_emission++) {

	  Real dtau_species_temp, species_T_temp;
	  Real dtau_absorber_temp, abs_temp;
	  Real sourcefn_temp;
	  
	  if (n_subsamples == 0) {
	    dtau_species_temp  = emissions[i_emission].dtau_species[current_voxel];
	    species_T_temp     = emissions[i_emission].species_T[current_voxel];
	    dtau_absorber_temp = emissions[i_emission].dtau_absorber[current_voxel];
	    abs_temp           = emissions[i_emission].abs[current_voxel];
	    sourcefn_temp      = emissions[i_emission].sourcefn[current_voxel];
	  } else {
	    dtau_species_temp  = interp_vals.dtau_species_interp[i_emission];
	    species_T_temp     = interp_vals.species_T_interp[i_emission];
	    dtau_absorber_temp = interp_vals.dtau_absorber_interp[i_emission];
	    abs_temp           = interp_vals.abs_interp[i_emission];
	    sourcefn_temp      = interp_vals.sourcefn_interp[i_emission];
	  } 

	  //should really rewrite this whole function to make changes
	  //to stepper (?) to capture interpolation and use voxel traverse
	  los.line[i_emission].update_start(species_T_temp, emissions[i_emission].species_T_ref,
					    dtau_species_temp,
					    dtau_absorber_temp,
					    abs_temp,
					    d_step);

	  //bishop formulation
	  sourcefn_temp = sourcefn_temp*emissions[i_emission].branching_ratio;

	  los.brightness[i_emission] += (sourcefn_temp
					 * los.line[i_emission].holstein_T_int
					 / emissions[i_emission].species_sigma_T_ref); 
	  //                                                     ^^^^^^^^^^^^^^^^^^^
	  //                                                     sigma was multiplied
	  //                                                     out of S0 so we must
	  //                                                     divide by it here
	  assert(!isnan(los.brightness[i_emission]) && "brightness must be a real number");
	  assert(los.brightness[i_emission]>=0 && "brightness must be positive");
	  
	  los.line[i_emission].update_end();
	}
      }
    }

    if (stepper.exits_bottom) {
      for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	los.line[i_emission].tau_absorber_final = -1;
      }
    }

    //convert to kR
    for (int i_emission=0;i_emission<n_emissions;i_emission++) {
      los.brightness[i_emission] *= emissions[i_emission].branching_ratio;
      los.brightness[i_emission] *= g[i_emission]/1e9; //megaphoton/cm2/s * 1e-3 = kR, see C&H pg 280-282
    }
  }

  void brightness(observation<n_emissions> &obs, const int n_subsamples=20) const {
    assert(obs.size()>0 && "there must be at least one observation to simulate!");
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      assert(obs.emission_g_factors[i_emission] != 0. && "set emission g factors before simulating brightness");

    my_clock clk;
    clk.start();

    obs.reset_output();

#pragma omp parallel for shared(obs) firstprivate(n_subsamples) default(none)
    for(int i=0; i<obs.size(); i++)
      brightness(obs.get_vec(i),obs.emission_g_factors,
		 obs.los[i],
		 n_subsamples);
    
    clk.stop();
    clk.print_elapsed("brightness calculation takes ");
  }

  void brightness_nointerp(observation<n_emissions> &obs) const {
    brightness(obs,0);
  }

  //hooks for porting to gpu
  void brightness_gpu(observation<n_emissions> &obs, const int n_subsamples=20);
  
};

//import the CUDA code if NVCC is the compiler
#ifdef __CUDACC__
#include "RT_gpu.cu"
#endif

#endif
