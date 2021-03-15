//RT_grid.h -- basic RT grid definitions
//             specifics for different geometries are in other files

#ifndef __RT_grid_H
#define __RT_grid_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include <iostream> // for file output and dialog
#include <cmath>    // for cos and sin
#include "my_clock.hpp"
#include "atmo_vec.hpp"
#include "boundaries.hpp"
#include "observation.hpp"
#include <string>
#include <cassert>
#include <type_traits>

//structure to hold the atmosphere grid
template<typename emission_type, int N_EMISSIONS, typename grid_type>
struct RT_grid {
  static const int n_emissions = N_EMISSIONS; //number of emissions to evaluate at each point in the grid
  emission_type *emissions[n_emissions];// all emissions must be the same type

  grid_type grid;//this stores all of the geometrical info

  //initialization parameters
  bool all_emissions_init;

  //GPU interface
  typedef RT_grid<emission_type,
		  N_EMISSIONS,
		  grid_type> RT_grid_type;
  RT_grid_type *d_RT=NULL; //pointer to the GPU partner of this object
  void RT_to_device();
  void RT_to_device_brightness();
  void RT_to_device_influence();
  void emissions_to_device_influence();
  void emissions_to_device_brightness();
  void emissions_influence_to_host();
  void emissions_solved_to_host();
  void device_clear();
  
  RT_grid(grid_type gridd,
	  emission_type *emissionss[n_emissions]) :
    grid(gridd) {

    all_emissions_init = true;
    for (int i_emission=0;i_emission<n_emissions;i_emission++) {
      emissions[i_emission] = emissionss[i_emission];
      all_emissions_init &= emissions[i_emission]->init();
    }
  }

  ~RT_grid() {
#ifdef __CUDACC__
    device_clear();
#endif
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
			typename emission_type::influence_tracker (&temp_influence)[n_emissions]) {
    //update the influence matrix for each emission
    
    for (int i_emission=0; i_emission < n_emissions; i_emission++)
      emissions[i_emission]->update_tracker_influence(stepper.current_voxel,
						      stepper.pathlength,
						      stepper.vec.ray.domega,
						      temp_influence[i_emission]);
  }

  CUDA_CALLABLE_MEMBER
  void get_single_scattering_optical_depths(boundary_intersection_stepper<grid_type::n_dimensions,
					                                  grid_type::n_max_intersections>& stepper,
					    typename emission_type::influence_tracker (&temp_influence)[n_emissions])
  {
    for (int i_emission=0; i_emission < n_emissions; i_emission++) {
      //update influence functions for this voxel
      emissions[i_emission]->update_tracker_start(stepper.current_voxel,
						  stepper.pathlength,
						  temp_influence[i_emission]);
      emissions[i_emission]->update_tracker_end(temp_influence[i_emission]);
    }
  }
  
  CUDA_CALLABLE_MEMBER
  void get_single_scattering(const atmo_point &pt, typename emission_type::influence_tracker (&temp_influence)[n_emissions]) {
    if (pt.z<0&&pt.x*pt.x+pt.y*pt.y<grid.rmin*grid.rmin) {
      //if the point is behind the planet, no single scattering
      for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	temp_influence[i_emission].tau_species_final  = REAL(-1.0);
	temp_influence[i_emission].tau_absorber_final = REAL(-1.0);
	temp_influence[i_emission].holstein_T_final   = REAL(0.0);
	emissions[i_emission]->compute_single_scattering(pt.i_voxel, temp_influence[i_emission]);
      }
    } else {
      // compute the RT towards the sun and pass to the emissions
      atmo_vector vec;
      vec.ptvec(pt, grid.sun_direction);
      voxel_traverse(vec,
		     &RT_grid::get_single_scattering_optical_depths,
		     temp_influence);
      
      for (int i_emission=0;i_emission<n_emissions;i_emission++)
	emissions[i_emission]->compute_single_scattering(pt.i_voxel, temp_influence[i_emission]);
    }
  }

  void solve() {
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      emissions[i_emission]->solve();
  }
  //gpu stuff is defined in RT_gpu.cu
  void solve_gpu();


  //generate source functions on the grid
  void generate_S() {
  
    //start timing
    my_clock clk;
    clk.start();

    atmo_vector vec;

    Real max_tau_species = 0;
  
#pragma omp parallel firstprivate(vec) shared(max_tau_species,emissions) default(none)
    {
      typename emission_type::influence_tracker temp_influence[n_emissions];
      for (int i_emission = 0; i_emission < n_emissions; i_emission++)
	temp_influence[i_emission].init();
      
#pragma omp for
      for (int i_vox = 0; i_vox < grid.n_voxels; i_vox++) {
	
	Real omega = 0.0; // make sure sum(domega) = 4*pi
	
	//now integrate outward along the ray grid:
	for (int i_ray=0; i_ray < grid.n_rays; i_ray++) {
	  
	  // reset vec and temp_influence for this ray
	  vec.ptray(grid.voxels[i_vox].pt, grid.rays[i_ray]);
	  omega += vec.ray.domega;
	  for (int i_emission=0;i_emission<n_emissions;i_emission++)
	    emissions[i_emission]->reset_tracker(i_vox, temp_influence[i_emission]);
	  
	  // accumulate influence along the ray
	  voxel_traverse(vec, &RT_grid::influence_update, temp_influence);
	  
	  // pack the contributions back into emissions
	  for (int i_emission = 0; i_emission < n_emissions; i_emission++) {
	    emissions[i_emission]->accumulate_influence(i_vox, temp_influence[i_emission]);
	    if (temp_influence[i_emission].max_tau_species > max_tau_species)
#pragma omp atomic write
	      max_tau_species = temp_influence[i_emission].max_tau_species;
	  }
	}
	
	assert(std::abs(omega - 1.0) < ABS && "omega must = 4*pi\n");
	
	// now compute the single scattering function:
	for (int i_emission = 0; i_emission < n_emissions; i_emission++)
	  emissions[i_emission]->reset_tracker(i_vox, temp_influence[i_emission]);
	get_single_scattering(grid.voxels[i_vox].pt, temp_influence);
	for (int i_emission = 0; i_emission < n_emissions; i_emission++)
	  if (temp_influence[i_emission].max_tau_species > max_tau_species)
	    max_tau_species = temp_influence[i_emission].max_tau_species;
      }
      
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
  
  void save_influence(const string fname = "test/influence_matrix.dat") const {
    std::ofstream file(fname);
    if (file.is_open())
      for (int i_emission=0;i_emission<n_emissions;i_emission++)
	emissions[i_emission]->save_influence(file);
  }

  void save_S(const string fname) const {
    grid.save_S(fname, emissions, n_emissions);
  }

  //interpolated brightness routine
  CUDA_CALLABLE_MEMBER
  void brightness(const atmo_vector &vec, const Real (&g)[n_emissions],
		  typename emission_type::brightness_tracker* (&los)[n_emissions], // array of pointers to los trackers
		  const int n_subsamples=5) const {
    assert(n_subsamples!=1 && "choose either 0 or n>1 voxel subsamples.");
    
    boundary_intersection_stepper<grid_type::n_dimensions,
				  grid_type::n_max_intersections> stepper;
    grid.ray_voxel_intersections(vec, stepper);

    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      emissions[i_emission]->reset_tracker(0/*voxel number doesn't matter
					      for brightness calculation*/,
					   *los[i_emission]);
    atmo_point pt;

    // interpolation stuff
    int indices[grid_type::n_interp_points];
    Real weights[grid_type::n_interp_points];
    int indices_1d[2*grid_type::n_dimensions];
    Real weights_1d[grid_type::n_dimensions];
    
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
      d_start += REAL(0.5)*eps*d_step;
      d_step *= REAL(1.0)-eps;

      int current_voxel = stepper.boundaries[i_bound-1].entering;

      for (int i_step=1;i_step<n_subsamples_distance;i_step++) {

	pt = vec.extend(d_start+i_step*d_step);
	
	if (n_subsamples!=0)
	  grid.interp_weights(current_voxel,pt,indices,weights,indices_1d,weights_1d);
	
	for (int i_emission=0;i_emission<n_emissions;i_emission++)
	  if (n_subsamples == 0)
	    emissions[i_emission]->update_tracker_brightness_nointerp(current_voxel,
								      d_step,
								      *los[i_emission]);
	  else
	    emissions[i_emission]->update_tracker_brightness_interp(grid_type::n_interp_points,
								    indices,
								    weights,
								    d_step,
								    *los[i_emission]);
      }
    }

    if (stepper.exits_bottom)
      for (int i_emission=0;i_emission<n_emissions;i_emission++)
	los[i_emission]->tau_absorber_final = -1;

    // convert to kR
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      los[i_emission]->brightness *= g[i_emission]/REAL(1e9); //megaphoton/cm2/s * 1e-3 = kR, see C&H pg 280-282
  }

  void brightness(observation<emission_type, n_emissions> &obs, const int n_subsamples=5) const {
    assert(obs.size()>0 && "there must be at least one observation to simulate!");
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      assert(obs.emission_g_factors[i_emission] != 0. && "set emission g factors before simulating brightness");

    my_clock clk;
    clk.start();
    
#pragma omp parallel for shared(obs) firstprivate(n_subsamples) default(none)
    for(int i=0; i<obs.size(); i++) {
      typename emission_type::brightness_tracker *los[n_emissions];
      for (int i_emission=0;i_emission<n_emissions;i_emission++)
	los[i_emission] = &obs.los[i_emission][i];
      brightness(obs.get_vec(i), obs.emission_g_factors,
		 los,
		 n_subsamples);
    }
    clk.stop();
    clk.print_elapsed("brightness calculation takes ");
  }
  
  void brightness_nointerp(observation<emission_type, n_emissions> &obs) const {
    brightness(obs,0);
  }
  
  //hooks for porting to gpu
  void brightness_gpu(observation<emission_type, N_EMISSIONS> &obs, const int n_subsamples=5);
};

//import the CUDA code if NVCC is the compiler
#ifdef __CUDACC__
#include "RT_gpu.cu"
#endif

#endif
