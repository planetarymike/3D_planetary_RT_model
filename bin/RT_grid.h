 //RT_grid.h -- basic RT grid definitions
//             specifics for different geometries are in other files

#ifndef __RT_grid_H
#define __RT_grid_H

#include <iostream> // for file output and dialog
#include <cmath>    // for cos and sin
#include "constants.h" // basic parameter definitions
#include "my_clock.h"
#include "atmo_vec.h"
#include "boundaries.h"
#include "influence.h"
#include "emission.h"
#include "observation.h"
#include <Eigen/Dense> //matrix solvers
#include <Eigen/StdVector>
#include <string>
#include <cassert>

//structure to hold the atmosphere grid
template<typename grid_type, typename influence_type>
struct RT_grid {
  int n_emissions;//number of emissions to evaluate at each point in the grid
  vector<emission> emissions;

  grid_type grid;//this stores all of the geometrical info
  static const int n_dimensions = grid_type::n_dimensions;

  //need to refactor so this lives inside each emission
  //the influence function
  influence_type transmission;

  //initialization parameters
  bool all_emissions_init;
  
  RT_grid(const vector<string> &emission_names,
	  const grid_type &gridd,
	  const influence_type &transmissionn) :
    n_emissions(emission_names.size()), grid(gridd), transmission(transmissionn)
  {
    all_emissions_init = false;
    
    emissions.resize(n_emissions);
    for (int i_emission=0; i_emission < n_emissions; i_emission++) {
      emissions[i_emission].name = emission_names[i_emission];
      emissions[i_emission].resize(grid.n_voxels);
    }
  }
  
  template<typename C>
  void define_emission(string emission_name,
		       double emission_branching_ratio,
		       C &atmosphere,
		       double (C::*species_density_function)(const atmo_point) ,
		       double (C::*species_sigma_function)(const atmo_point),
		       double (C::*absorber_density_function)(const atmo_point),
		       double (C::*absorber_sigma_function)(const atmo_point)) {


    struct find_name : std::unary_function<emission, bool> {
      string name;
      find_name(string name): name(name) { }
      bool operator()(emission const& e) const {
        return e.name == name;
      }
    };
    
    std::vector<emission>::iterator it = std::find_if(emissions.begin(), emissions.end(), 
						    find_name(emission_name));
    if (it == emissions.end()) {
      assert(false && "can't find emission name in define_emission");
    } else {
      int n = std::distance(emissions.begin(), it);
      assert((0<=n&&n<n_emissions) && "attempt to set invalid emission in define_emitter");

      emissions[n].define(emission_branching_ratio,
			  atmosphere,
			  species_density_function,
			  species_sigma_function,
			  absorber_density_function,
			  absorber_sigma_function,
			  grid.pts);
      
      all_emissions_init = true;
      for (auto&& emiss: emissions)
	if (!emiss.init)
	  all_emissions_init = false;
    }
  }

  
  template <typename R>
  void voxel_traverse(const atmo_vector &v,
		      void (RT_grid::*function)(boundary_intersection_stepper<n_dimensions>& ,R& ),
		      R &retval)
  {
    assert(all_emissions_init && "!nitialize grid and influence function before calling voxel_traverse.");
  
    //get boundary intersections
    static thread_local boundary_intersection_stepper<n_dimensions> stepper;
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


  struct tau_tracker {
    int n_emissions;
    vector<double> tau_species_initial;
    vector<double> tau_species_final;
    vector<double> tau_absorber_initial;
    vector<double> tau_absorber_final;
    double max_tau_species;

    tau_tracker(int n_emissionss)
      : n_emissions(n_emissionss)
    {
      tau_species_initial.resize(n_emissions,0.);
      tau_species_final.resize(n_emissions,0.);
      tau_absorber_initial.resize(n_emissions,0.);
      tau_absorber_final.resize(n_emissions,0.);
      
      max_tau_species = 0;
    }

    void reset() {
      std::fill(tau_species_initial.begin(), tau_species_initial.end(), 0);
      std::fill(tau_species_final.begin(), tau_species_final.end(), 0);
      std::fill(tau_absorber_initial.begin(), tau_absorber_initial.end(), 0);
      std::fill(tau_absorber_final.begin(), tau_absorber_final.end(), 0);
    }
    
    void check_max_tau() {
      for (int i_emission=0; i_emission < n_emissions; i_emission++) {
	if (tau_species_final[i_emission] > max_tau_species)
	  max_tau_species = tau_species_final[i_emission];
      }
    }
    
    void update_start(boundary_intersection_stepper<n_dimensions>& stepper,
		      vector<emission> &emissions) {
      for (int i_emission=0; i_emission < n_emissions; i_emission++) {
	tau_species_final[i_emission] = ( tau_species_initial[i_emission]
					  + (emissions[i_emission].dtau_species(stepper.current_voxel)
					     * stepper.pathlength));
	tau_absorber_final[i_emission] = ( tau_absorber_initial[i_emission]
					   + (emissions[i_emission].dtau_absorber(stepper.current_voxel)
					      * stepper.pathlength)); 
      }
    }
    
    void update_end() {
      for (int i_emission=0; i_emission < n_emissions; i_emission++) {
	tau_species_initial[i_emission]=tau_species_final[i_emission];
	tau_absorber_initial[i_emission]=tau_absorber_final[i_emission];
      }
      check_max_tau();
    }
  };

  void influence_update(boundary_intersection_stepper<n_dimensions>& stepper, tau_tracker& los) {
    //update the influence matrix for each emission

    los.update_start(stepper,emissions);
      
    for (int i_emission=0; i_emission < n_emissions; i_emission++) {
      
      //see Bishop1999 for derivation of this formula
      double coef = stepper.vec.ray.domega;

      //bishop formulation
      coef *= ((transmission.T(los.tau_species_initial[i_emission])
		- transmission.T(los.tau_species_final[i_emission]))
	       
      	       *exp(-0.5*(los.tau_absorber_initial[i_emission]
      			  +los.tau_absorber_final[i_emission])));
      
      emissions[i_emission].influence_matrix(stepper.start_voxel,
					     stepper.current_voxel) += coef;
    
    }

    los.update_end();

  }

  void get_single_scattering_optical_depths(boundary_intersection_stepper<n_dimensions>& stepper,
					    double& max_tau_species)
  {
    for (auto&& emiss: emissions) {
      emiss.tau_species_single_scattering(stepper.start_voxel) += (emiss.dtau_species(stepper.current_voxel)
								   * stepper.pathlength);
      emiss.tau_absorber_single_scattering(stepper.start_voxel) += (emiss.dtau_absorber(stepper.current_voxel)
								    * stepper.pathlength);
      
      double tscomp = emiss.tau_species_single_scattering(stepper.start_voxel);
      max_tau_species = tscomp > max_tau_species ? tscomp : max_tau_species;    
    }
  }
  

  void get_single_scattering(const atmo_point &pt, double &max_tau_species) {

    if (pt.z<0&&pt.x*pt.x+pt.y*pt.y<grid.rmin*grid.rmin) {
      //if the point is behind the planet, no single scattering
      for (auto&& emiss: emissions)
	emiss.singlescat(pt.i_voxel)=0.0;
    } else {
      atmo_vector vec = atmo_vector(pt, grid.sun_direction);
      voxel_traverse(vec,
		     &RT_grid::get_single_scattering_optical_depths,
		     max_tau_species);
      
      for (auto&& emiss: emissions)
	emiss.singlescat(pt.i_voxel) = ( (transmission.T(emiss.tau_species_single_scattering(pt.i_voxel))
					  * exp(-emiss.tau_absorber_single_scattering(pt.i_voxel)) )
					 / emiss.species_sigma(pt.i_voxel));
    }
  }

  void solve() {
    for (auto&& emiss: emissions) {
      emiss.influence_matrix *= emiss.branching_ratio;

      MatrixXd kernel = MatrixXd::Identity(grid.n_voxels,grid.n_voxels);
      kernel -= emiss.influence_matrix;

      emiss.sourcefn=kernel.partialPivLu().solve(emiss.singlescat);//partialPivLu has multithreading support

      // // iterative solution.
      // double err = 1;
      // int it = 0;
      // VectorXd sourcefn_old(grid.n_voxels);
      // sourcefn_old = emiss.singlescat;
      // while (err > 1e-6 && it < 500) {
      // 	emiss.sourcefn = emiss.singlescat + emiss.influence_matrix * sourcefn_old;

      // 	err=((emiss.sourcefn-sourcefn_old).array().abs()/sourcefn_old.array()).maxCoeff();
      // 	sourcefn_old = emiss.sourcefn;
      // 	it++;
      // }
      // std::cout << "For " << emiss.name << std::endl;
      // std::cout << "  Scattering up to order: " << it << " included.\n";
      // std::cout << "  Error at final order is: " << err << " .\n";
      
      for (int i=0;i<grid.n_voxels;i++)
	emiss.log_sourcefn(i) = emiss.sourcefn(i) == 0 ? -1e5 : log(emiss.sourcefn(i));
    emiss.solved=true;

    }
  }

  //generate source functions on the grid
  void generate_S() {
  
    //start timing
    my_clock clk;
    clk.start();

    for (auto&& emiss: emissions) {
      emiss.influence_matrix.setZero();
      emiss.tau_species_single_scattering.setZero();
      emiss.tau_absorber_single_scattering.setZero();
    }
    
    atmo_vector vec;
    double max_tau_species = 0;
  
#pragma omp parallel for firstprivate(vec) shared(max_tau_species,std::cout) default(none)
    for (int i_pt = 0; i_pt < grid.n_voxels; i_pt++) {
      double omega = 0.0; // make sure sum(domega) = 4*pi
    
      //now integrate outward along the ray grid:
      for (int i_ray=0; i_ray < grid.n_rays; i_ray++) {
	vec = atmo_vector(grid.pts[i_pt], grid.rays[i_ray]);
	omega += vec.ray.domega;

	tau_tracker los(n_emissions);
      
	voxel_traverse(vec, &RT_grid::influence_update, los);

	if (los.max_tau_species > max_tau_species)
	  max_tau_species = los.max_tau_species;
      }

      if ((omega - 1.0 > 1e-6) || (omega - 1.0 < -1e-6)) {
	std::cout << "omega != 4*pi\n";
	throw(10);
      }
    
      //now compute the single scattering function:
      get_single_scattering(grid.pts[i_pt], max_tau_species);
    
    }
  
    //solve for the source function
    solve();

    //std::cout << "max_tau_species = " << max_tau_species << std::endl;
  
    // print time elapsed
    clk.stop();
    clk.print_elapsed();
    
    return;
  }

  void save_influence(const string fname = "test/influence_matrix.dat") {
    std::ofstream file(fname);
    if (file.is_open())
      for (auto&& emiss: emissions)
	file << "Here is the influence matrix for " << emiss.name <<":\n" 
	     << emiss.influence_matrix << "\n\n";
  }

  void save_S(const string fname) {
    grid.save_S(fname,emissions);
  }





  //interpolation support
  class interpolated_values {
  public:
    vector<double> dtau_species_interp;
    vector<double> dtau_absorber_interp;
    vector<double> abs_interp;
    vector<double> sourcefn_interp;
    
    interpolated_values(const int n_emissions) {
      dtau_species_interp.resize(n_emissions);
      dtau_absorber_interp.resize(n_emissions);
      abs_interp.resize(n_emissions);
      sourcefn_interp.resize(n_emissions);      
    }
  };

  double interp_array(const int *indices, const double *weights, const VectorXd &arr) const {
    double retval=0;
    for (int i=0;i<grid.n_interp_points;i++)
      retval+=weights[i]*arr(indices[i]);
    return retval;
  }

  void interp(const int &ivoxel, const atmo_point &pt, interpolated_values &retval) const {

    int indices[grid.n_interp_points];
    double weights[grid.n_interp_points];

    grid.interp_weights(ivoxel,pt,indices,weights);
    
    for (int i_emission=0;i_emission<n_emissions;i_emission++) {
      retval.dtau_species_interp[i_emission]  = exp(interp_array(indices, weights,
								 emissions[i_emission].log_dtau_species));
      retval.dtau_absorber_interp[i_emission] = exp(interp_array(indices, weights,
								 emissions[i_emission].log_dtau_absorber));
      retval.abs_interp[i_emission]           = exp(interp_array(indices, weights,
								 emissions[i_emission].log_abs));
      retval.sourcefn_interp[i_emission]      = exp(interp_array(indices, weights,
								 emissions[i_emission].log_sourcefn));
    }
  };
  
  struct brightness_tracker : tau_tracker {
    vector<double> brightness;

    brightness_tracker(int n_emissionss) : tau_tracker(n_emissionss) {
      brightness.resize(n_emissionss,0.);
    }

    void reset() {
      tau_tracker::reset();
      std::fill(brightness.begin(), brightness.end(), 0);
    }
  };


  //interpolated brightness routine
  brightness_tracker brightness(const atmo_vector &vec, const vector<double> &g, const int n_subsamples=5) const {
    assert(n_subsamples!=1 && "choose either 0 or n>1 voxel subsamples.");
    
    static thread_local boundary_intersection_stepper<n_dimensions> stepper;
    grid.ray_voxel_intersections(vec, stepper);
    static thread_local brightness_tracker los(n_emissions);
    los.reset();

    static thread_local atmo_point pt;
    static thread_local interpolated_values interp_vals(n_emissions);
    
    //brightness is zero if we do not intersect the grid
    if (stepper.boundaries.size() == 0)
      return los;

    int n_subsamples_distance = n_subsamples;
    if (n_subsamples==0)
      n_subsamples_distance=2;
    
    for(unsigned int i_bound=1;i_bound<stepper.boundaries.size();i_bound++) {
      double d_start=stepper.boundaries[i_bound-1].distance;
      //first factor here accounts for small rounding errors in boundary crossings
      double d_step=(1-1e-6)*(stepper.boundaries[i_bound].distance-d_start)/(n_subsamples_distance-1);
      int current_voxel = stepper.boundaries[i_bound-1].entering;

      
      for (int i_step=1;i_step<n_subsamples_distance;i_step++) {

	pt = vec.extend(d_start+i_step*d_step);

	if (n_subsamples!=0) {
	  interp(current_voxel,pt,interp_vals);
	}
	for (int i_emission=0;i_emission<n_emissions;i_emission++) {

	  double dtau_species_temp, dtau_absorber_temp, sourcefn_temp;
	  //double abs_temp;
	  
	  if (n_subsamples == 0) {
	    dtau_species_temp  = emissions[i_emission].dtau_species(current_voxel);
	    dtau_absorber_temp = emissions[i_emission].dtau_absorber(current_voxel);
	    //abs_temp           = emissions[i_emission].abs(current_voxel);
	    sourcefn_temp      = emissions[i_emission].sourcefn(current_voxel);
	  } else {
	    dtau_species_temp  = interp_vals.dtau_species_interp[i_emission];
	    dtau_absorber_temp = interp_vals.dtau_absorber_interp[i_emission];
	    //abs_temp           = interp_vals.abs_interp[i_emission];
	    sourcefn_temp      = interp_vals.sourcefn_interp[i_emission];
	  } 


	  los.tau_species_final[i_emission] = ( los.tau_species_initial[i_emission]
						+ dtau_species_temp*d_step);
	  los.tau_absorber_final[i_emission] = ( los.tau_absorber_initial[i_emission]
						 + dtau_absorber_temp*d_step);


	  //bishop formulation
	  los.brightness[i_emission] += (sourcefn_temp
					 
	  				 *emissions[i_emission].branching_ratio
					 
	  				 *(transmission.Tint(los.tau_species_final[i_emission])
	  				   - transmission.Tint(los.tau_species_initial[i_emission]))
					 
	  				 *exp(-0.5*(los.tau_absorber_initial[i_emission]
	  					    +los.tau_absorber_final[i_emission])));

	  los.tau_species_initial[i_emission]=los.tau_species_final[i_emission];
	  los.tau_absorber_initial[i_emission]=los.tau_absorber_final[i_emission];
	}
      }
    }

    if (stepper.exits_bottom) 
      los.tau_absorber_final = vector<double>(n_emissions,-1);

    //convert to kR
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      los.brightness[i_emission] *= g[i_emission]/1e9; //megaphoton/cm2/s * 1e-3 = kR, see C&H pg 280-282

    return los;
  }

  vector<brightness_tracker> brightness(const vector<atmo_vector> &vecs,
					const vector<double> &g,
					const int n_subsamples=5) const {
    vector<brightness_tracker> retval;
    
    retval.resize(vecs.size(),brightness_tracker(n_emissions));

#pragma omp parallel for shared(retval) firstprivate(vecs,g,n_subsamples) default(none)
    for(unsigned int i=0; i<vecs.size(); i++)
      retval[i] = brightness(vecs[i],g,n_subsamples);
    
    return retval;
  }

  void brightness(observation &obs, const int n_subsamples=5) const {
    assert(obs.size()>0 && "there must be at least one observation to simulate!");
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      assert(obs.emission_g_factors[i_emission] != 0. && "set emission g factors before simulating brightness");
    
    vector<brightness_tracker> los = brightness(obs.get_vecs(), obs.emission_g_factors, n_subsamples);
    for (int i=0;i<obs.size();i++) {
      obs.brightness[i]   = los[i].brightness;
      obs.tau_species[i]  = los[i].tau_species_final;
      obs.tau_absorber[i] = los[i].tau_absorber_final;
    }
  }

  void brightness_nointerp(observation &obs) const {
    brightness(obs,0);
  }

  //hooks for porting to gpu
  void traverse_gpu(observation &obs, const int n_subsamples=5);
  
};


#endif
