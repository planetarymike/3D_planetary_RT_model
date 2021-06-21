//generate_source_function.cpp -- program to generate a source
//function for comparison with analytic solutions and other models

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "atm/temperature.hpp"
#include "atm/chamb_diff_1d.hpp"
#include "RT_grid.hpp"
#include "emission/singlet_CFR.hpp"
#include "grid_plane_parallel.hpp"
#include "grid_spherical_azimuthally_symmetric.hpp"

int main(__attribute__((unused)) int argc, __attribute__((unused)) char* argv[]) {

  // read input if provided
  Real input_exobase_density, input_exobase_temp;
  if (argc == 1) {
    input_exobase_density = 5e5;// cm-3
    input_exobase_temp = 200;//K
  } else {
    input_exobase_density = atof(argv[1]);
    input_exobase_temp = atof(argv[2]);
  }
  
  //define the physical atmosphere
  Real exobase_temp = input_exobase_temp;//K
  Real H_T_ref = exobase_temp;//K
  krasnopolsky_temperature temp(exobase_temp);

  Real H_exobase_density = input_exobase_density;//cm-3
  Real CO2_exobase_density = 2e8;//cm-3
  hydrogen_density_parameters H_thermosphere;
  chamb_diff_1d atm(H_exobase_density,
		    CO2_exobase_density,
		    &temp,
		    &H_thermosphere);
  // // fix temperature to the exobase temp, eliminate CO2 absorption to compare with JY
  // atm.temp_dependent_sH=false;
  // atm.constant_temp_sH=exobase_temp;
  //atm.no_CO2_absorption = true;
  atm.save("test/test_atmosphere.dat");
  


  // define the geometry of the grid
  //#define PLANE_PARALLEL
#ifdef PLANE_PARALLEL
  static const int n_radial_boundaries = 40;
  static const int n_rays_theta = 6;
  static const int n_rays_phi = 12;
  typedef plane_parallel_grid<n_radial_boundaries,
			      n_rays_theta> grid_type;
  atm.spherical = false;
  
  grid_type grid;
  grid.rmethod = grid.rmethod_log_n_species;
#else
  static const int n_radial_boundaries = 40;
  static const int n_sza_boundaries = 20;/*20 for 10 deg increments with szamethod_uniform*/
  static const int n_rays_theta = 6;
  static const int n_rays_phi = 12;
  typedef spherical_azimuthally_symmetric_grid<n_radial_boundaries,
  					       n_sza_boundaries,
  					       n_rays_phi,
  					       n_rays_theta> grid_type;

  grid_type grid;
  
  //grid.rmethod = grid.rmethod_altitude;
  grid.rmethod = grid.rmethod_log_n_species;
  //grid.rmethod = grid.rmethod_log_n_species_tau_absorber;
  //grid.szamethod = grid.szamethod_uniform; //requires CONEABS = 1e-2 in Real.hpp
  grid.szamethod = grid.szamethod_uniform_cos;

  grid.raymethod_theta = grid.raymethod_theta_uniform;
#endif
  
  grid.setup_voxels(atm);
  grid.setup_rays();


  //define the emissions to be solved for
  static const int n_emissions = 2;

  //solve for H lyman alpha
  typedef singlet_CFR<grid_type::n_voxels> emission_type;
  emission_type lyman_alpha;
  lyman_alpha.define("H Lyman alpha",
		     /*emission branching ratio = */1.0,
		     H_T_ref, atm.sH_lya(H_T_ref),
		     atm,
		     &chamb_diff_1d::n_species_voxel_avg,   &chamb_diff_1d::Temp_voxel_avg,
		     &chamb_diff_1d::n_absorber_voxel_avg,  &chamb_diff_1d::sCO2_lya,
		     grid.voxels);
  //solve for H lyman beta
  emission_type lyman_beta;
  lyman_beta.define("H Lyman beta",
		    /*emission branching ratio = */lyman_beta_branching_ratio,
		    H_T_ref, atm.sH_lyb(H_T_ref),
		    atm,
		    &chamb_diff_1d::n_species_voxel_avg,   &chamb_diff_1d::Temp_voxel_avg,
		    &chamb_diff_1d::n_absorber_voxel_avg,  &chamb_diff_1d::sCO2_lyb,
		    grid.voxels);

  emission_type *emissions[n_emissions] = {&lyman_alpha, &lyman_beta};

  // now set up the RT object
  RT_grid<emission_type, n_emissions, grid_type> RT(grid, emissions);
  std::cout << "size of RT grid: " << sizeof(RT) << "\n";

  //solve for the source function
#ifndef __CUDACC__
  RT.generate_S();
  // RT.save_influence();
#else
  //RT.generate_S();
  RT.generate_S_gpu();
  RT.generate_S_gpu();
#endif
  //now print out the output
  RT.save_S("test/test_source_function.dat");

  string sfn_name = "test/test_source_function_";
  sfn_name += std::to_string(n_radial_boundaries) + "x";
#ifndef PLANE_PARALLEL
  sfn_name += std::to_string(n_sza_boundaries) + "x";
#endif
  sfn_name += std::to_string(n_rays_phi) + "x";
  sfn_name += std::to_string(n_rays_theta) + ".dat";
  RT.save_S(sfn_name);
  
#ifndef PLANE_PARALLEL
  //simulate a fake observation
  observation<emission_type, n_emissions> obs(emissions);
  observation<emission_type, n_emissions> obs_nointerp(emissions);

  lyman_alpha.set_emission_g_factor(lyman_alpha_typical_g_factor);
  lyman_beta.set_emission_g_factor(lyman_beta_typical_g_factor);
  
  std::cout << "lyman alpha g factor is: " << lyman_alpha_typical_g_factor << std::endl;
  std::cout << "lyman alpha line center cross section coef is: "
  	    <<  lyman_alpha_line_center_cross_section_coef << std::endl;
  std::cout << "lyman alpha tau=1 brightness at " << exobase_temp << " K : "
  	    <<   (lyman_alpha_typical_g_factor/
  		  (lyman_alpha_line_center_cross_section_coef/std::sqrt(exobase_temp))
  		  /1e9)
  	    << " kR" << std::endl;
  
  std::cout << "lyman beta g factor is: " << lyman_beta_typical_g_factor << std::endl;
  std::cout << "lyman beta line center cross section coef is: "
  	    <<  lyman_beta_line_center_cross_section_coef << std::endl;
  std::cout << "lyman beta tau=1 brightness at " << exobase_temp << " K : "
  	    <<   (lyman_beta_typical_g_factor/
  		  (lyman_beta_line_center_cross_section_coef/std::sqrt(exobase_temp))
  		  /1e6)
  	    << " R" << std::endl;
  std::cout << std::endl;

  Real dist = 30*rMars;

#ifndef __CUDACC__
  //CPU-only code
  int size = 600;
  obs.fake(dist, 30, size);
  obs_nointerp.fake(dist, 30, size);
  
  std::cout << "Performing brightness calculation without interpolation...\n";
  RT.brightness_nointerp(obs_nointerp);
  obs_nointerp.save_brightness("test/test_brightness_nointerp.dat");
  std::cout << "Performing interpolated brightness calculation...\n";
  RT.brightness(obs);
  obs.save_brightness("test/test_brightness.dat");
#else
  //GPU code
  vector<int> sizes = {10,100,300,600/*,1200,2400*/};

  for (auto&& size: sizes) {
    std::cout << "simulating image size "<< size << "x" << size << ":" << std::endl;
    obs.fake(dist,30,size);
    RT.brightness_gpu(obs);

    my_clock save_clk;
    save_clk.start();
    string fname = "test/test_brightness_gpu" + std::to_string(size) + "x" + std::to_string(size) + ".dat";
    obs.save_brightness(fname);  
    save_clk.stop();
    save_clk.print_elapsed("saving file takes ");
    
    std::cout << std::endl;
  }
#endif
#endif
  
  return 0; 
}
