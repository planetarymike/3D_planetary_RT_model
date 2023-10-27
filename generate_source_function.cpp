//generate_source_function.cpp -- program to generate a source
//function for comparison with analytic solutions and other models

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "atm/temperature.hpp"
#include "atm/chamb_diff_1d.hpp"
#include "grid_plane_parallel.hpp"
#include "grid_spherical_azimuthally_symmetric.hpp"
#include "RT_grid.hpp"
#include "emission/singlet_CFR.hpp"
#include "emission/O_1026.hpp"
#include "emission/H_lyman_multiplet_test.hpp"
#include "emission/H_lyman_multiplet.hpp"

//#define GENERATE_O_1026
//#define H_LYMAN_MULTIPLET_TEST
//#define H_LYMAN_SINGLET_TEST
//#define NO_SIM_BRIGHTNESS

int main(__attribute__((unused)) int argc, __attribute__((unused)) char* argv[]) {

  // read input if provided
  Real input_exobase_density, input_exobase_temp;
  if (argc == 1) {
#ifdef GENERATE_O_1026
    input_exobase_density = 2e7;// cm-3
#else
    input_exobase_density = 5e5;// cm-3
#endif
    input_exobase_temp = 200;//K
  } else {
    input_exobase_density = atof(argv[1]);
    input_exobase_temp = atof(argv[2]);
  }
  
  //define the physical atmosphere
  Real exobase_temp = input_exobase_temp;//K
  krasnopolsky_temperature temp(exobase_temp);

  Real species_exobase_density = input_exobase_density;//cm-3
  Real CO2_exobase_density = 2e8;//cm-3
#ifdef GENERATE_O_1026
  oxygen_density_parameters species_thermosphere;
  Real min_exosphere_density = 10;
  Real rmin = rMars+100e5;
  int method = thermosphere_exosphere::method_nspmin_nCO2exo;
  //  // to specify max altitude instead of min density:
  //  Real rmax = rMars+1500e5;
  //  int method = thermosphere_exosphere::method_rmax_nCO2rmin;
  //  min_exosphere_density = rmax;
#else
  hydrogen_density_parameters species_thermosphere;
  Real min_exosphere_density = 10;
  Real rmin = rMars+80e5;
  int method = thermosphere_exosphere::method_nspmin_nCO2exo;
#endif

  chamb_diff_1d atm(/* rmin = */ rmin,
		    /* rexo = */ rMars+200e5,
		    /* rmaxx_or_nspmin = */ min_exosphere_density,
		    /* rmindifussion = */ rmin,
		    species_exobase_density,
		    CO2_exobase_density,
		    &temp,
		    &species_thermosphere,
		    method);
  // // fix temperature to the exobase temp, eliminate CO2 absorption to compare with JY
  // atm.temp_dependent_sH=false;
  // atm.constant_temp_sH=exobase_temp;
  // atm.no_CO2_absorption = true;

#if defined(GENERATE_O_1026)
  atm.save("test/test_atmosphere_O_1026.dat");
#elif defined(H_LYMAN_MULTIPLET_TEST)
  atm.save("test/test_atmosphere_H_multiplet.dat");
#elif defined(H_LYMAN_SINGLET_TEST)
  atm.save("test/test_atmosphere_H_singlet.dat");
#else
  atm.save("test/test_atmosphere.dat");
#endif

  // define the geometry of the grid
  //#define PLANE_PARALLEL
#ifdef PLANE_PARALLEL
  static const int n_radial_boundaries = 40;
  static const int n_rays_theta = 6;
  typedef plane_parallel_grid<n_radial_boundaries,
			      n_rays_theta> grid_type;
  atm.spherical = false;
  
  grid_type grid;
  grid.rmethod = grid.rmethod_log_n_species;
#else
  static const int n_radial_boundaries = 40;
  static const int n_sza_boundaries = 20;/*20 for 10 deg increments with szamethod_uniform*/
  static const int n_rays_theta = 7;
  static const int n_rays_phi = 12;
  typedef spherical_azimuthally_symmetric_grid<n_radial_boundaries,
  					       n_sza_boundaries,
  					       n_rays_theta,
  					       n_rays_phi> grid_type;

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
#if defined(GENERATE_O_1026)
  static const int n_emissions = 1;

  //solve for O 1026
  typedef O_1026_emission<grid_type::n_voxels> emission_type;
  emission_type oxygen_1026;
  oxygen_1026.define("O_1026",
		     atm,
		     &chamb_diff_1d::n_species_voxel_avg,   
		     &chamb_diff_1d::Temp_voxel_avg,
		     &chamb_diff_1d::n_absorber_voxel_avg,
		     grid.voxels);
  oxygen_1026.set_solar_brightness(1.69e-3); /* ph/cm2/s/Hz, 
						solar line center brightness at Lyman beta, 
						solar minimum */

  emission_type *emissions[n_emissions] = {&oxygen_1026};
#elif defined(H_LYMAN_MULTIPLET_TEST)

  static const int n_emissions = 1;

  // true doublets for Lyman alpha and beta
  typedef H_lyman_multiplet<grid_type::n_voxels> emission_type;
  emission_type lyman_emission;
  if (atm.no_CO2_absorption)
    lyman_emission.set_CO2_absorption_off();
  if (not atm.temp_dependent_sH)
    lyman_emission.set_constant_temp_RT(exobase_temp);
  lyman_emission.define("H Lyman alpha and beta",
			atm,
			&chamb_diff_1d::n_species_voxel_avg,
			&chamb_diff_1d::Temp_voxel_avg,
			&chamb_diff_1d::n_absorber_voxel_avg,
			grid.voxels);
  lyman_emission.set_solar_brightness(lyman_alpha_flux_Mars_typical, /* ph/cm2/s/Hz, solar brightness */
				      lyman_beta_flux_Mars_typical); 

  emission_type *emissions[n_emissions] = {&lyman_emission};

#elif defined(H_LYMAN_SINGLET_TEST)

  static const int n_emissions = 1;

  // singlet treatment in standard multiplet code
  typedef H_lyman_singlet<grid_type::n_voxels> emission_type;
  emission_type lyman_emission;
  if (atm.no_CO2_absorption)
    lyman_emission.set_CO2_absorption_off();
  if (not atm.temp_dependent_sH)
    lyman_emission.set_constant_temp_RT(exobase_temp);
  lyman_emission.define("H Lyman alpha and beta",
			atm,
			&chamb_diff_1d::n_species_voxel_avg,
			&chamb_diff_1d::Temp_voxel_avg,
			&chamb_diff_1d::n_absorber_voxel_avg,
			grid.voxels);
  lyman_emission.set_solar_brightness(lyman_alpha_flux_Mars_typical, /* ph/cm2/s/Hz, solar brightness */
				      lyman_beta_flux_Mars_typical); 

  emission_type *emissions[n_emissions] = {&lyman_emission};


#else
  // normal singlet CFR treatment
  
  //static const int n_emissions = 1;
  static const int n_emissions = 2;

  //solve for H lyman alpha
  typedef singlet_CFR<grid_type::n_voxels> emission_type;
  emission_type lyman_alpha;
  lyman_alpha.define("H Lyman alpha",
		     /*emission branching ratio = */1.0,
		     exobase_temp, atm.sH_lya(exobase_temp),
		     atm,
		     &chamb_diff_1d::n_species_voxel_avg,   &chamb_diff_1d::Temp_voxel_avg,
		     &chamb_diff_1d::n_absorber_voxel_avg,  &chamb_diff_1d::sCO2_lya,
		     grid.voxels);
  //solve for H lyman beta
  emission_type lyman_beta;
  lyman_beta.define("H Lyman beta",
  		    /*emission branching ratio = */lyman_beta_branching_ratio,
  		    exobase_temp, atm.sH_lyb(exobase_temp),
  		    atm,
  		    &chamb_diff_1d::n_species_voxel_avg,   &chamb_diff_1d::Temp_voxel_avg,
  		    &chamb_diff_1d::n_absorber_voxel_avg,  &chamb_diff_1d::sCO2_lyb,
		    grid.voxels);

  //emission_type *emissions[n_emissions] = {&lyman_alpha};
  emission_type *emissions[n_emissions] = {&lyman_alpha, &lyman_beta};
#endif
  
  // now set up the RT object
  RT_grid<emission_type, n_emissions, grid_type> RT(grid, emissions);
  std::cout << "size of RT grid: " << sizeof(RT) << "\n";

  //solve for the source function
#ifndef __CUDACC__
  RT.generate_S();
#else
  //RT.generate_S();
  RT.generate_S_gpu();
  RT.generate_S_gpu();
#endif
  //now print out the output

#if defined(GENERATE_O_1026)
  string sfn_name_tag = "_O_1026";
#elif defined(H_LYMAN_MULTIPLET_TEST)
  string sfn_name_tag = "_H_multiplet";
#elif defined(H_LYMAN_SINGLET_TEST)
  string sfn_name_tag = "_H_singlet";
#else
  string sfn_name_tag = "";
#endif
  
  RT.save_S("test/test_source_function"+sfn_name_tag+".dat");

  string sfn_name_tag_dims = sfn_name_tag+"_";
  sfn_name_tag_dims += std::to_string(n_radial_boundaries) + "x";
#ifndef PLANE_PARALLEL
  sfn_name_tag_dims += std::to_string(n_sza_boundaries) + "x";
#endif
  sfn_name_tag_dims += std::to_string(n_rays_theta);
#ifndef PLANE_PARALLEL
  sfn_name_tag_dims += "x";
  sfn_name_tag_dims += std::to_string(n_rays_phi);
#endif
  RT.save_S("test/test_source_function"+sfn_name_tag_dims+".dat");

  RT.save_influence("test/influence_matrix"+sfn_name_tag+".dat");
  RT.save_influence("test/influence_matrix"+sfn_name_tag_dims+".dat");

#if !defined(NO_SIM_BRIGHTNESS) && !defined(PLANE_PARALLEL)
  //simulate a fake observation
  observation<emission_type, n_emissions> obs(emissions);
  observation<emission_type, n_emissions> obs_nointerp(emissions);

  Real dist = 30*rMars;

#ifdef GENERATE_O_1026
  Real angle = grid.rmax / dist * 180/pi;
  Vector3 loc = {0.,-1.,0.};
#else
  Real angle = 30;
  Vector3 loc = {0.,-1.0,0.};
#endif

#if !defined(GENERATE_O_1026) && !defined(H_LYMAN_MULTIPLET_TEST) && !defined(H_LYMAN_SINGLET_TEST)
  lyman_alpha.set_emission_g_factor(lyman_alpha_typical_g_factor);
  lyman_beta.set_emission_g_factor(lyman_beta_typical_g_factor);
  
  std::cout << "lyman alpha solar brightness is: " << lyman_alpha_flux_Mars_typical << " ph/cm2/s/Hz" <<std::endl;
  std::cout << "lyman alpha g factor is: " << lyman_alpha_typical_g_factor << std::endl;
  std::cout << "lyman alpha line center cross section coef is: "
  	    <<  lyman_alpha_line_center_cross_section_coef << std::endl;
  std::cout << "lyman alpha tau=1 brightness at " << exobase_temp << " K : "
  	    <<   (lyman_alpha_typical_g_factor/
  		  (lyman_alpha_line_center_cross_section_coef/std::sqrt(exobase_temp))
  		  /1e9)
  	    << " kR" << std::endl;
  
  std::cout << "lyman beta solar brightness is: " << lyman_beta_flux_Mars_typical << " ph/cm2/s/Hz" <<std::endl;
  std::cout << "lyman beta g factor is: " << lyman_beta_typical_g_factor << std::endl;
  std::cout << "lyman beta line center cross section coef is: "
  	    <<  lyman_beta_line_center_cross_section_coef << std::endl;
  std::cout << "lyman beta tau=1 brightness at " << exobase_temp << " K : "
  	    <<   (lyman_beta_typical_g_factor/
  		  (lyman_beta_line_center_cross_section_coef/std::sqrt(exobase_temp))
  		  /1e6)
  	    << " R" << std::endl;
  std::cout << std::endl;
#endif

#ifndef __CUDACC__
  //CPU-only code
  int size = 600;

  obs.fake(dist, angle, size, loc);
  obs_nointerp.fake(dist, angle, size, loc);
  
  std::cout << "Performing brightness calculation without interpolation...\n";
  RT.brightness_nointerp(obs_nointerp);
  obs_nointerp.save_brightness("test/test_brightness"+sfn_name_tag+"_nointerp.dat");
  std::cout << "Performing interpolated brightness calculation...\n";
  RT.brightness(obs);
  obs.save_brightness("test/test_brightness"+sfn_name_tag+".dat");
#else
  //GPU code
  vector<int> sizes = {10,100,300,600/*,1200,2400*/};

  for (auto&& size: sizes) {
    std::cout << "simulating image size "<< size << "x" << size << ":" << std::endl;
    obs.fake(dist,angle,size,loc);
    RT.brightness_gpu(obs);

    my_clock save_clk;
    save_clk.start();
    string fname = "test/test_brightness_gpu" + sfn_name_tag + "_"+ std::to_string(size) + "x" + std::to_string(size) + ".dat";
    obs.save_brightness(fname);  
    save_clk.stop();
    save_clk.print_elapsed("saving file takes ");
    
    std::cout << std::endl;
  }
#endif
#endif
  return 0; 
}
