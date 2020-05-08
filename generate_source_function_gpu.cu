//Generate_source_function.cpp -- program to generate a source
//function for comparison with analytic solutions and other models

#include "Real_is_float.h"
#include "cuda_compatibility.h"
#include "atmosphere.h"
#include "RT_grid.h"
#include "RT_gpu.h"
#include "grid_plane_parallel.h"
#include "grid_spherical_azimuthally_symmetric.h"

int main(int argc, char* argv[]) {

  //define the physical atmosphere
  Real exobase_temp = 200;//K
  krasnopolsky_temperature temp(exobase_temp);

  Real H_exobase_density = 5e5;// cm-3
  Real CO2_exobase_density = 2e8;//cm-3
  chamb_diff_1d atm(H_exobase_density,
		    CO2_exobase_density,
		    temp);

  //use holstein functions to compute influence integrals
  holstein_approx hol; //this is the most time consuming part of startup ~0.4s

  //define the RT grid
  static const int n_emissions=2;
  vector<string> emission_names = {"H Lyman alpha", "H Lyman beta"};

  // static const int n_radial_boundaries = 40;
  // static const int n_rays_phi = 6;
  // static const int n_rays_theta = 12;
  // plane_parallel_grid<n_radial_boundaries,n_rays_phi,n_rays_theta> grid;

  static const int n_radial_boundaries = 40;
  static const int n_sza_boundaries = 20;/*20 for 10 deg increments with szamethod_uniform*/
  static const int n_rays_phi = 6;
  static const int n_rays_theta = 12;
  spherical_azimuthally_symmetric_grid<n_radial_boundaries,n_sza_boundaries,n_rays_phi,n_rays_theta> grid;
  //grid.save_intersections = true;
  //grid.rmethod = grid.rmethod_altitude;
  grid.rmethod = grid.rmethod_lognH;
  //grid.szamethod = grid.szamethod_uniform;
  grid.szamethod = grid.szamethod_uniform_cos;
  
  grid.setup_voxels(atm);
  grid.setup_rays();
  


  RT_grid<n_emissions,
	  typeof(grid),
	  holstein_approx> RT(emission_names, grid, hol);

  //solve for H lyman alpha
  RT.define_emission("H Lyman alpha",
		       1.0,
		       atm,
		       &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lya,
		       &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lya);
  //solve for H lyman beta
  RT.define_emission("H Lyman beta",
  		       lyman_beta_branching_ratio,
  		       atm,
  		       &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lyb,
  		       &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lyb);

  //solve for the source function
  //  RT.save_influence = true;
  RT.generate_S();
  
  //now print out the output
  //RT.save_S("test/test_source_function.dat");

  //simulate a fake observation
  observation obs(emission_names);

  vector<Real> g = {lyman_alpha_typical_g_factor, lyman_beta_typical_g_factor};
  obs.emission_g_factors = g;
  
  // std::cout << "lyman alpha g factor is: " << lyman_alpha_typical_g_factor << std::endl;
  // std::cout << "lyman alpha tau=1 brightness at " << exobase_temp << " K : "
  // 	    <<   (lyman_alpha_typical_g_factor/
  // 		  (lyman_alpha_line_center_cross_section_coef/std::sqrt(exobase_temp))
  // 		  /1e9)
  // 	    << " kR" << std::endl;
  
  // std::cout << "lyman beta g factor is: " << lyman_beta_typical_g_factor << std::endl;
  // std::cout << "lyman beta tau=1 brightness at " << exobase_temp << " K : "
  // 	    <<   (lyman_beta_typical_g_factor/
  // 		  (lyman_beta_line_center_cross_section_coef/std::sqrt(exobase_temp))
  // 		  /1e6)
  // 	    << " R" << std::endl;

  
  vector<int> sizes = {10,100,300,600,1200,2400};
  //vector<int> sizes = {600};

  for (auto&& size: sizes) {
    std::cout << "simulating image size "<< size << "x" << size << ":" << std::endl;
    Real dist = 30*rMars;
    obs.fake(dist,30,size);
    observation obs_nointerp = obs;
    RT.brightness_gpu(obs);

    my_clock save_clk;
    save_clk.start();
    string fname = "test/test_brightness_gpu" + std::to_string(size) + "x" + std::to_string(size) + ".dat";
    obs.save_brightness(fname);  
    save_clk.stop();
    save_clk.print_elapsed("saving file takes ");
    
    std::cout << std::endl;
}

  
  
  return 0; 
}

