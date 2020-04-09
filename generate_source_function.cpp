//generate_source_function.cpp -- program to generate a source
//function for comparison with analytic solutions and other models

#include "atmosphere.h"
#include "RT_plane_parallel.h"
#include "RT_spherical_azimuthally_symmetric.h"

int main(int argc, char* argv[]) {

  //define the physical atmosphere
  double exobase_temp = 200;//K
  krasnopolsky_temperature temp(exobase_temp);

  double H_exobase_density = 5e5;// cm-3
  double CO2_exobase_density = 2e8;//cm-3
  chamb_diff_1d atm(H_exobase_density,
		    CO2_exobase_density,
		    temp);

  //use holstein functions to compute influence integrals
  holstein_approx hol;

  //define the RT grid
  vector<string> emission_names = {"H Lyman alpha", "H Lyman beta"};

  // plane_parallel_grid grid(n_emissions, hol);
  // grid.setup_voxels(160, atm);
  // grid.setup_rays(10);

  spherical_azimuthally_symmetric_grid grid(emission_names, hol);

  //grid.save_intersections = true;
  //grid.rmethod = grid.rmethod_altitude;
  grid.rmethod = grid.rmethod_lognH;
  //grid.szamethod = grid.szamethod_uniform;
  grid.szamethod = grid.szamethod_uniform_cos;
  

  // grid.setup_voxels(10, 4/*20 for 10 deg increments with sza_uniform*/, atm);
  // grid.setup_rays(6, 12);
  grid.setup_voxels(40, 20/*20 for 10 deg increments with sza_uniform*/, atm);
  grid.setup_rays(6, 12);
  
  //solve for H lyman alpha
  grid.define_emission("H Lyman alpha",
		       1.0,
		       atm,
		       &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lya,
		       &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lya);
  //solve for H lyman beta
  grid.define_emission("H Lyman beta",
		       lyman_beta_branching_ratio,
		       atm,
		       &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lyb,
		       &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lyb);

  //solve for the source function
  //  grid.save_influence = true;
  grid.generate_S();
  
  //now print out the output
  grid.save_S("test/test_source_function.dat");

  //simulate a fake observation
  observation obs(grid.emission_names);

  vector<double> g = {lyman_alpha_typical_g_factor, lyman_beta_typical_g_factor};
  obs.emission_g_factors = g;
  
  // std::cout << "lyman alpha g factor is:" << lyman_alpha_typical_g_factor << std::endl;
  // std::cout << "lyman beta g factor is:" << lyman_beta_typical_g_factor << std::endl;


  double dist = 30*rMars;
  obs.fake(dist,30,300);
  observation obs_nointerp = obs;

#pragma omp parallel
  {
    observation local_obs = obs;
    observation local_obs_nointerp = obs_nointerp;
    spherical_azimuthally_symmetric_grid local_grid = grid;

    
    local_grid.brightness(local_obs);
    local_grid.brightness_nointerp(local_obs_nointerp);

    #pragma omp critical
    {
      obs += local_obs;
      obs_nointerp += local_obs_nointerp;
    }
  }

  obs.save_brightness("test/test_brightness.dat");
  obs_nointerp.save_brightness("test/test_brightness_nointerp.dat");

  
  return 0; 
}

