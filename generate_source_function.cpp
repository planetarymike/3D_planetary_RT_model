//generate_source_function.cpp -- program to generate a source
//function for comparison with analytic solutions and other models

#include "atmosphere.h"
#include "RT_plane_parallel.h"
#include "RT_spherical_azimuthally_symmetric.h"


int main(int argc, char* argv[]) {

  //define the physical atmosphere
  double exobase_temp = 200;//K
  krasnopolsky_temperature temp(exobase_temp);

  double H_exobase_density = 1e6;// cm-3
  double CO2_exobase_density = 2e8;//cm-3
  chamb_diff_1d atm(H_exobase_density,
		    CO2_exobase_density,
		    temp);

  //use holstein functions to compute influence integrals
  holstein_approx hol;

  //define the RT grid
  int n_emissions = 1;
  // plane_parallel_grid grid(n_emissions, hol);
  // grid.setup_voxels(160, atm);
  // grid.setup_rays(10);
  spherical_azimuthally_symmetric_grid grid(n_emissions, hol);

  //grid.save_intersections = true;

  grid.rmethod = grid.rmethod_altitude;
  //grid.rmethod = grid.rmethod_lognH;

  //grid.szamethod = grid.szamethod_uniform;
  grid.szamethod = grid.szamethod_uniform_cos;
  
  
  grid.setup_voxels(20, 20/*20 for 10 deg increments*/, atm);
  grid.setup_rays(6, 12);

  
  //solve for H lyman alpha
  grid.define_emission(0,
		       "H Lyman alpha",
		       atm,
		       &chamb_diff_1d::nH,   &chamb_diff_1d::sH_lya,
		       &chamb_diff_1d::nCO2, &chamb_diff_1d::sCO2_lya);

  //solve for the source function
  grid.generate_S();
  
  //now print out the output
  grid.save_S("test/test_source_function.dat");

  //get some fake coordinates
  observation obs(n_emissions, grid.emission_names);
  obs.fake(30,600);
  grid.brightness(obs);
  obs.save_brightness("test/test_brightness.dat");
  grid.interp_brightness(obs);
  obs.save_brightness("test/test_interp_brightness.dat");

  
  


  
  return 0; 
}

