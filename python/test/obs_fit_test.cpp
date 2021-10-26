// obs_fit_test.cpp --- test the observation_fit class outside of the python wrapper

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "constants.hpp"
#include "observation_fit.hpp"

int main(__attribute__((unused)) int argc, __attribute__((unused)) char* argv[]) {

  observation_fit obsfit;//the object we're testing


  // grid info

  static const int n_radial_boundaries = 40;
  static const int n_sza_boundaries = 20;/*20 for 10 deg increments with szamethod_uniform*/
  static const int n_rays_theta = 6;
  static const int n_rays_phi = 12;
  typedef spherical_azimuthally_symmetric_grid<n_radial_boundaries,
  					       n_sza_boundaries,
  					       n_rays_theta,
					       n_rays_phi> grid_type;
  
  grid_type grid;

  // emission info
  static const int n_emissions = 2;

  //fake H emissions
  typedef singlet_CFR<grid_type::n_voxels> emission_type;
  emission_type lyman_alpha, lyman_beta;
  emission_type *emissions[n_emissions] = {&lyman_alpha, &lyman_beta};
  
  observation<emission_type, n_emissions> obs(emissions);
  Real dist = 30*rMars;
  Real angle = 30;
  int obs_size = 600;
  Vector3 loc = {0.,-1.0,0.};
  obs.fake(dist, angle, obs_size, loc);

  // get the geometry out of the observation object
  vector<vector<Real>> locations, directions;
  locations.resize(obs_size*obs_size);
  directions.resize(obs_size*obs_size);
  for (int i=0; i<obs_size*obs_size; i++) {
    locations[i].resize(3);
    directions[i].resize(3);
    atmo_vector vec = obs.get_vec(i);
    locations[i] = {vec.pt.x, vec.pt.y, vec.pt.z};
    directions[i] = {vec.line_x, vec.line_y, vec.line_z};
  }
  obsfit.add_observation(locations,directions);
  vector<Real> g = {lyman_alpha_typical_g_factor, lyman_beta_typical_g_factor};
  obsfit.set_g_factor(g);

  vector<Real> ra;
  vector<Real> dec;
  ra.resize(obs_size*obs_size);
  dec.resize(obs_size*obs_size);
  for (int i=0; i<obs_size*obs_size; i++) {
    ra[i] = 360.0*(i%obs_size)/(1.0*obs_size);
    dec[i] = 180.0*(i/obs_size)/(1.0*obs_size)-90;
  }

  vector<Real> marspos = {std::sqrt(2.0f),0,0};

  //obsfit.add_observation_ra_dec(marspos, ra, dec);

  // obsfit.generate_source_function_temp_asym(1e5,300,100);
  // brightness = obsfit.brightness();

  obsfit.generate_source_function(5e5,200);
  obsfit.save_influence_matrix("/home/mike/Documents/Mars/3D_planetary_RT_model/test/py_influence_matrix_cpu.dat");
  
  vector<vector<Real>> brightness;
  obsfit.O_1026_generate_source_function(2e7,200,1.69e-3);
  brightness = obsfit.O_1026_brightness();

  return 0; 
}

