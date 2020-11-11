//generate_source_function.cpp -- program to generate a source
//function for comparison with analytic solutions and other models

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "constants.hpp"
#include "observation_fit.hpp"

int main(__attribute__((unused)) int argc, __attribute__((unused)) char* argv[]) {

  observation_fit obsfit;//the object we're testing

  //a fake observation to get some geometry
  static const int n_emissions=1;
  string emission_names[n_emissions] = {"fake"};
  observation<n_emissions> obs(emission_names);

  int n_fake = 600;
  Real dist = 30*rMars;
  obs.fake(dist,30,n_fake);

  vector<vector<Real>> locations, directions;
  locations.resize(n_fake*n_fake);
  directions.resize(n_fake*n_fake);
  for (int i=0; i<n_fake*n_fake; i++) {
    locations[i].resize(3);
    directions[i].resize(3);
    atmo_vector vec = obs.get_vec(i);
    locations[i] = {vec.pt.x, vec.pt.y, vec.pt.z};
    directions[i] = {vec.line_x, vec.line_y, vec.line_z};
  }
  obsfit.add_observation(locations,directions);
  vector<Real> g = {lyman_alpha_typical_g_factor, lyman_beta_typical_g_factor};
  obsfit.set_g_factor(g);

  obsfit.generate_source_function_temp_asym(1e5,300,100);

  vector<vector<Real>> brightness;
  brightness = obsfit.brightness();

  return 0; 
}

