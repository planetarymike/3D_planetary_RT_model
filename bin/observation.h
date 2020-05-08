//observation.h -- structures for storing atmosphere observations

#ifndef __observation_h
#define __observation_h

#include "atmo_vec.h"
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <iostream>

using std::string;
typedef Eigen::Matrix<Real, 3, 1> Vector3;
typedef Eigen::Matrix<Real,3,3> Matrix3;
typedef Eigen::AngleAxis<Real> AngleAxis;

struct observation {
protected:
  const int n_emissions;
  const vector<string> emission_names;
  int n_obs;

  // there are two coordinate systems,
  // MSO, and model
  // model_x = MSO_z
  // model_y = -MSO_y
  // model_z = MSO_x

  vector<vector<Real>> location_MSO;
  vector<vector<Real>> location_model;

  vector<vector<Real>> direction_MSO;
  vector<vector<Real>> direction_model;

  vector<atmo_vector> obs_vecs;
  
  void add_MSO_observation(vector<Real> loc, vector<Real> dir) {
    assert(loc.size() == 3 && "only three-vectors allowed");
    assert(dir.size() == 3 && "only three-vectors allowed");

    
    vector<Real> loc_MSO = { loc[0] , loc[1], loc[2] };
    location_MSO.push_back(loc_MSO);
    vector<Real> loc_model = { loc[2] , -loc[1], loc[0] };
    location_model.push_back(loc_model);

    vector<Real> dir_MSO = { dir[0] , dir[1], dir[2] };
    direction_MSO.push_back(dir_MSO);
    vector<Real> dir_model = { dir[2] , -dir[1], dir[0] };
    direction_model.push_back(dir_model);

    atmo_point pt;
    pt.xyz(loc_model[0],loc_model[1],loc_model[2]);
    obs_vecs.push_back(atmo_vector(pt, dir_model[0], dir_model[1], dir_model[2]));

    n_obs++;
  }

  void resize_input(int n_obss) {
    n_obs=n_obss;
    location_MSO.resize(n_obss);
    location_model.resize(n_obss);
    direction_MSO.resize(n_obss);
    direction_model.resize(n_obss);
    obs_vecs.resize(n_obss);
  }

public:
  vector<Real> emission_g_factors;
  
  vector<vector<Real>> brightness;
  vector<vector<Real>> tau_species;
  vector<vector<Real>> tau_absorber;

  observation(const vector<string> &emission_namess)
    : n_emissions(emission_namess.size()),
      emission_names(emission_namess),
      n_obs(0),
      emission_g_factors(vector<Real>(n_emissions,0.))
  { }

  observation& operator+=(const observation &O) {
    assert(O.emission_names == emission_names &&
	   O.location_MSO == location_MSO &&
	   O.direction_MSO == direction_MSO &&
	   O.emission_g_factors == emission_g_factors &&
	   "only add observations of the same kind");

    for (int iobs=0;iobs<n_obs;iobs++) {
      for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	brightness[iobs][i_emission] = brightness[iobs][i_emission] + O.brightness[iobs][i_emission];
	tau_species[iobs][i_emission] = tau_species[iobs][i_emission] + O.tau_species[iobs][i_emission];
	tau_absorber[iobs][i_emission] = tau_absorber[iobs][i_emission] + O.tau_absorber[iobs][i_emission];
      }
    }
    
    return *this;
  }
  
  void reset_output() {
    brightness.resize(n_obs,vector<Real>(n_emissions,0.));
    tau_species.resize(n_obs,vector<Real>(n_emissions,0.));
    tau_absorber.resize(n_obs,vector<Real>(n_emissions,0.));
  }

  void add_MSO_observation(vector<vector<Real>> &locations, vector<vector<Real>> &directions) {
    assert(locations.size() == directions.size() && "location and look direction must have the same length.");

    resize_input(0);
    
    for (unsigned int i=0;i<locations.size();i++)
      add_MSO_observation(locations[i],directions[i]);

    reset_output();
  }


  // void add_MSO_observation(Real *location_array, Real *direction_array, int n) {
  //   resize_input(n);

  //   for (int i=0;i<n;i++) {
  //     vector<Real> loc_MSO = {location_array[3*i],
  // 				location_array[3*i+1],
  // 				location_array[3*i+2]};
  //     vector<Real> dir_MSO = {direction_array[3*i],
  // 				direction_array[3*i+1],
  // 				direction_array[3*i+2]};
  //     add_MSO_observation(loc_MSO,dir_MSO);
  //   }

  //   reset_output();
  // }


  
  int size() const {
    return n_obs;
  }
  
  atmo_vector get_vec(int i) const {
    return obs_vecs[i];
  }

  vector<atmo_vector> get_vecs() const {
    return obs_vecs;
  }

  void save_brightness(string fname) {
    std::ofstream file(fname.c_str());
    if (file.is_open())
      {
	for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	  VectorX brightness_write_out;
	  brightness_write_out.resize(brightness.size());
	  for (unsigned int i=0;i<brightness.size();i++)
	    brightness_write_out[i] = brightness[i][i_emission];

	  file << emission_names[i_emission] << " brightness [kR]: " << brightness_write_out.transpose() << "\n";
	}
      }
  }

  void fake(Real dist,
	    Real angle_deg = 30,
	    int nsamples = 300) {

    vector<Real> loc = {0.,-dist,0.};

    Real angle_rad = M_PI/180. * angle_deg;
    Real dangle_rad = 2*angle_rad/(nsamples-1);
    
    n_obs = nsamples*nsamples;
    
    location_model.resize(n_obs, loc);

    atmo_point pt;
    pt.xyz(loc[0],loc[1],loc[2]);

    direction_model.resize(n_obs);
    obs_vecs.resize(n_obs);
    
    for (int i=0;i<nsamples;i++) {
      for (int j=0;j<nsamples;j++) {
	Matrix3 r;
	r = (AngleAxis(-angle_rad+i*dangle_rad, Vector3::UnitZ())
	     * AngleAxis(-angle_rad+j*dangle_rad,  Vector3::UnitX()));
	Vector3 dir = r * Vector3::UnitY();

	int iobs = i*nsamples+j;
	
	direction_model[iobs] = {dir[0], dir[1], dir[2]};
	obs_vecs[iobs] = atmo_vector(pt, dir[0], dir[1], dir[2]);
      }
    }

    reset_output();
  }

};








#endif
