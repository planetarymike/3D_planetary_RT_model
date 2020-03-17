//observation.h -- structures for storing atmosphere observations

#ifndef __observation_h
#define __observation_h

#include "atmo_vec.h"
#include <Eigen/Dense>
#include <string>
#include <fstream>

using std::string;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::AngleAxisd;

struct observation {
protected:
  const int &n_emissions;
  const vector<string> emission_names;
  int n_obs;

  // there are two coordinate systems,
  // MSO, and model
  // model_x = MSO_z
  // model_y = -MSO_y
  // model_z = MSO_x
  vector<vector<double>> location_MSO;
  vector<vector<double>> location_model;

  vector<vector<double>> direction_MSO;
  vector<vector<double>> direction_model;

  vector<atmo_vector> obs_vecs;

  void add_MSO_observation(vector<double> &loc, vector<double> &dir) {
    assert(loc.size() == 3 && "only three-vectors allowed");
    assert(dir.size() == 3 && "only three-vectors allowed");
    
    location_MSO.push_back(loc);
    vector<double> loc_model = { loc[2] , -loc[1], loc[0] };
    location_model.push_back(loc_model);

    direction_MSO.push_back(dir);
    vector<double> dir_model = { dir[2] , -dir[1], dir[0] };
    direction_model.push_back(dir_model);

    atmo_point pt;
    pt.xyz(loc_model[0],loc_model[1],loc_model[2]);
    obs_vecs.push_back(atmo_vector(pt, dir_model[0], dir_model[1], dir_model[2]));

    n_obs++;
  }

  void resize_input(int n_obss) {
    n_obs=n_obss;
    location_MSO.resize(n_obss,vector<double>(3,0.));
    location_model.resize(n_obss,vector<double>(3,0.));
    direction_MSO.resize(n_obss,vector<double>(3,0.));
    direction_model.resize(n_obss,vector<double>(3,0.));
  }
  void resize_output() {
    brightness.resize(n_obs,vector<double>(n_emissions,0.));
    tau_species.resize(n_obs,vector<double>(n_emissions,0.));
    tau_absorber.resize(n_obs,vector<double>(n_emissions,0.));
  }

public:
  vector<double> emission_g_factors;
  
  vector<vector<double>> brightness;
  vector<vector<double>> tau_species;
  vector<vector<double>> tau_absorber;

  observation(const vector<string> &emission_namess)
    : n_emissions(emission_namess.size()),
      emission_names(emission_namess),
      n_obs(0),
      emission_g_factors(vector<double>(n_emissions,0.))
  { }

  void add_MSO_observation(vector<vector<double>> &locations, vector<vector<double>> &directions) {
    assert(locations.size() == directions.size() && "location and look direction must have the same length.");

    resize_input(0);
    
    for (unsigned int i=0;i<locations.size();i++)
      add_MSO_observation(locations[i],directions[i]);

    resize_output();
  }

  int size() const {
    return n_obs;
  }
  
  vector<atmo_vector> get_obs_vecs() const {
    return obs_vecs;
  }

  void save_brightness(string fname) {
    std::ofstream file(fname.c_str());
    if (file.is_open())
      {
	for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	  VectorXd brightness_write_out;
	  brightness_write_out.resize(brightness.size());
	  for (unsigned int i=0;i<brightness.size();i++)
	    brightness_write_out[i] = brightness[i][i_emission];

	  file << emission_names[i_emission] << " brightness [kR]: " << brightness_write_out.transpose() << "\n";
	}
      }
  }

  void fake(vector<double> loc = {0.,-1e10,0.},
	    double angle_deg = 30,
	    int nsamples = 300) {

    double angle_rad = M_PI/180. * angle_deg;
    double dangle_rad = 2*angle_rad/(nsamples-1);
    
    n_obs = nsamples*nsamples;
    
    location_model.resize(n_obs, loc);

    atmo_point pt;
    pt.xyz(loc[0],loc[1],loc[2]);

    direction_model.resize(n_obs);
    obs_vecs.resize(n_obs);
    
    for (int i=0;i<nsamples;i++) {
      for (int j=0;j<nsamples;j++) {
	Matrix3d r;
	r = (AngleAxisd(-angle_rad+i*dangle_rad, Vector3d::UnitZ())
	     * AngleAxisd(-angle_rad+j*dangle_rad,  Vector3d::UnitX()));
	Vector3d dir = r * Vector3d::UnitY();

	int iobs = i*nsamples+j;
	
	direction_model[iobs] = {dir[0], dir[1], dir[2]};
	obs_vecs[iobs] = atmo_vector(pt, dir[0], dir[1], dir[2]);
      }
    }

    resize_output();
  }

};








#endif
