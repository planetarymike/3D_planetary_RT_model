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

  int &n_emissions;
  vector<string> emission_names;
  vector<double> emission_g_factors;
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

  vector<vector<double>> brightness;

  observation(int &n_emissionss, vector<string> &emission_namess, vector<double> &emission_g_factorss)
    : n_emissions(n_emissionss), emission_names(emission_namess), emission_g_factors(emission_g_factorss), n_obs(0)
  { }
  
  void resize(int n_obss) {
    location_MSO.resize(n_obs,vector<double>(3,0.));
    location_model.resize(n_obs,vector<double>(3,0.));
    direction_MSO.resize(n_obs,vector<double>(3,0.));
    direction_model.resize(n_obs,vector<double>(3,0.));
    
    brightness.resize(n_obs,vector<double>(n_emissions,0.));
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

  void fake(vector<double> loc = {0.,-30*rMars,0.},
	    double angle_deg = 30,
	    int nsamples = 300) {

    double angle_rad = M_PI/180. * angle_deg;
    double dangle_rad = 2*angle_rad/(nsamples-1);
    
    n_obs = nsamples*nsamples;
    
    location_model.resize(n_obs,loc);

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
  }


};








#endif
