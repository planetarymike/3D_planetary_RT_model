//observation.h -- structures for storing atmosphere observations

#ifndef __observation_h
#define __observation_h

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "atmo_vec.hpp"
#include "los_tracker.hpp"
#include "gpu_vector.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include <Eigen/Dense>
typedef Eigen::Matrix<Real, 3, 1> Vector3;
typedef Eigen::Matrix<Real,3,3> Matrix3;
typedef Eigen::AngleAxis<Real> AngleAxis;

template<int N_EMISS>
struct observation {
protected:
  static const int n_emissions = N_EMISS;
  std::string emission_names[n_emissions];

  int n_obs;

  gpu_vector<atmo_vector> obs_vecs;
  
  void add_MSO_observation(const vector<Real> &loc, const vector<Real> &dir, const int i) {
    // there are two coordinate systems,
    // MSO, and model
    // model_x = MSO_z
    // model_y = -MSO_y
    // model_z = MSO_x
    
    assert(loc.size() == 3 && "only three-vectors allowed");
    assert(dir.size() == 3 && "only three-vectors allowed");
    
    //vector<Real> loc_MSO = { loc[0] , loc[1], loc[2] };
    vector<Real> loc_model = { loc[2] , -loc[1], loc[0] };

    //vector<Real> dir_MSO = { dir[0] , dir[1], dir[2] };
    vector<Real> dir_model = { dir[2] , -dir[1], dir[0] };

    atmo_point pt;
    pt.xyz(loc_model[0],loc_model[1],loc_model[2]);
    obs_vecs[i].ptxyz(pt, dir_model[0], dir_model[1], dir_model[2]);
  }

  void resize_input(int n_obss) {
    n_obs=n_obss;

    obs_vecs.resize(n_obs);
    los_observed.resize(n_obs);
  }

public:

  gpu_vector<brightness_tracker<n_emissions>> los;
  gpu_vector<observed<n_emissions>> los_observed;
  Real log_likelihood;

  Real emission_g_factors[n_emissions];

  //CUDA device pointers for dynamic arrays
  observation<n_emissions> *d_obs=NULL;
  
  observation()
    : n_obs(0)
  {
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      emission_g_factors[i_emission] = 0;
  }
  observation(const std::string (&emission_names)[n_emissions])
  : observation()
  {
    set_names(emission_names);
  }
  ~observation() {
#ifdef __CUDACC__
    device_clear();
#endif
  }

  void set_names(const string (&emission_namess)[n_emissions])
  {
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      emission_names[i_emission] = emission_namess[i_emission];
  }
  
  void set_emission_g_factors(Real (&g)[n_emissions]) {
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      emission_g_factors[i_emission] = g[i_emission];
  }
  void get_emission_g_factors(Real (&g)[n_emissions]) const {
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      g[i_emission] = emission_g_factors[i_emission];
  }
  
  void reset_output() {
    los.resize(n_obs,brightness_tracker<n_emissions>());
    for (int i_obs=0;i_obs<n_obs;i_obs++)
      los[i_obs].init();
  }

  void add_MSO_observation(const vector<vector<Real>> &locations, const vector<vector<Real>> &directions) {
    assert(locations.size() == directions.size() && "location and look direction must have the same length.");

    resize_input(locations.size());
    
    for (unsigned int i=0;i<locations.size();i++)
      add_MSO_observation(locations[i],directions[i],i);

    reset_output();
  }

  CUDA_CALLABLE_MEMBER
  int size() const {
    return n_obs;
  }

  CUDA_CALLABLE_MEMBER
  atmo_vector get_vec(int i) const {
    return obs_vecs[i];
  }

  void save_brightness(string fname) {
    std::ofstream file(fname.c_str());
    if (file.is_open())
      {
	for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	  VectorX brightness_write_out;
	  brightness_write_out.resize(n_obs);
	  for (int i=0;i<n_obs;i++)
	    brightness_write_out[i] = los[i].brightness[i_emission];

	  file << emission_names[i_emission] << " brightness [kR]: " << brightness_write_out.transpose() << "\n";
	}
      }
  }

  void fake(Real dist,
	    Real angle_deg = 30,
	    int nsamples = 300) {

    Vector3 loc = {0.,-dist,0.};

    Real angle_rad = pi/180. * angle_deg;
    Real dangle_rad = 2*angle_rad/(nsamples-1);
    
    resize_input(nsamples*nsamples);
    
    atmo_point pt;
    pt.xyz(loc[0],loc[1],loc[2]);
    
    for (int i=0;i<nsamples;i++) {
      for (int j=0;j<nsamples;j++) {
	Matrix3 r;
	r = (AngleAxis(-angle_rad+i*dangle_rad, Vector3::UnitZ())
	     * AngleAxis(-angle_rad+j*dangle_rad,  Vector3::UnitX()));
	Vector3 dir = r * Vector3::UnitY();

	int iobs = i*nsamples+j;
	
	//direction_model[iobs] = {dir[0], dir[1], dir[2]};
	obs_vecs[iobs].ptxyz(pt, dir[0], dir[1], dir[2]);
      }
    }

    reset_output();
  }

#ifdef __CUDACC__
  void to_device() {
    if (d_obs == NULL) {
      //allocate space on GPU
      checkCudaErrors(
		      cudaMalloc((void **) &d_obs,
				 sizeof(observation<n_emissions>))
		      );
    }
    //copy observation over
    checkCudaErrors(
		    cudaMemcpy(d_obs,
			       this,
			       sizeof(observation<n_emissions>),
			       cudaMemcpyHostToDevice)
		    );
    //static arrays are copied automatically by CUDA, including emission_g_factors
    
    //move geometry
    obs_vecs.to_device();
    //point the object device pointer to the same location
    checkCudaErrors(
		    cudaMemcpy(&((d_obs->obs_vecs).v),
			       &obs_vecs.d_v,
			       sizeof(atmo_vector*),
			       cudaMemcpyHostToDevice)
		    );

    //move brightness and sigma
    los_observed.to_device();
    checkCudaErrors(
		    cudaMemcpy(&((d_obs->los_observed).v),
			       &los_observed.d_v,
			       sizeof(atmo_vector*),
			       cudaMemcpyHostToDevice)
		    );
    
    //prepare simulated brightness storage
    los.to_device();
    //point the object device pointer to the same location
    checkCudaErrors(
		    cudaMemcpy(&((d_obs->los).v),
			       &los.d_v,
			       sizeof(brightness_tracker<n_emissions>*),
			       cudaMemcpyHostToDevice)
		    );
  }

  void to_host() {
    los.to_host();
    // checkCudaErrors(
    // 		    cudaMemcpy(&(log_likelihood),
    // 			       &(d_obs.log_likelihood),
    // 			       sizeof(Real),
    // 			       cudaMemcpyDeviceToHost)
    // 		    );
  }

  void device_clear() {
    obs_vecs.device_clear();
    los.device_clear();
    los_observed.device_clear();
    if (d_obs != NULL)
      checkCudaErrors(cudaFree(d_obs));
    d_obs=NULL;
  }
  
#endif




};








#endif
