//observation.h -- structures for storing atmosphere observations

#ifndef __observation_h
#define __observation_h

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "atmo_vec.hpp"
#include "los_tracker.hpp"
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

  atmo_vector *obs_vecs = NULL;
  
  void add_MSO_observation(const Vector3 &loc, const Vector3 &dir) {
    // there are two coordinate systems,
    // MSO, and model
    // model_x = MSO_z
    // model_y = -MSO_y
    // model_z = MSO_x
    
    assert(loc.size() == 3 && "only three-vectors allowed");
    assert(dir.size() == 3 && "only three-vectors allowed");
    
    //Vector3 loc_MSO = { loc[0] , loc[1], loc[2] };
    Vector3 loc_model = { loc[2] , -loc[1], loc[0] };

    //Vector3 dir_MSO = { dir[0] , dir[1], dir[2] };
    Vector3 dir_model = { dir[2] , -dir[1], dir[0] };

    atmo_point pt;
    pt.xyz(loc_model[0],loc_model[1],loc_model[2]);
    obs_vecs[n_obs] = atmo_vector(pt, dir_model[0], dir_model[1], dir_model[2]);

    n_obs++;
  }

  void resize_input(int n_obss) {
    n_obs=n_obss;

    if (obs_vecs != NULL)
      delete [] obs_vecs;
    obs_vecs = new atmo_vector[n_obs];
  }


public:
  brightness_tracker<n_emissions> *los = NULL;

  Real emission_g_factors[n_emissions];

  //CUDA device pointers for dynamic arrays
  atmo_vector *d_obs_vecs = NULL;
  brightness_tracker<n_emissions> *d_los = NULL;
  observation<n_emissions> *d_obs=NULL;
  
  observation(const std::string (&emission_namess)[n_emissions])
  : n_obs(0)
  {
    for (int i_emission=0;i_emission<n_emissions;i_emission++) {
      emission_names[i_emission] = emission_namess[i_emission];
      emission_g_factors[i_emission] = 0;
    }
  }
  ~observation() {
    if (obs_vecs != NULL)
      delete [] obs_vecs;
    if (los != NULL)
      delete [] los;

#ifdef __CUDACC__
    if (d_los != NULL)
      checkCudaErrors(cudaFree(d_los));
    if (d_obs_vecs != NULL)
      checkCudaErrors(cudaFree(d_obs_vecs));
    if (d_obs != NULL)
      checkCudaErrors(cudaFree(d_obs));
#endif
    
  }
  observation(const observation &copy) {
    assert(n_emissions == copy.n_emissions);
    for (int i = 0; i<n_emissions; i++) {
      assert(emission_names[i] == copy.emission_names[i]);
      emission_g_factors[i] = copy.emission_g_factors[i];
    }
    n_obs = copy.n_obs;
    resize_input(n_obs);
    reset_output();
    for (int i = 0; i < n_obs; i++) {
      obs_vecs[i] = copy.obs_vecs[i];
      los[i] = copy.los[i];
    }

    d_obs = copy.d_obs;
    d_obs_vecs = copy.d_obs_vecs;
    d_los = copy.d_los;
  }
  observation& operator=(const observation &rhs) {
    if(this == &rhs) return *this;

    assert(n_emissions == rhs.n_emissions);
    for (int i = 0; i<n_emissions; i++) {
      assert(emission_names[i] == rhs.emission_names[i]);
      emission_g_factors[i] = rhs.emission_g_factors[i];
    }
    
    n_obs = rhs.n_obs;
    resize_input(n_obs);
    reset_output();
    for (int i = 0; i < n_obs; i++) {
      obs_vecs[i] = rhs.obs_vecs[i];
      los[i] = rhs.los[i];
    }

    d_obs = rhs.d_obs;
    d_obs_vecs = rhs.d_obs_vecs;
    d_los = rhs.d_los;

    return *this;
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
    if (los != NULL)
      delete [] los;
    los = new brightness_tracker<n_emissions>[n_obs];
  }

  void add_MSO_observation(vector<Vector3> &locations, vector<Vector3> &directions) {
    assert(locations.size() == directions.size() && "location and look direction must have the same length.");

    resize_input(0);
    
    for (unsigned int i=0;i<locations.size();i++)
      add_MSO_observation(locations[i],directions[i]);

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

    Real angle_rad = M_PI/180. * angle_deg;
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
	obs_vecs[iobs] = atmo_vector(pt, dir[0], dir[1], dir[2]);
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
    if (d_obs_vecs != NULL)
      checkCudaErrors(cudaFree(d_obs_vecs));
    checkCudaErrors(
		    cudaMalloc((void**) &d_obs_vecs,
			       n_obs*sizeof(atmo_vector))
		    );
    checkCudaErrors(
		    cudaMemcpy(d_obs_vecs,
			       obs_vecs,
			       n_obs*sizeof(atmo_vector),
			       cudaMemcpyHostToDevice)
		    );
    //point the device pointer at the same location we just moved memory to
    checkCudaErrors(
		    cudaMemcpy(&(d_obs->obs_vecs),
			       &d_obs_vecs,
			       sizeof(atmo_vector*),
			       cudaMemcpyHostToDevice)
		    );
    
    //prepare brightness object
    if (d_los != NULL)
      checkCudaErrors(cudaFree(d_los));
    checkCudaErrors(
		    cudaMalloc((void**) &d_los,
			       n_obs*sizeof(brightness_tracker<n_emissions>) )
		    );
    //point the device pointer at the same location we just moved memory to
    checkCudaErrors(
		    cudaMemcpy(&(d_obs->los),
			       &d_los,
			       sizeof(brightness_tracker<n_emissions>*),
			       cudaMemcpyHostToDevice)
		    );
  }

  void to_host() {
    reset_output();
    checkCudaErrors(
		    cudaMemcpy(los,
			       d_los,
			       n_obs*sizeof(brightness_tracker<n_emissions>),
			       cudaMemcpyDeviceToHost)
		    );
  }

#endif




};








#endif
