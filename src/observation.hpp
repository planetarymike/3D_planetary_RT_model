//observation.h -- structures for storing atmosphere observations

#ifndef __observation_h
#define __observation_h

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "atmo_vec.hpp"
#include "los_tracker.hpp"
#include "gpu_vector.hpp"
#include "emission.hpp"
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

template<typename emission_type, int N_EMISS>
struct observation {
protected:
  static const int n_emissions = N_EMISS;
  emission_type* emissions[n_emissions];
  
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
    vector<Real> dir_model = { dir[2] , -dir[1], dir[0] };\

    atmo_point pt;
    pt.xyz(loc_model[0],loc_model[1],loc_model[2]);
    obs_vecs[i].ptxyz(pt, dir_model[0], dir_model[1], dir_model[2]);
  }

  void reset_output() {
    // need to resize with type belonging to each emission
    for (int i_emission = 0; i_emission < n_emissions; i_emission++) {
      los[i_emission].resize(n_obs, typename emission_type::brightness_tracker());
      for (int i_obs = 0; i_obs < n_obs; i_obs++)
        los[i_emission][i_obs].init();
    }
  }

  void resize_input(int n_obss) {
    n_obs=n_obss;
    obs_vecs.resize(n_obs);
    reset_output();
  }

public:
  gpu_vector<typename emission_type::brightness_tracker> los[n_emissions];

  Real emission_g_factors[n_emissions];

  vector<vector<Real>> iph_brightness_unextincted;
  vector<vector<Real>> iph_brightness_observed;
  vector<Real> ra;
  vector<Real> dec;
  vector<Real> mars_ecliptic_pos;

  //CUDA device pointers for dynamic arrays
  observation<emission_type, n_emissions> *d_obs=NULL;

  observation(emission_type* (&emissionss)[n_emissions])
    : n_obs(0)
  {
    for (int i_emission=0;i_emission<n_emissions;i_emission++) {
      emissions[i_emission] = emissionss[i_emission];
      emission_g_factors[i_emission] = 0;
    }
  }
  ~observation() {
#if defined(__CUDACC__) and not defined(__CUDA_ARCH__)
    device_clear();
#endif
  }
  observation(const observation<emission_type, n_emissions> &copy) = delete;
  observation<emission_type, n_emissions>& operator=(const observation<emission_type, n_emissions> & rhs) = delete;
  
  void set_emission_g_factors(Real (&g)[n_emissions]) {
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      emission_g_factors[i_emission] = g[i_emission];
  }
  void get_emission_g_factors(Real (&g)[n_emissions]) const {
    for (int i_emission=0;i_emission<n_emissions;i_emission++)
      g[i_emission] = emission_g_factors[i_emission];
  }
  
  void add_MSO_observation(const vector<vector<Real>> &locations, const vector<vector<Real>> &directions) {
    assert(locations.size() == directions.size() && "location and look direction must have the same length.");

    resize_input(locations.size());
    
    for (unsigned int i=0;i<locations.size();i++)
      add_MSO_observation(locations[i],directions[i],i);
  }

  void add_observation_ra_dec(const std::vector<Real> &mars_ecliptic_coords,
			      const std::vector<Real> &RAA,
			      const std::vector<Real> &Decc) {
    assert(n_obs == RAA.size() && RAA.size() == Decc.size() &&
	   "IPH coordinates must have the same dimensions as locations and directions");

    assert(mars_ecliptic_coords.size()==3 && "mars ecliptic coords must be a 3D position.");
    mars_ecliptic_pos.resize(3);
    for (int i=0;i<3;i++)
      mars_ecliptic_pos[i] = mars_ecliptic_coords[i];

    ra.resize(n_obs);
    dec.resize(n_obs);
    for (unsigned int i = 0; i < n_obs; i++) {
      ra[i] = RAA[i];
      dec[i] = Decc[i];
    }

    iph_brightness_unextincted.resize(n_obs);
    iph_brightness_observed.resize(n_obs);
    for (int i_obs=0;i_obs<n_obs;i_obs++) {
      iph_brightness_unextincted[i_obs].resize(n_emissions);
      iph_brightness_observed[i_obs].resize(n_emissions);
    }
  }

  void update_iph_extinction() {
    for (int i_obs=0;i_obs<n_obs;i_obs++) {
      for (int i_emission=0;i_emission<n_emissions;i_emission++) {
	if (los[i_emission][i_obs]->tau_absorber_final != -1)
	  iph_brightness_observed[i_obs][i_emission] = (iph_brightness_unextincted[i_obs][i_emission]
							*std::exp(-los[i_emission][i_obs]->tau_absorber_final));
	else
	  iph_brightness_observed[i_obs][i_emission] = 0.0;
      }
    }
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
	    brightness_write_out[i] = los[i_emission][i].brightness;

	  file << emissions[i_emission]->name() << " brightness [kR]: " << brightness_write_out.transpose() << "\n";
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
  }

#ifdef __CUDACC__
  void to_device() {
    if (d_obs == NULL) {
      //allocate space on GPU
      checkCudaErrors(
		      cudaMalloc((void **) &d_obs,
				 sizeof(*this))
		      );
    }
    //copy observation over
    checkCudaErrors(
		    cudaMemcpy(d_obs,
			       this,
			       sizeof(*this),
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

    for (int i_emission=0; i_emission<n_emissions; i_emission++) {
      //prepare simulated brightness storage
      los[i_emission].to_device();
      //point the object device pointer to the same location
      checkCudaErrors(
		      cudaMemcpy(&((d_obs->los[i_emission]).v),
				 &los[i_emission].d_v,
				 sizeof(typename emission_type::brightness_tracker*),
				 cudaMemcpyHostToDevice)
		      );
    }
  }

  void to_host() {
    for (int i_emission=0; i_emission<n_emissions; i_emission++)
      los[i_emission].to_host();
  }

  void device_clear() {
    obs_vecs.device_clear();

    for (int i_emission=0; i_emission<n_emissions; i_emission++)
      los[i_emission].device_clear();

    if (d_obs != NULL)
      checkCudaErrors(cudaFree(d_obs));
    d_obs=NULL;
  }
  
#endif

};


#endif
