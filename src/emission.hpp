//emission.h -- structure defining an atmospheric emission

#ifndef __emission_h_
#define __emission_h_

#include "Real.hpp"
#include <string>
#include <Eigen/Dense>
#include "atmo_vec.hpp"
#include "cuda_compatibility.hpp"
#include "voxel_vector.hpp"

using std::string;

template <int N_VOXELS>
struct emission {
  static const unsigned int n_voxels = N_VOXELS;
  string name;
  bool init;
  bool solved;
  
  Real branching_ratio;
  
  //these store physical atmospheric parameters on the grid (dimension n_voxels)
  //dynamic arrays are required for dim>~32 for Eigen
  //the vectors point to the Eigen objects so that these can be used interchangably
  voxel_vector species_density; //densities of scatterers and absorbers on the tabulated grid
  voxel_vector absorber_density; 
  voxel_vector species_sigma;//scatterer and absorber cross section on the tabulated grid
  voxel_vector absorber_sigma;

  voxel_vector dtau_species;
  voxel_vector log_dtau_species; //quantities that need to be interpolated are also stored as log
  voxel_vector dtau_absorber;
  voxel_vector log_dtau_absorber;
  voxel_vector abs; //ratio of dtau_abs to dtau_species
  voxel_vector log_abs; 

  //Radiative transfer parameters
  voxel_matrix influence_matrix; //influence matrix has dimensions n_voxels, n_voxels)

  //vectors to compute the single scattering have dimensions (n_voxels)  
  voxel_vector tau_species_single_scattering;
  voxel_vector tau_absorber_single_scattering;
  voxel_vector singlescat; 

  voxel_vector sourcefn; //have dimensions (n_voxels)  
  voxel_vector log_sourcefn; 

  CUDA_CALLABLE_MEMBER
  emission() :
    init(false)
  { }
  CUDA_CALLABLE_MEMBER
  ~emission() { };

  void resize() {
    species_density.resize(n_voxels);
    absorber_density.resize(n_voxels);
    species_sigma.resize(n_voxels);
    absorber_sigma.resize(n_voxels);

    dtau_species.resize(n_voxels);
    log_dtau_species.resize(n_voxels);
    dtau_absorber.resize(n_voxels);
    log_dtau_absorber.resize(n_voxels);
    abs.resize(n_voxels);
    log_abs.resize(n_voxels);

    influence_matrix.resize(n_voxels);

    singlescat.resize(n_voxels);
    sourcefn.resize(n_voxels);
    log_sourcefn.resize(n_voxels);

    tau_species_single_scattering.resize(n_voxels);
    tau_absorber_single_scattering.resize(n_voxels);
  }
  
  template<typename C, typename V>
  void define(Real emission_branching_ratio,
	      C &atmosphere,
	      Real (C::*species_density_function)(const atmo_point),
	      Real (C::*species_sigma_function)(const atmo_point),
	      Real (C::*absorber_density_function)(const atmo_point),
	      Real (C::*absorber_sigma_function)(const atmo_point),
	      const V &pts) {
    
    branching_ratio = emission_branching_ratio;
    
    for (unsigned int i_voxel=0;i_voxel<n_voxels;i_voxel++) {
      species_density[i_voxel]=(atmosphere.*species_density_function)(pts[i_voxel]);
      species_sigma[i_voxel]=(atmosphere.*species_sigma_function)(pts[i_voxel]);
      absorber_density[i_voxel]=(atmosphere.*absorber_density_function)(pts[i_voxel]);
      absorber_sigma[i_voxel]=(atmosphere.*absorber_sigma_function)(pts[i_voxel]);
    }
    
    //define differential optical depths by coefficientwise multiplication
    dtau_species = species_density.array() * species_sigma.array();
    dtau_absorber = absorber_density.array() * absorber_sigma.array();
    abs = dtau_absorber.array() / dtau_species.array();
    
    for (unsigned int i=0;i<log_dtau_species.size();i++) {
      log_dtau_species(i) = dtau_species(i) == 0 ? -1e5 : log(dtau_species(i));
      log_dtau_absorber(i) = dtau_absorber(i) == 0 ? -1e5 : log(dtau_species(i));
      log_abs(i) = abs(i) == 0 ? -1e5 : log(abs(i));
    }
    
    init=true;
  }

  //methods to transfer objects to device
  void copy_to_device_influence(emission<N_VOXELS> *device_emission);
  void copy_to_device_brightness(emission<N_VOXELS> *device_emission);
  void vector_to_device(voxel_vector & device_vec, voxel_vector & host_vec, bool transfer = true);
  void matrix_to_device(voxel_matrix & device_vec, voxel_matrix & host_vec, bool transfer = true);

  void copy_influence_to_host();
  void vector_to_host(voxel_vector & host_vec);
  void matrix_to_host(voxel_matrix & host_vec);
};

#ifdef __CUDACC__
#include "emission_gpu.cu"
#endif

#endif
