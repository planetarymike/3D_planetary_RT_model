//emission_gpu.cu --- routines to copy emission object to gpu

#include "Real.hpp"
#include "emission.hpp"
#include "helper_cuda.h"

template <int N_VOXELS>
void emission<N_VOXELS>::vector_to_device(voxel_vector<N_VOXELS> & device_vec,
					  voxel_vector<N_VOXELS> & host_vec,
					  bool transfer/*=true*/) {
  //if transfer = false vector is allocated on device but not copied

  host_vec.to_device(transfer);
  //point the device pointer at the same location we just moved memory to
  checkCudaErrors(
		  cudaMemcpy(&device_vec.vec,
			     &host_vec.d_vec,
			     sizeof(Real*),
			     cudaMemcpyHostToDevice)
		  );
}

template <int N_VOXELS>
void emission<N_VOXELS>::matrix_to_device(voxel_matrix<N_VOXELS> & device_mat,
					  voxel_matrix<N_VOXELS> & host_mat,
					  bool transfer/*=true*/) {
  //if transfer = false vector is allocated on device but not copied

  host_mat.to_device(transfer);
  //point the device pointer at the same location we just moved memory to
  checkCudaErrors(
		  cudaMemcpy(&device_mat.mat,
			     &host_mat.d_mat,
			     sizeof(Real*),
			     cudaMemcpyHostToDevice)
		  );
}

template <int N_VOXELS>
void emission<N_VOXELS>::copy_to_device_influence(emission<N_VOXELS> *device_emission) {
  vector_to_device(device_emission->species_T, species_T);

  vector_to_device(device_emission->dtau_species, dtau_species);
  vector_to_device(device_emission->dtau_absorber, dtau_absorber);
  vector_to_device(device_emission->abs, abs);

  matrix_to_device(device_emission->influence_matrix, influence_matrix, false);

  vector_to_device(device_emission->tau_species_single_scattering, tau_species_single_scattering, false);
  vector_to_device(device_emission->tau_absorber_single_scattering, tau_absorber_single_scattering, false);
  vector_to_device(device_emission->singlescat, singlescat, false);
  vector_to_device(device_emission->sourcefn, sourcefn, false);
  //  vector_to_device(device_emission->log_sourcefn, log_sourcefn, false);
}

template <int N_VOXELS>
void emission<N_VOXELS>::copy_to_device_brightness(emission<N_VOXELS> *device_emission) {
  //free some of the influnce stuff if we used it
  species_T.free_d_vec();
  dtau_species.free_d_vec();
  dtau_absorber.free_d_vec();
  abs.free_d_vec();
  influence_matrix.free_dmat();
  tau_species_single_scattering.free_d_vec();
  tau_absorber_single_scattering.free_d_vec();
  singlescat.free_d_vec();

  //  vector_to_device(device_emission->dtau_species, dtau_species);
  //  vector_to_device(device_emission->log_dtau_species, log_dtau_species);
  vector_to_device(device_emission->dtau_species_pt, dtau_species_pt);
  //  vector_to_device(device_emission->log_dtau_species_pt, log_dtau_species_pt);
  vector_to_device(device_emission->species_T_pt, species_T_pt);  


  //  vector_to_device(device_emission->dtau_absorber, dtau_absorber);
  //  vector_to_device(device_emission->log_dtau_absorber, log_dtau_absorber);
  vector_to_device(device_emission->dtau_absorber_pt, dtau_absorber_pt);
  //  vector_to_device(device_emission->log_dtau_absorber_pt, log_dtau_absorber_pt);
  vector_to_device(device_emission->abs_pt, abs_pt);


  vector_to_device(device_emission->sourcefn, sourcefn); 
  //  vector_to_device(device_emission->log_sourcefn, log_sourcefn); 
}

#ifdef __CUDACC__
template <int N_VOXELS>
__device__
void emission<N_VOXELS>::copy_to_shared_brightness() {
  //move important parameters to shared memory
  dtau_species.to_shared();
  dtau_species_pt.to_shared();
  species_T_pt.to_shared();
  abs_pt.to_shared();
  sourcefn.to_shared();
}

template <int N_VOXELS>
__device__
void emission<N_VOXELS>::from_shared_brightness() {
  //move important parameters to shared memory
  dtau_species.from_shared();
  dtau_species_pt.from_shared();
  species_T_pt.from_shared();
  abs_pt.from_shared();
  sourcefn.from_shared();
}
#endif


template <int N_VOXELS>
void emission<N_VOXELS>::copy_influence_to_host() {
  influence_matrix.to_host();
  
  tau_species_single_scattering.to_host();
  tau_absorber_single_scattering.to_host();
  singlescat.to_host();
}

template <int N_VOXELS> void emission<N_VOXELS>::copy_solved_to_host() {

  tau_species_single_scattering.to_host();
  tau_absorber_single_scattering.to_host();
  singlescat.to_host();

  sourcefn.to_host();
  //  log_sourcefn.to_host();
}

