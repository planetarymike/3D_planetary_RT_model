//emission_gpu.cu --- routines to copy emission object to gpu

#include "Real.hpp"
#include "emission.hpp"
#include "helper_cuda.h"

template <int N_VOXELS>
void emission<N_VOXELS>::vector_to_device(voxel_vector & device_vec, voxel_vector & host_vec, bool transfer/*=true*/) {
  //if transfer = false vector is allocated on device but not copied
  
  //allocate the host's d_vec to point at device memory
  checkCudaErrors(
		  cudaMalloc((void **) &host_vec.d_vec,
			     n_voxels*sizeof(Real))
		  );
  if (transfer)
    //copy from host to host's device pointer
    checkCudaErrors(
		    cudaMemcpy(host_vec.d_vec,
			       host_vec.vec,
			       n_voxels*sizeof(Real),
			       cudaMemcpyHostToDevice)
		  );
  //point the device pointer at the same location we just moved memory to
  checkCudaErrors(
		  cudaMemcpy(&device_vec.vec,
			     &host_vec.d_vec,
			     sizeof(Real*),
			     cudaMemcpyHostToDevice)
		  );
}

template <int N_VOXELS>
void emission<N_VOXELS>::matrix_to_device(voxel_matrix & device_mat, voxel_matrix & host_mat, bool transfer/*=true*/) {
  //if transfer = false vector is allocated on device but not copied

  //allocate the host's d_vec to point at device memory
  checkCudaErrors(
		  cudaMalloc((void **) &host_mat.d_mat,
			     n_voxels*n_voxels*sizeof(Real))
		  );
  if (transfer)
    //copy from host to host's device pointer
    checkCudaErrors(
		    cudaMemcpy(host_mat.d_mat,
			       host_mat.mat,
			       n_voxels*n_voxels*sizeof(Real),
			       cudaMemcpyHostToDevice)
		    );
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

  vector_to_device(device_emission->species_sigma, species_sigma);

  vector_to_device(device_emission->dtau_species, dtau_species);
  vector_to_device(device_emission->dtau_absorber, dtau_absorber);

  matrix_to_device(device_emission->influence_matrix, influence_matrix, false);

  vector_to_device(device_emission->tau_species_single_scattering, tau_species_single_scattering, false);
  vector_to_device(device_emission->tau_absorber_single_scattering, tau_absorber_single_scattering, false);
  vector_to_device(device_emission->singlescat, singlescat, false);
}

template <int N_VOXELS>
void emission<N_VOXELS>::copy_to_device_brightness(emission<N_VOXELS> *device_emission) {
  vector_to_device(device_emission->log_dtau_species, log_dtau_species);
  vector_to_device(device_emission->log_dtau_absorber, log_dtau_absorber);
  vector_to_device(device_emission->log_sourcefn, log_sourcefn); 
}

template <int N_VOXELS>
void emission<N_VOXELS>::vector_to_host(voxel_vector & host_vec) {

  // host_vec already has a pointer to the device location
  // (allocated in vector_to_device)
  // copy back from there
  checkCudaErrors(
		  cudaMemcpy(host_vec.vec,
			     host_vec.d_vec,
			     n_voxels*sizeof(Real),
			     cudaMemcpyDeviceToHost)
		  );
}

template <int N_VOXELS>
void emission<N_VOXELS>::matrix_to_host(voxel_matrix & host_mat) {

  // host_mat already has a pointer to the device location
  // (allocated in matrix_to_device)
  // copy back from there
  checkCudaErrors(
		  cudaMemcpy(host_mat.mat,
			     host_mat.d_mat,
			     n_voxels*n_voxels*sizeof(Real),
			     cudaMemcpyDeviceToHost)
		  );
}

template <int N_VOXELS>
void emission<N_VOXELS>::copy_influence_to_host() {
  matrix_to_host(influence_matrix);

  vector_to_host(tau_species_single_scattering);
  vector_to_host(tau_absorber_single_scattering);
  vector_to_host(singlescat);
}
