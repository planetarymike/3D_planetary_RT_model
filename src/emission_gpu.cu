//emission_gpu.cu --- routines to copy emission object to gpu

#include "Real.hpp"
#include "emission.hpp"
#include "helper_cuda.h"

template <int N_VOXELS>
void emission<N_VOXELS>::vector_to_device(voxel_vector<N_VOXELS> & device_vec,
					  voxel_vector<N_VOXELS> & host_vec,
					  const int n_dim,
					  const int *dim,
					  const bool transfer/*=true*/)
{
  //if transfer = false vector is allocated on device but not copied
  //texture memory is used when USE_CUDA_TEXTURES flag is set
  if (transfer) {
#ifdef USE_CUDA_TEXTURES
    host_vec.to_device_read_only(&device_vec, n_dim, dim);
#else
    host_vec.to_device(/*transfer = */true);
#endif
  } else
    host_vec.to_device(/*transfer = */false);
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
					  const bool transfer/*=true*/) {
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
void emission<N_VOXELS>::copy_to_device_influence(emission<N_VOXELS> *device_emission, const int n_dim, const int*dim) {

  //copy the read-only atmosphere arrays
  bool transfer = true;
  
  vector_to_device(device_emission->species_T_ratio, species_T_ratio, n_dim, dim, transfer);

  vector_to_device(device_emission->dtau_species, dtau_species, n_dim, dim, transfer);
  vector_to_device(device_emission->dtau_absorber, dtau_absorber, n_dim, dim, transfer);
  vector_to_device(device_emission->abs, abs, n_dim, dim, transfer);


  //now the arrays we populate on the device
  transfer = false;

  matrix_to_device(device_emission->influence_matrix, influence_matrix, transfer);

  vector_to_device(device_emission->tau_species_single_scattering, tau_species_single_scattering, n_dim, dim, transfer);
  vector_to_device(device_emission->tau_absorber_single_scattering, tau_absorber_single_scattering, n_dim, dim, transfer);
  vector_to_device(device_emission->singlescat, singlescat, n_dim, dim, transfer);
  vector_to_device(device_emission->sourcefn, sourcefn, n_dim, dim, transfer);
}

template <int N_VOXELS>
void emission<N_VOXELS>::copy_to_device_brightness(emission<N_VOXELS> *device_emission,
						   const int n_dim, const int *dim)
{
  //free some of the influnce stuff if we used it
  species_T_ratio.free_d_vec();
  dtau_species.free_d_vec();
  dtau_absorber.free_d_vec();
  abs.free_d_vec();
  influence_matrix.free_dmat();
  tau_species_single_scattering.free_d_vec();
  tau_absorber_single_scattering.free_d_vec();
  singlescat.free_d_vec();

  bool transfer = true;
  
  vector_to_device(device_emission->dtau_species_pt, dtau_species_pt, n_dim, dim, transfer);
  vector_to_device(device_emission->species_T_ratio_pt, species_T_ratio_pt, n_dim, dim, transfer);  

  vector_to_device(device_emission->dtau_absorber_pt, dtau_absorber_pt, n_dim, dim, transfer);
  vector_to_device(device_emission->abs_pt, abs_pt, n_dim, dim, transfer);

  vector_to_device(device_emission->sourcefn, sourcefn, n_dim, dim, transfer); 
}

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

