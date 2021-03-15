//emission.h -- structure defining an atmospheric emission

#ifndef __emission_h_
#define __emission_h_

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "los_tracker.hpp"

using std::string;
using std::isnan;

template<typename emission_type,
	 template<bool, int> class los_tracker_type>
struct emission {
// base class with basic methods that other classes point to. Storage
// for RT parameters is defined in derived template classes.
protected:
  typedef emission<emission_type, los_tracker_type> this_emission_type;
  
  static const int name_length = 100;
  char internal_name[name_length];
  bool internal_init;
  bool internal_solved;

public:
  //pointer to device copy of this object
  emission_type* device_emission = NULL;

  CUDA_CALLABLE_MEMBER
  emission()
    : internal_init(false), internal_solved(false)
  {
    internal_name[0] = '\0';
  }

  CUDA_CALLABLE_MEMBER
  ~emission() {
#if defined(__CUDACC__) and not defined(__CUDA_ARCH__)
    device_clear();
#endif
  }

  
  // protected member read operations
  CUDA_CALLABLE_MEMBER
  const char* name() const {
    return internal_name;
  }

  CUDA_CALLABLE_MEMBER
  bool init() const {
    return internal_init;
  }

  CUDA_CALLABLE_MEMBER
  bool solved() const {
    return internal_solved;
  }

  // radiative transfer methods
  // override all functions with static_cast to define a new emission type
  template <bool transmission, int n_voxels>
  using los = los_tracker_type<transmission, n_voxels>;

  typedef los<false, 0> brightness_tracker;

  template<int n_voxels>
  using influence_tracker = los<true, n_voxels>;
  
  // TODO: change these to accept stepper instead of voxel number?
  // will require use of std::any, possible performance hit
  template<bool influence, int n_voxels>
  CUDA_CALLABLE_MEMBER
  void reset_tracker(const int &start_voxel,
		     los<influence,n_voxels> &tracker) const {
    static_cast<const emission_type*>(this)->reset_tracker(start_voxel,
							   tracker);
  }
  
  // update the tracker on voxel entry
  template<bool influence, int n_voxels>
  CUDA_CALLABLE_MEMBER
  void update_tracker_start(const int &current_voxel,
			    const Real & pathlength,
			    los<influence,n_voxels> &tracker) const {
    static_cast<const emission_type*>(this)->update_tracker_start(current_voxel,
							    pathlength,
							    tracker);
  }
  
  // update the tracker using an interpolation point inside a voxel 
  CUDA_CALLABLE_MEMBER
  void update_tracker_start_interp(const int &n_interp_points,
				   const int *indices,
				   const Real *weights,
				   const Real &pathlength,
				   brightness_tracker &tracker) const {
    static_cast<const emission_type*>(this)->update_tracker_start(n_interp_points,
								  indices,
								  weights,
								  pathlength,
								  tracker);
  }
  
  template<bool influence, int n_voxels>
  CUDA_CALLABLE_MEMBER
  void update_tracker_end(los<influence,n_voxels> &tracker) const {
    tracker.update_end();
  }
  
  // update the influence tracker with the contribution from this voxel
  template<int n_voxels>
  CUDA_CALLABLE_MEMBER
  void update_tracker_influence(const int &current_voxel,
				const Real &pathlength,
				const Real &domega,
				influence_tracker<n_voxels> &tracker) const {
    static_cast<const emission_type*>(this)->update_tracker_influence(current_voxel,
								      pathlength,
								      domega,
								      tracker);
  }
  
  // compute the single scattering from the input tracker
  template<bool influence, int n_voxels>
  CUDA_CALLABLE_MEMBER
  void compute_single_scattering(const int &start_voxel,
				 los<influence,n_voxels> &tracker) {
    static_cast<emission_type*>(this)->compute_single_scattering(start_voxel,
								 tracker);
  }
  
  template<int n_voxels>
  CUDA_CALLABLE_MEMBER
  void accumulate_influence(const int &start_voxel,
			    influence_tracker<n_voxels> &tracker) {
    static_cast<emission_type*>(this)->compute_single_scattering(start_voxel,
								 tracker);
  }

  CUDA_CALLABLE_MEMBER
  void reset_solution() {
    static_cast<emission_type*>(this)->reset_solution();
    internal_solved=false;
  }
  
  void solve() {
    static_cast<emission_type*>(this)->solve();
  }
  
  // update the brightness using an interpolation point inside a voxel
  CUDA_CALLABLE_MEMBER
  void update_tracker_brightness_interp(const int n_interp_points,
					const int *indices,
					const Real *weights,
					const Real &pathlength,
					brightness_tracker &tracker) const {
    static_cast<const emission_type*>(this)->update_tracker_brightness_interp(n_interp_points,
									      indices,
									      weights,
									      pathlength,
									      tracker);
  }
  
  // update the brightness with the contribution from this voxel
  CUDA_CALLABLE_MEMBER
  void update_tracker_brightness_nointerp(const int &current_voxel,
					  const Real &pathlength,
					  brightness_tracker &tracker) const {
    static_cast<const emission_type*>(this)->update_tracker_brightness_nointerp(current_voxel,
										pathlength,
										tracker);
  }

  // save routines
  void save(std::ostream &file, VectorX (*function)(VectorX, int), int i) const {
    static_cast<const emission_type*>(this)->save(file, function, i);
  }
  void save_influence(std::ostream &file) const {
    static_cast<const emission_type*>(this)->save_influence(file);
  }

#ifdef __CUDACC__
protected:
  template <typename T>
  void copy_trivial_member_to_device(T& host_val, T& device_loc, int size=1) {
    checkCudaErrors(
		    cudaMemcpy(&device_loc,
			       &host_val,
			       size*sizeof(T),
			       cudaMemcpyHostToDevice)
		    );
    checkCudaErrors( cudaPeekAtLastError() );
    checkCudaErrors( cudaDeviceSynchronize() );
  }
  template <typename T>
  void copy_trivial_member_to_host(T& host_val, T& device_loc, int size=1) {
    checkCudaErrors(
		    cudaMemcpy(&host_val,
			       &device_loc,
			       size*sizeof(T),
			       cudaMemcpyDeviceToHost)
		    );
    checkCudaErrors( cudaPeekAtLastError() );
    checkCudaErrors( cudaDeviceSynchronize() );
  }

public:
  // routines to copy emission to/from device
  // override these in derived classes to copy non-trivial members
  void allocate_device_emission() {
    if (device_emission == NULL) {
      //move grid to GPU
      checkCudaErrors(
		      cudaMalloc((void **)&device_emission,
				 sizeof(emission_type))
		      );
      checkCudaErrors( cudaPeekAtLastError() );
      checkCudaErrors( cudaDeviceSynchronize() );
    }
  }

  void copy_trivial_members_to_device() {
    this_emission_type* typed_device_emission = static_cast<this_emission_type*>(device_emission);
    
    copy_trivial_member_to_device(internal_name[0], typed_device_emission->internal_name[0], name_length);
    copy_trivial_member_to_device(internal_init, typed_device_emission->internal_init);
    copy_trivial_member_to_device(internal_solved, typed_device_emission->internal_solved);
  }
  void copy_to_device_influence() {
    copy_trivial_members_to_device();
    static_cast<emission_type*>(this)->copy_to_device_influence();
  }
  void copy_to_device_brightness() {
    copy_trivial_members_to_device();
    static_cast<emission_type*>(this)->copy_to_device_brightness();
  }

  void copy_trivial_members_to_host() {
    this_emission_type* typed_device_emission = static_cast<this_emission_type*>(device_emission);

    copy_trivial_member_to_host(internal_name[0], typed_device_emission->internal_name[0], name_length);
    copy_trivial_member_to_host(internal_init, typed_device_emission->internal_init);
    copy_trivial_member_to_host(internal_solved, typed_device_emission->internal_solved);
  }
  void copy_influence_to_host() {
    copy_trivial_members_to_host();
    static_cast<emission_type*>(this)->copy_influence_to_host();
  }
  void copy_solved_to_host() {
    copy_trivial_members_to_host();
    static_cast<emission_type*>(this)->copy_solved_to_host();
  }

  void device_clear() {
    // do not call derived class methods here, this is used in the
    // class destructor
    if (device_emission != NULL)
      checkCudaErrors(cudaFree(device_emission));
    device_emission=NULL;
    checkCudaErrors( cudaPeekAtLastError() );
    checkCudaErrors( cudaDeviceSynchronize() );
  }
#endif
};

#endif
