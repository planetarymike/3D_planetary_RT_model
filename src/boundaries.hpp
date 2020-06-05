// boundaries.hpp -- structure to contain boundary crossing information for ray-grid intersections

#ifndef __boundaries_H
#define __boundaries_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "atmo_vec.hpp"
#include <algorithm>
#include <limits>
#include <iostream>
#include <stdexcept>


template <int NDIM>
struct boundary {
  int entering; //index of voxel the ray is entering
  static const int n_dimensions = NDIM;
  int entering_indices[NDIM]; //indices in each dimension of the voxel
  Real distance; //distance to the boundary crossing

  CUDA_CALLABLE_MEMBER
  boundary() {
    entering = -2;
    for (int i=0;i<n_dimensions;i++)
      entering_indices[i] = -2;
    distance = -1;
  }
  CUDA_CALLABLE_MEMBER
  ~boundary() { }
  CUDA_CALLABLE_MEMBER
  boundary(const boundary<NDIM> &copy) {
    entering = copy.entering;
    distance = copy.distance;
    for (int i=0;i<NDIM;i++)
      entering_indices[i] = copy.entering_indices[i];
  }
  CUDA_CALLABLE_MEMBER
  boundary & operator=(const boundary<NDIM> &rhs) {
    if(this == &rhs) return *this;
    entering = rhs.entering;
    distance = rhs.distance;
    for (int i=0;i<NDIM;i++)
      entering_indices[i] = rhs.entering_indices[i];

    return *this;
  }

  CUDA_CALLABLE_MEMBER
  void reset() {
    entering = -2;
    for (int i=0;i<n_dimensions;i++)
      entering_indices[i] = -2;
    distance = -1;
  }

  CUDA_CALLABLE_MEMBER
  bool operator<(const boundary<NDIM> &rhs) const { return distance < rhs.distance; }
  CUDA_CALLABLE_MEMBER
  bool operator<=(const boundary<NDIM> &rhs) const { return distance <= rhs.distance; }
};








template <int NDIM, unsigned int MAX_SIZE>
class boundary_set {
private:
  unsigned int begin;
  unsigned int internal_size;
  static const unsigned int internal_max_size = MAX_SIZE;
  boundary<NDIM> boundaries[MAX_SIZE];
public:
  static const unsigned int n_dimensions = NDIM;

  CUDA_CALLABLE_MEMBER
  boundary_set()
  { 
    internal_size=0;begin=0;
  }
  CUDA_CALLABLE_MEMBER
  ~boundary_set() { }
  CUDA_CALLABLE_MEMBER
  boundary_set(const boundary_set<NDIM,MAX_SIZE> &copy) {
    //    assert(n_dimensions == copy.n_dimensions && "only copy same dimension boundaries");

    begin = copy.begin;
    internal_size = copy.internal_size;

    for (unsigned int i=begin;i<begin+internal_size;i++)
      boundaries[i] = copy.boundaries[i];
  }
  CUDA_CALLABLE_MEMBER
  boundary_set &operator=(const boundary_set<NDIM,MAX_SIZE> &rhs) {
    //    assert(n_dimensions == rhs.n_dimensions && "only copy same dimension boundaries");
    if(this == &rhs) return *this;

    begin = rhs.begin;
    internal_size = rhs.internal_size;

    for (unsigned int i=begin;i<begin+internal_size;i++)
      boundaries[i] = rhs.boundaries[i];

    return *this;
  }

  CUDA_CALLABLE_MEMBER
  boundary<NDIM> operator[](int i) const {
    return boundaries[begin+i];
  }

  CUDA_CALLABLE_MEMBER
  boundary<NDIM> back() const {
    return boundaries[begin+internal_size-1];
  }
  
  CUDA_CALLABLE_MEMBER
  unsigned int size() const {
    return internal_size;
  }

  CUDA_CALLABLE_MEMBER
  void reset() {
    internal_size=0;begin=0;
  }
  
  CUDA_CALLABLE_MEMBER
  void sort() {
    //basic insertion sort, replace with something better when arrays get large
    boundary<NDIM> key;
    int i, j;  
    for (i = begin+1; i < (int) (begin+internal_size); i++) 
      {  
        key = boundaries[i];  
        j = i - 1;  
	
        /* Move elements of arr[0..i-1], that are  
	   greater than key, to one position ahead  
	   of their current position */
        while (j >= 0 && key < boundaries[j]) 
	  {  
            boundaries[j + 1] = boundaries[j];  
            j = j - 1;  
	  }  
        boundaries[j + 1] = key;  
      }  
  }

  CUDA_CALLABLE_MEMBER
  void append(boundary<NDIM> b) {
    assert(internal_size < internal_max_size);
    boundaries[begin+internal_size] = b;
    internal_size++;
  }

  CUDA_CALLABLE_MEMBER
  void add_intersections(const Real start, const int dim, 
			 const int idx, const Real &coordinate,
			 const Real (&distances)[2], const int n_hits) {
    boundary<NDIM> new_boundary;
    new_boundary.reset();

    assert(0 <= n_hits && n_hits <=2 && "only 0, 1, or 2 intersections are possible with a convex shell");
    
    if (n_hits == 0) {
      return;
    } else if (n_hits == 1) {
      bool above = (start > coordinate);
      

      new_boundary.entering_indices[dim] = above ? idx-1 : idx;
      new_boundary.distance=distances[0];

      append(new_boundary);
    } else { // n_hits == 2

      bool above = (start > coordinate);
      bool in_order = (distances[1] > distances[0]);

      new_boundary.entering_indices[dim] = above ? idx-1 : idx  ;
      new_boundary.distance = in_order ? distances[0] : distances[1];
      append(new_boundary);
      
      new_boundary.entering_indices[dim] = above ? idx   : idx-1;
      //^^ indices are reversed to come back out of the region we just entered
      new_boundary.distance = in_order ? distances[1] : distances[0];
      append(new_boundary);
    }
  }

  CUDA_CALLABLE_MEMBER
  void propagate_indices() {
    //propogate voxel indices in each coordinate
    for (unsigned int i=1;i<size();i++) {
      for (unsigned int j=0;j<n_dimensions;j++)
	if (boundaries[i].entering_indices[j]==-2)
	  boundaries[i].entering_indices[j] = boundaries[i-1].entering_indices[j];
    }
  }

  template <class C>
  CUDA_CALLABLE_MEMBER
  void assign_voxel_indices(C *obj) {
    for (unsigned int i=1;i<size();i++) 
      (obj->indices_to_voxel)(boundaries[i].entering_indices, boundaries[i].entering);
  }

  
  CUDA_CALLABLE_MEMBER
  void trim() {
    //trim the list of intersections to those inside the grid
    while (boundaries[begin].entering == -1 && begin < size()-1)
      begin++;

    if (begin==size()-1) {
      begin=internal_size=0;
    } else {
      unsigned int end = begin;
      do {
	end++;
      } while (boundaries[end].entering != -1 && end < size()-1);

      assert(end < size() && "end must be inside list of boundary intersections");
      assert(boundaries[end].entering == -1 && "trim error in boundary_set: ray does not exit grid");
    
      internal_size=end-begin+1;
    }
  }

  
  CUDA_CALLABLE_MEMBER
  bool check(const int (&n_bounds)[NDIM], const int &n_voxels) const {
    if (size() > 0) {
      assert(size() > 1 && "there must be more than one boundary crossing for each ray");
      
      for (unsigned int i=begin+1;i<begin+internal_size;i++) {
	int n_dims_changing = 0;
	for (unsigned int j=0;j<n_dimensions;j++) {
	  assert(boundaries[i].entering_indices[j] >= -1
		 && "entering index must be greater than -1.");
	  assert(boundaries[i].entering_indices[j] <= n_bounds[j]
		 && "entering index must be less than the number of voxels in this dimension.");
	  
	  int diff = boundaries[i].entering_indices[j] - boundaries[i-1].entering_indices[j];

	  if (diff==0)
	    continue;

	  //we only get here if this boundary index is changing
	  n_dims_changing++;
	  assert(n_dims_changing <=1 && "only one boundary can be crossed at a time.");
	  assert((diff==1||diff==-1) && "changes in each dimension must be continuous.");
	}	

	assert(boundaries[i].entering < n_voxels
	       && boundaries[i].entering >= -1
	       && "voxel index must be in bounds.");

      }
      
      assert((back().entering_indices[0] == -1 ||
	      back().entering_indices[0] == n_bounds[0])
	     && "ray must exit grid radially via the top or bottom");
    }
    return true;
  }
};





template <int NDIM, int MAX_SIZE>
struct boundary_intersection_stepper {
  bool init;
  atmo_vector vec;
  boundary_set<NDIM, MAX_SIZE> boundaries;

  bool inside;
  bool exits_bottom;
  bool exits_top;
  unsigned int i_boundary;

  int start_voxel;
  int current_voxel;
  Real pathlength;

  CUDA_CALLABLE_MEMBER
  boundary_intersection_stepper() : init(false) { }

  CUDA_CALLABLE_MEMBER
  ~boundary_intersection_stepper() { }

  CUDA_CALLABLE_MEMBER
  boundary_intersection_stepper(const boundary_intersection_stepper &copy) {
    init = copy.init;
    vec = copy.vec;
    boundaries = copy.boundaries;

    inside = copy.inside;
    exits_bottom = copy.exits_bottom;
    exits_top = copy.exits_top;
    i_boundary = copy.i_boundary;

    start_voxel = copy.start_voxel;
    current_voxel = copy.current_voxel;
    pathlength = copy.pathlength;
  }

  CUDA_CALLABLE_MEMBER
  boundary_intersection_stepper &operator=(const boundary_intersection_stepper &rhs) {
    if(this == &rhs) return *this;

    init = rhs.init;
    vec = rhs.vec;
    boundaries = rhs.boundaries;

    inside = rhs.inside;
    exits_bottom = rhs.exits_bottom;
    exits_top = rhs.exits_top;
    i_boundary = rhs.i_boundary;

    start_voxel = rhs.start_voxel;
    current_voxel = rhs.current_voxel;
    pathlength = rhs.pathlength;

    return *this;
  }

  CUDA_CALLABLE_MEMBER
  void init_stepper()
  {
    if (boundaries.size() > 0) {
      start_voxel = boundaries[0].entering;
      
      if (boundaries.back().entering_indices[0] == -1) {
	exits_bottom = true;
	exits_top = false;
      } else {
	exits_bottom = false;
	exits_top = true;
      }
    }
    init=true;
  }
  
  CUDA_CALLABLE_MEMBER
  boundary_intersection_stepper(atmo_vector vecc, boundary_set<NDIM,MAX_SIZE> boundariess)
    : vec(vecc), boundaries(boundariess)
  {
    init_stepper();
  }
  
  CUDA_CALLABLE_MEMBER
  void origin() {
    assert(init && "boundary_stepper must be initialized before calling origin().");
    
    inside = true;
    i_boundary = 1;
    
    current_voxel = start_voxel;
    pathlength = boundaries[1].distance - boundaries[0].distance;
  }
  
  CUDA_CALLABLE_MEMBER
  void next() {
    current_voxel = boundaries[i_boundary].entering;
    i_boundary++;
    if (i_boundary > boundaries.size()-1)
      inside = false;
    else
      pathlength = boundaries[i_boundary].distance - boundaries[i_boundary-1].distance;
  }

};
  
#endif
