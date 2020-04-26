// boundaries.h -- structure to contain boundary crossing information for ray-grid intersections

#ifndef __boundaries_H
#define __boundaries_H

#include "cuda_compatibility.h"
#include "atmo_vec.h"
#include <algorithm>
#include <limits>
#include <iostream>
#include <stdexcept>


template <int NDIM>
struct boundary {
  int entering;//index of voxel the ray is entering
  static const int n_dimensions = NDIM;
  int entering_indices[NDIM];//indices in each dimension of the voxel
  double distance;//distance to the boundary crossing

  boundary() {
    entering = -2;
    for (int i=0;i<n_dimensions;i++)
      entering_indices[i] = -2;
    distance = -1;
  }

  void reset() {
    entering = -2;
    for (int i=0;i<n_dimensions;i++)
      entering_indices[i] = -2;
    distance = -1;
  }

  bool operator<(const boundary<NDIM> &rhs) const { return distance < rhs.distance; }
  bool operator<=(const boundary<NDIM> &rhs) const { return distance <= rhs.distance; }
};








template <int NDIM>
class boundary_set {
private:
  unsigned int begin;
  unsigned int internal_size;
  unsigned int internal_max_size;
  boundary<NDIM> *boundaries;
public:
  static const unsigned int n_dimensions = NDIM;

  boundary_set(const int max_size = 0)
  { 
    internal_size=0;begin=0;
    internal_max_size=max_size;
    
    boundaries = new boundary<NDIM>[internal_max_size];
  }
  ~boundary_set() {
    delete []boundaries;
  }
  boundary_set(const boundary_set &copy) {
    assert(n_dimensions == copy.n_dimensions && "only copy same dimension boundaries");

    begin = copy.begin;
    internal_size = copy.internal_size;
    internal_max_size = copy.internal_max_size;
    
    boundaries = new boundary<NDIM>[internal_max_size];

    for (unsigned int i=begin;i<begin+internal_size;i++)
      boundaries[i] = copy.boundaries[i];
  }
  boundary_set &operator=(const boundary_set &rhs) {
    assert(n_dimensions == rhs.n_dimensions && "only copy same dimension boundaries");
    if(this == &rhs) return *this;

    begin = rhs.begin;
    internal_size = rhs.internal_size;
    internal_max_size = rhs.internal_max_size;
    
    boundaries = new boundary<NDIM>[internal_max_size];

    for (unsigned int i=begin;i<begin+internal_size;i++)
      boundaries[i] = rhs.boundaries[i];

    return *this;
  }


  boundary<NDIM> operator[](int i) {
    return boundaries[begin+i];
  }

  boundary<NDIM> back() {
    return boundaries[begin+internal_size-1];
  }
  
  unsigned int size() {
    return internal_size;
  }

  void reset() {
    internal_size=0;begin=0;
  }
  
  void sort() {
    vector<boundary<NDIM>> vecbound(internal_size);
    for (unsigned int i=0;i<internal_size;i++)
      vecbound[i] = boundaries[begin+i];

    std::sort(vecbound.begin(),vecbound.end());

    for (unsigned int i=0;i<internal_size;i++)
      boundaries[begin+i] = vecbound[i];
  }

  void append(boundary<NDIM> b) {
    assert(internal_size < internal_max_size);
    boundaries[begin+internal_size] = b;
    internal_size++;
  }

  void add_intersections(const double start, const int dim, 
			 const int idx, const double &coordinate,
			 const double (&distances)[2], const int n_hits) {
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

  void propagate_indices() {
    //propogate voxel indices in each coordinate
    for (unsigned int i=1;i<size();i++) {
      for (unsigned int j=0;j<n_dimensions;j++)
	if (boundaries[i].entering_indices[j]==-2)
	  boundaries[i].entering_indices[j] = boundaries[i-1].entering_indices[j];
    }
  }


  template<typename C>
  void assign_voxel_indices(C *obj, void (C::*indices_to_voxel)(const int (&/*indices*/)[NDIM] , int &/*voxel*/) ) {  
    for (unsigned int i=1;i<size();i++) 
      (obj->indices_to_voxel)(boundaries[i].entering_indices, boundaries[i].entering);
  }

  
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

  
  bool check(const vector<int> &n_bounds, const int &n_voxels) {
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





template <int NDIM>
struct boundary_intersection_stepper {
  bool init;
  atmo_vector vec;
  boundary_set<NDIM> boundaries;

  bool inside;
  bool exits_bottom;
  bool exits_top;
  unsigned int i_boundary;

  int start_voxel;
  int current_voxel;
  double pathlength;

  boundary_intersection_stepper() : init(false) { }

  boundary_intersection_stepper(atmo_vector vecc, boundary_set<NDIM> boundariess)
    : vec(vecc), boundaries(boundariess)
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
  
  void origin() {
    assert(init && "boundary_stepper must be initialized before calling origin().");
    
    inside = true;
    i_boundary = 1;
    
    current_voxel = start_voxel;
    pathlength = boundaries[1].distance - boundaries[0].distance;
  }
  
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
