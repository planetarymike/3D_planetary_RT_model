// boundaries.h -- structure to contain boundary crossing information for ray-grid intersections

#ifndef __boundaries_H
#define __boundaries_H

#include "cuda_compatibility.h"
#include "atmo_vec.h"
#include <algorithm>
#include <vector>
#include <limits>
#include <iostream>
#include <stdexcept>

using std::vector;

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

  void set_entering_indices(const vector<int> &rhs) {
    assert(rhs.size() == n_dimensions && "vector indices must have same dimensionality as grid.");
    for (int i=0;i<n_dimensions;i++)
      entering_indices[i] = rhs[i];
  }
  
  bool operator<(const boundary<NDIM> &rhs) const { return distance < rhs.distance; }
  bool operator<=(const boundary<NDIM> &rhs) const { return distance <= rhs.distance; }
};








template <int NDIM>
class boundary_set {
  //some savings to be found here in avoiding consistent reallocation of vector<boundary>
  //use a fixed-size array and navigate with internal pointers begin, end, and internal_size
private:
  unsigned int begin;
  unsigned int end;
  int internal_size;
  int internal_max_size;
public:
  static const unsigned int n_dimensions = NDIM;
  vector<boundary<NDIM>> boundaries;

  boundary_set(const int max_size = 200)
  { 
    internal_size=0;begin=0;end=0;
    
    boundaries.reserve(max_size);
  }

  boundary<NDIM> operator[](int i) {
    return boundaries[begin+i];
  }

  boundary<NDIM> back() {
    return boundaries.back();
  }
  
  unsigned int size() {
    return internal_size;
  }

  void reset() {
    internal_size=0;begin=0;end=0;
    boundaries.clear();
  }
  
  void sort() {
    std::sort(boundaries.begin(),boundaries.end());
  }

  void append(boundary<NDIM> b) {
    boundaries.push_back(b);
    internal_size++;
  }

  void add_intersections(const double start, const int dim, 
			 const int idx, const double &coordinate,
			 const double (&distances)[2], const int n_hits) {
    static thread_local boundary<NDIM> new_boundary;
    new_boundary.reset();

    assert( 0 <= n_hits && n_hits <=2 && "only 0, 1, or 2 intersections are possible with a convex shell");
    
    if (n_hits == 0) {
      return;
    } else if (n_hits == 1) {
      bool above = (start > coordinate);
      

      new_boundary.entering_indices[dim] = above ? idx-1 : idx;
      new_boundary.distance=distances[0];

      boundaries.push_back(new_boundary);
      internal_size++;
    } else { // n_hits == 2

      bool above = (start > coordinate);
      bool in_order = (distances[1] > distances[0]);

      new_boundary.entering_indices[dim] = above ? idx-1 : idx;
      new_boundary.distance = in_order ? distances[0] : distances[1];
      boundaries.push_back(new_boundary);
      internal_size++;
      
      new_boundary.entering_indices[dim] = above ? idx : idx-1;
      //^^ indices are reversed to come back out of the region we just entered
      new_boundary.distance = in_order ? distances[1] : distances[0];
      boundaries.push_back(new_boundary);
      internal_size++;
    }
  }

  void propagate_indices() {
    //propogate voxel indices in each coordinate
    for (unsigned int i=1;i<boundaries.size();i++) {
      for (unsigned int j=0;j<n_dimensions;j++)
	if (boundaries[i].entering_indices[j]==-2)
	  boundaries[i].entering_indices[j] = boundaries[i-1].entering_indices[j];
    }
  }


  template<typename C>
  void assign_voxel_indices(C *obj, int (C::*indices_to_voxel)(const vector<int>& ) ) {  
    for (unsigned int i=1;i<boundaries.size();i++) 
      boundaries[i].entering = (obj->indices_to_voxel)(boundaries[i].entering_indices);
  }

  
  void trim() {
    //trim the list of intersections to those inside the grid
    while (boundaries[begin].entering == -1 && begin < boundaries.size()-1)
      begin++;

    if (begin==boundaries.size()-1) {
      begin=end=internal_size=0;
    } else {
      unsigned int end = begin;
      do {
	end++;
      } while (boundaries[end].entering != -1 && end < boundaries.size()-1);

      assert(end < boundaries.size() && "end must be inside list of boundary intersections");
      assert(boundaries[end].entering == -1 && "trim error in boundary_set: ray does not exit grid");
    
      internal_size=end-begin+1;
    }
  }

  
  bool check(const vector<int> &n_bounds, const int &n_voxels) {
    if (boundaries.size() > 0) {
      assert(boundaries.size() > 1 && "there must be more than one boundary crossing for each ray");
      
      for (unsigned int i=1;i<boundaries.size();i++) {
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

	assert((boundaries.back().entering_indices[0] == -1 ||
		boundaries.back().entering_indices[0] == n_bounds[0])
	       && "ray must exit grid radially via the top or bottom");


	assert(boundaries[i].entering < n_voxels
	       && boundaries[i].entering >= -1
	       && "voxel index must be in bounds.");

      }
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
