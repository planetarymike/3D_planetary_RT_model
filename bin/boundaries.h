// boundaries.h -- structure to contain boundary crossing information for ray-grid intersections

#ifndef __boundaries_H
#define __boundaries_H

#include "atmo_vec.h"
#include <algorithm>
#include <vector>
#include <limits>
#include <iostream>
#include <stdexcept>

using std::vector;


struct boundary {
  int entering;//index of voxel the ray is entering
  vector<int> entering_indices;//indices in each dimension of the voxel
  double distance;//distance to the boundary crossing

  boundary(int n_dimensions) {
    entering_indices.resize(n_dimensions,-2);
  }

  bool operator<(const boundary &rhs) const { return distance < rhs.distance; }
  bool operator<=(const boundary &rhs) const { return distance <= rhs.distance; }
};










struct boundary_set {
  unsigned int n_dimensions;
  vector<boundary> boundaries;

  boundary_set() : n_dimensions(0) { }

  boundary_set(int n_dim) : n_dimensions(n_dim) { }
  
  boundary operator[](int i) {
    return boundaries[i];
  }

  boundary back() {
    return boundaries.back();
  }
  
  unsigned int size() {
    return boundaries.size();
  }
  
  void sort() {
    std::sort(boundaries.begin(),boundaries.end());
  }

  void append(boundary b) {
    assert(b.entering_indices.size() == n_dimensions && "boundary must match dimensions of set.");
    boundaries.push_back(b);
  }

  template<class T>
  void add_intersections(const double start, const int dim, 
			 const int idx, const double &coordinate, T distances) {
    boundary new_boundary(n_dimensions);

    std::sort(distances.begin(),distances.end());
    
    if (start > coordinate) 
      new_boundary.entering_indices[dim]=idx-1;
    else 
      new_boundary.entering_indices[dim]=idx;
    
    new_boundary.distance=distances[0];

    boundaries.push_back(new_boundary);

    if (distances.size()>1) {
      if (start > coordinate) 
	new_boundary.entering_indices[dim]=idx;
      else 
	new_boundary.entering_indices[dim]=idx-1;
      new_boundary.distance=distances[1];
      boundaries.push_back(new_boundary);
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
    unsigned int begin = 0;
    while (boundaries[begin].entering == -1 && begin < boundaries.size()-1)
      begin++;

    if (begin==boundaries.size()-1) {
      boundaries.clear();
    } else {
      unsigned int end = begin;
      do {
	end++;
      } while (boundaries[end].entering != -1 && end < boundaries.size()-1);

      assert(end < boundaries.size() && "end must be inside list of boundary intersections");
      assert(boundaries[end].entering == -1 && "trim error in boundary_set: ray does not exit grid");
    
      boundaries = vector<boundary>(boundaries.begin()+begin, boundaries.begin()+end+1);
    }
  }

  
  void check(const vector<int> &n_bounds, const int &n_voxels) {
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
  }
};






struct boundary_intersection_stepper {
  bool init;
  atmo_vector vec;
  boundary_set boundaries;

  bool inside;
  bool exits_bottom;
  bool exits_top;
  unsigned int i_boundary;

  int start_voxel;
  int current_voxel;
  double pathlength;

  boundary_intersection_stepper() : init(false) { }

  boundary_intersection_stepper(atmo_vector vecc, boundary_set boundariess)
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
