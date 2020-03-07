// boundaries.h -- structure to contain boundary crossing information for ray-grid intersections

#ifndef __boundaries_H
#define __boundaries_H

#include "atmo_vec.h"
#include <algorithm>
#include <vector>
#include <limits>

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

  unsigned int size() {
    return boundaries.size();
  }
  
  void sort() {
    std::sort(boundaries.begin(),boundaries.end());
  }

  void append(boundary b) {
    if (b.entering_indices.size() == n_dimensions) {
      boundaries.push_back(b);
    } else {
      std::cout << "boundary does not match dimensions of set.\n";
      throw(99);
    }
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
    unsigned int end = begin;
    do {
      end++;
    } while (boundaries[end].entering != -1 && end < boundaries.size()-1);
    
    if (begin == boundaries.size() || end == begin || boundaries[end].entering != -1) {
      std::cout << "boundary errors in ray_voxel_intersections.\n";
      throw(99);
    }

    boundaries = vector<boundary>(boundaries.begin()+begin, boundaries.begin()+end+1);
  }
};






struct boundary_intersection_stepper {
  bool init;
  atmo_vector vec;
  boundary_set boundaries;

  bool inside;
  unsigned int i_boundary;

  int start_voxel;
  int current_voxel;
  double pathlength;

  int n_emissions;
  vector<double> tau_species_initial;
  vector<double> tau_species_final;
  vector<double> tau_absorber_initial;
  vector<double> tau_absorber_final;

  boundary_intersection_stepper() { init=false; }

  boundary_intersection_stepper(atmo_vector vecc, boundary_set boundariess, int n_emissionss)
    : vec(vecc), boundaries(boundariess), n_emissions(n_emissionss)
  {
    start_voxel = boundaries[0].entering;
    
    tau_species_initial.resize(n_emissions,0.);
    tau_species_final.resize(n_emissions,0.);
    tau_absorber_initial.resize(n_emissions,0.);
    tau_absorber_final.resize(n_emissions,0.);

    init=true;
  }

  void origin() {
    if (init) {
      inside = true;
      i_boundary = 1;
      
      current_voxel = start_voxel;
      pathlength = boundaries[1].distance - boundaries[0].distance;
      
      tau_species_initial.resize(n_emissions,0.);
      tau_species_final.resize(n_emissions,0.);
      tau_absorber_initial.resize(n_emissions,0.);
      tau_absorber_final.resize(n_emissions,0.);
    } else {
      throw(99);
    }
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
