//intersections.h -- routines for computing intersections between lines and geometric primitives

#ifndef __INTERSECTIONS_H
#define __INTERSECTIONS_H

#include "cuda_compatibility.h"
#include "atmo_vec.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

using std::vector;

class geom_primitive {
 protected:
  CUDA_CALLABLE_MEMBER
  bool samesign(double a, double b) const {
    if ((a>0&&b>0) || (a<0&&b<0) || (a==0&&b==0))
      return true;
    else
      return false;
  }
  CUDA_CALLABLE_MEMBER
  bool is_zero(double a) const {
    double tol = 1e-10;
    if (a > tol || a < -tol)
      return false;
    else
      return true;
  }
};


class plane : geom_primitive {
 public:
  double z;

  CUDA_CALLABLE_MEMBER
  void intersections(const atmo_vector & vec, double (&distances)[2], int &n_hits) const {
    n_hits = 0;

    if (vec.line_z != 0) {
      double d = (z-vec.pt.z)/vec.line_z;
      if (d > 0) {
	distances[n_hits]=d;
	n_hits++;
      }
    }
    
    for (int i=0;i<n_hits;i++) 
      assert(is_zero(vec.extend(distances[i]).z/z-1.0) && "vector must intersect plane at specified distance.");
  }
  
};



class sphere : geom_primitive {
protected:
  double r;
  double r2;

public:
  void set_radius(const double &rr) {
    r=rr;
    r2=r*r;
  }

  CUDA_CALLABLE_MEMBER
  double get_r() {
    return r;
  }
  
  CUDA_CALLABLE_MEMBER
  void intersections(const atmo_vector & vec, double (&distances)[2], int &n_hits) const {
    n_hits = 0;

    // we solve the quadratic t^2 - 2*(vec.r * vec.cost)*t + (vec.r^2 - r^2) == 0
    double B = vec.pt.r * vec.ray.cost;
    double C = vec.pt.r*vec.pt.r - r2;

    double discr = B*B-C;

    //we can ignore cases where the ray misses or is tangent
    if (discr > 0) {
      discr = std::sqrt(discr);
      
      double d0 = -B - discr;
      if (d0>0) {
	distances[n_hits]=d0;
	n_hits++;
      }
      double d1 = -B + discr;
      if (d1>0) {
	distances[n_hits]=d1;
	n_hits++;
      }
    }

    for (int i=0;i<n_hits;i++) 
      assert(is_zero(vec.extend(distances[i]).r/r-1.0) && "vector must intersect sphere at specified distance.");
  }

};


class cone : geom_primitive {
protected:  
  double angle;
  double cosangle;
  double cosangle2;

public:
  
  void set_angle(const double &a) {
    angle=a;
    cosangle=std::cos(angle);
    cosangle2=cosangle*cosangle;
    assert(!(is_zero(cosangle)) && "problems occur if there is a cone with an opening angle of pi/2 degrees.");
  }

  CUDA_CALLABLE_MEMBER
  void intersections(const atmo_vector & vec, double (&distances)[2], int &n_hits) const {
    n_hits = 0;

    double A = vec.line_z * vec.line_z - cosangle2;
    double B = vec.pt.z * vec.line_z - vec.pt.r * vec.ray.cost * cosangle2;
    double C = vec.pt.z * vec.pt.z - vec.pt.r * vec.pt.r * cosangle2;
    
    if (!is_zero(A)) {
      double discr = B*B-A*C;

      if (discr > 0) {
	discr = std::sqrt(discr);
	
	double d0 = (-B - discr)/A;
	if (d0 > 0 && samesign(vec.pt.z + d0*vec.line_z, cosangle) ) {
	  distances[n_hits]=d0;
	  n_hits++;
	}
	double d1 = (-B + discr)/A;
	if (d1 > 0 && samesign(vec.pt.z + d1*vec.line_z, cosangle)) {
	  distances[n_hits]=d1;
	  n_hits++;
	}
      }
    } else {
      double d = -C/(2*B);
      if (d>0 && samesign(vec.pt.z + d*vec.line_z, cosangle)) {
	distances[n_hits]=d;
	n_hits++;
      }
    }

    for (int i=0;i<n_hits;i++) 
      assert(is_zero(vec.extend(distances[i]).t/angle-1.0) && "vector must intersect cone at specified distance.");
  }
};


// struct intersection_writer {
//   string fname;
//   std::ofstream file;

//   intersection_writer() { }
  
//   intersection_writer(string fnamee) : fname(fnamee) { }

//   ~intersection_writer() {
//     if (file.is_open()) {
//       file.close();
//     }
//   }
  
//   void save_coordinates(vector<double> radial_boundaries,
// 			vector<double> sza_boundaries) {
//     file.open(fname.c_str());
//     if (file.is_open()) {
//       VectorXd r_boundaries_write_out = Eigen::Map<VectorXd>(radial_boundaries.data(),
// 							     radial_boundaries.size());
    
//       file << "radial boundaries [cm]:\n" << r_boundaries_write_out.transpose() << "\n";
    
//       VectorXd sza_boundaries_write_out = Eigen::Map<VectorXd>(sza_boundaries.data(),
// 							       sza_boundaries.size());
    
//       file << "sza boundaries [rad]:\n" << sza_boundaries_write_out.transpose() << "\n\n";
//     }
//   }

//   template <int NDIM>
//   void append_intersections(const atmo_vector &vec,
// 			    const boundary_set<NDIM> &boundaries) {
//     if (file.is_open()) {
      
//       file << "ray origin (xyz) [cm]:\n"
// 	   << vec.pt.x << "  "
// 	   << vec.pt.y << "  "
// 	   << vec.pt.z << "\n";
      
//       file << "ray direction (xyz):\n"
// 	   << vec.line_x << "  "
// 	   << vec.line_y << "  "
// 	   << vec.line_z << "\n";
      
      
//       file << "Intersections [t, r_bound, sza_bound]\n";
//       for (unsigned int i=0;i<boundaries.size();i++) {
// 	file << boundaries[i].distance << "  "
// 	     << boundaries[i].entering_indices[0] << "  "
// 	     << boundaries[i].entering_indices[1] << "\n";
//       }
//       file << "\n";
//     }
//   }
// };
  


#endif
