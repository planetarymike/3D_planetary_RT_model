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
  static constexpr Real scale = 1e9;

  CUDA_CALLABLE_MEMBER
  bool samesign(Real a, Real b) const {
    if ((a>0&&b>0) || (a<0&&b<0) || (a==0&&b==0))
      return true;
    else
      return false;
  }
  CUDA_CALLABLE_MEMBER
  bool is_zero(Real a, Real tol=STRICTABS) const {
    if (a > tol || a < -tol)
      return false;
    else
      return true;
  }
};


class plane : geom_primitive {
 public:
  Real z;

  CUDA_CALLABLE_MEMBER
  void intersections(const atmo_vector & vec, Real (&distances)[2], int &n_hits) const {
    n_hits = 0;

    if (vec.line_z != 0) {
      Real d = (z-vec.pt.z)/vec.line_z;
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
  Real r;
  Real r2;

public:
  void set_radius(const Real &rr) {
    r=rr/scale;
    r2=r*r;
  }

  CUDA_CALLABLE_MEMBER
  Real get_r() {
    return r;
  }
  
  CUDA_CALLABLE_MEMBER
  void intersections(const atmo_vector & vec, Real (&distances)[2], int &n_hits) const {
    n_hits = 0;

    // we solve the quadratic t^2 - 2*(vec.r * vec.cost)*t + (vec.r^2 - r^2) == 0
    const Real r_norm = vec.pt.r/scale;
    const Real B = r_norm * vec.ray.cost;
    const Real C = r_norm * r_norm - r2;

    const Real discr = B*B-C;

    //we can ignore cases where the ray misses or is tangent
    if (discr > 0) {

      const Real d0 = (B > 0) ? 
	-B - std::sqrt(discr) : 
	-B + std::sqrt(discr); 
      if (d0>0) {
	distances[n_hits]=d0*scale;
	n_hits++;
      }
      const Real d1 = C / d0; 
      if (d1>0) {
	distances[n_hits]=d1*scale;
	n_hits++;
      }
    }

    for (int i=0;i<n_hits;i++)
      assert(is_zero(vec.extend(distances[i]).r/r/scale-1.0,ABS) && "vector must intersect sphere at specified distance.");
  }

};


class cone : geom_primitive {
protected:  
  Real angle;
  Real cosangle;
  Real cosangle2;

public:
  
  void set_angle(const Real &a) {
    angle=a;
    cosangle=std::cos(angle);
    cosangle2=cosangle*cosangle;
    assert(!(is_zero(cosangle)) && "problems occur if there is a cone with an opening angle of pi/2 degrees.");
  }

  CUDA_CALLABLE_MEMBER
  void intersections(const atmo_vector & vec, Real (&distances)[2], int &n_hits) const {
    n_hits = 0;

    const Real z_norm = vec.pt.z/scale;
    const Real r_norm = vec.pt.r/scale;

    const Real A = vec.line_z * vec.line_z - cosangle2;
    const Real B = z_norm * vec.line_z - r_norm * vec.ray.cost * cosangle2;
    const Real C = z_norm * z_norm - r_norm * r_norm * cosangle2;
    
    if (!is_zero(A)) {
      const Real discr = B*B-A*C;

      if (discr > 0) {

	const Real q = (B > 0) ? 
	  -B - std::sqrt(discr) : 
	  -B + std::sqrt(discr); 

	const Real d0 = q/A;
	if (d0 > 0 && samesign(z_norm + d0*vec.line_z, cosangle) ) {
	  distances[n_hits]=d0*scale;
	  n_hits++;
	}
	const Real d1 = C/q;
	if (d1 > 0 && samesign(z_norm + d1*vec.line_z, cosangle)) {
	  distances[n_hits]=d1*scale;
	  n_hits++;
	}
      }
    } else {
      const Real d = -C/(2*B);
      if (d>0 && samesign(z_norm + d*vec.line_z, cosangle)) {
	distances[n_hits]=d*scale;
	n_hits++;
      }
    }

    for (int i=0;i<n_hits;i++) {
      //this test needs work, still needs very large error term to pass
      //std::cout << vec.extend(distances[i]).t-angle << std::endl;
      assert(is_zero(vec.extend(distances[i]).t-angle,CONEABS) && "vector must intersect cone at specified distance.");
    }
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
  
//   void save_coordinates(vector<Real> radial_boundaries,
// 			vector<Real> sza_boundaries) {
//     file.open(fname.c_str());
//     if (file.is_open()) {
//       VectorX r_boundaries_write_out = Eigen::Map<VectorX>(radial_boundaries.data(),
// 							     radial_boundaries.size());
    
//       file << "radial boundaries [cm]:\n" << r_boundaries_write_out.transpose() << "\n";
    
//       VectorX sza_boundaries_write_out = Eigen::Map<VectorX>(sza_boundaries.data(),
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
