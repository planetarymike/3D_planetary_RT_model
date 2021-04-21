//intersections.h -- routines for computing intersections between lines and geometric primitives

#ifndef __INTERSECTIONS_H
#define __INTERSECTIONS_H

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "atmo_vec.hpp"

class geom_primitive {
 protected:
  static constexpr Real scale = 1e9;

  CUDA_CALLABLE_MEMBER
  bool samesign(Real a, Real b) const;
  CUDA_CALLABLE_MEMBER
  bool is_zero(Real a, Real tol=STRICTEPS) const;
};


class plane : geom_primitive {
 public:
  Real z;

  CUDA_CALLABLE_MEMBER
  void intersections(const atmo_vector & vec, 
		     Real (&distances)[2], 
		     int &n_hits) const;  
};



class sphere : geom_primitive {
protected:
  Real r;
  Real r2;

public:
  void set_radius(const Real &rr);

  CUDA_CALLABLE_MEMBER
  void intersections(const atmo_vector & vec, 
		     Real (&distances)[2], 
		     int &n_hits) const;
};


class cone : geom_primitive {
protected:  
  Real angle;
  Real cosangle;
  Real cosangle2;

  Real rmin;
public:
  void set_angle(const Real &a);
  void set_rmin(const Real &rminn);
  
  CUDA_CALLABLE_MEMBER
  void intersections(const atmo_vector & vec,
		     Real (&distances)[2],
		     int &n_hits) const;
};

  


#endif
