//intersections.h -- routines for computing intersections between lines and geometric primitives

#include "intersections.hpp"
#include <cmath>

CUDA_CALLABLE_MEMBER
bool geom_primitive::samesign(Real a, Real b) const {
  if ((a>0&&b>0) || (a<0&&b<0) || (a==0&&b==0))
    return true;
  else
    return false;
}
CUDA_CALLABLE_MEMBER
bool geom_primitive::is_zero(Real a, Real tol) const {
  if (a > tol || a < -tol)
    return false;
  else
    return true;
}



//plane

CUDA_CALLABLE_MEMBER
void plane::intersections(const atmo_vector & vec,
			  Real (&distances)[2],
			  int &n_hits) const {
  n_hits = 0;

  if (vec.line_z != 0) {
    Real d = (z-vec.pt.z)/vec.line_z;
    if (d > 0) {
      distances[n_hits]=d;
      n_hits++;
    }
  }
    
  for (int i=0;i<n_hits;i++) 
    assert(is_zero(vec.extend(distances[i]).z/z-1.0)
	   && "vector must intersect plane at specified distance.");
}




//sphere

void sphere::set_radius(const Real &rr) {
  r=rr/scale;
  r2=r*r;
}

CUDA_CALLABLE_MEMBER
void sphere::intersections(const atmo_vector & vec,
			   Real (&distances)[2],
			   int &n_hits) const {
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
    assert(is_zero(vec.extend(distances[i]).r/r/scale-1.0,ABS)
	   && "vector must intersect sphere at specified distance.");
}





//cone

void cone::set_angle(const Real &a) {
  angle=a;
  cosangle=std::cos(angle);
  cosangle2=cosangle*cosangle;
  assert(!(is_zero(cosangle))
	 && "problems occur if there is a cone with an opening angle of pi/2 degrees.");
}

CUDA_CALLABLE_MEMBER
void cone::intersections(const atmo_vector & vec,
			 Real (&distances)[2],
			 int &n_hits) const {
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
    //this test needs work, still needs very large error term to pass on floats
    assert(is_zero(vec.extend(distances[i]).t-angle,CONEABS)
	   && "vector must intersect cone at specified distance.");
  }
}


