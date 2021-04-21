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
    
  for (int i=0;i<n_hits;i++) {
    atmo_point ipt = vec.extend(distances[i]);
    assert(is_zero(ipt.z/z-1.0)
	   && "vector must intersect plane at specified distance.");
  }
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

  for (int i=0;i<n_hits;i++) {
    atmo_point ipt = vec.extend(distances[i]);
    assert(is_zero(ipt.r/r/scale-1.0,EPS)
	   && "vector must intersect sphere at specified distance.");
  }
}





//cone

void cone::set_angle(const Real &a) {
  angle=a;
  cosangle=std::cos(angle);
  cosangle2=cosangle*cosangle;
  assert(!(is_zero(cosangle))
	 && "problems occur if there is a cone with an opening angle of pi/2 degrees.");
}

void cone::set_rmin(const Real &rminn) {
  //set the minimum radius at which to consider cone
  //intersections. This solves a floating point rounding issue when
  //intersections are near the origin and far from the original
  //location
  rmin = rminn;
}
  

CUDA_CALLABLE_MEMBER
void cone::intersections(const atmo_vector & vec,
			 Real (&distances)[2],
			 int &n_hits) const {
  n_hits = 0;

  const Real rscale = vec.pt.r;
  const Real z_norm = vec.pt.z/rscale;
  //  const Real r_norm = vec.pt.r/rscale;

  const Real A = vec.line_z * vec.line_z - cosangle2;
  const Real B = z_norm * vec.line_z - vec.ray.cost * cosangle2;
  const Real C = z_norm * z_norm - cosangle2;
    
  if (!is_zero(A)) {
    const Real discr = B*B-A*C;
    // const Real discr = cosangle2*(vec.line_z*vec.line_z
    // 				  +z_norm*z_norm
    // 				  -2*z_norm*vec.line_z*vec.ray.cost
    // 				  -cosangle2*vec.ray.sint*vec.ray.sint);
    
    if (discr > 0) {

      const Real q = (B > 0) ? 
	-B - std::sqrt(discr) : 
	-B + std::sqrt(discr); 

      const Real d0 = q/A;
      if (d0 > 0 && samesign(z_norm + d0*vec.line_z, cosangle) ) {
	distances[n_hits]=d0*rscale;
	n_hits++;
      }
      const Real d1 = C/q;
      if (d1 > 0 && samesign(z_norm + d1*vec.line_z, cosangle)) {
	distances[n_hits]=d1*rscale;
	n_hits++;
      }
    }
  } else {
    const Real d = -C/(2*B);
    if (d>0 && samesign(z_norm + d*vec.line_z, cosangle)) {
      distances[n_hits]=d*rscale;
      n_hits++;
    }
  }

  for (int i=0;i<n_hits;i++) {
    atmo_point ipt = vec.extend(distances[i]);
    if (ipt.r > rmin)
      //This conditional solves a floating point rounding issue when
      //intersections are near the origin and far from the original
      //location
      assert(is_zero(ipt.t/angle-1,CONEEPS)
	     && "vector must intersect cone at specified distance.");
    //we don't need to discard intersections with r<rmin because the
    //integration code does this automatically
  }
}


