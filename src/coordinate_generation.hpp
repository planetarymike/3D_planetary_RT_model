//coordinate_generation.h -- routines to get point and ray coordinates

#ifndef __coordinate_generation_h
#define __coordinate_generation_h

#include "Real.hpp"
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "constants.hpp"
#include "atmo_vec.hpp"
#include "atm/atmosphere_base.hpp"
#include "interp.hpp"
#include "gauss_legendre_quadrature.hpp" // numerical recipies gauss quadrature points.
// would switch to boost, but gauss weights are templated, not dynamic



using std::vector;
using std::exp;
using std::log;

template <class V>
void gauss_quadrature_points(V &pts,
			     V &wts,
			     Real start,
			     Real end,
			     int npts) {

  pts.resize(npts);
  wts.resize(npts);
  
  gauleg(start,end,pts,wts);
  
}

template <class V>
void uniform_quadrature_points(V &pts,
			       V &wts,
			       Real start,
			       Real end,
			       int npts,
			       bool cyclic = false,
			       Real offset = 0.0) {

  pts.resize(npts);
  wts.resize(npts);
  
  int ndivisions;
  if (cyclic) 
    ndivisions = npts+1;
  else 
    ndivisions = npts;

  Real step = ( end - start )/Real(ndivisions-1);

  for (int i = 0; i<npts; i++) {
    pts[i] = start + i * step;
    if (cyclic)
      pts[i]+=step*offset;
    wts[i] = step;
  }
  
}

template <class V>
void get_radial_log_linear_points(V &rpts,
				  int nrpts,
				  Real rminatm,
				  Real rexo,
				  Real rmaxatm) {
  // gets radial points that are split, with half linearly spaced below
  //   the exobase and half logarithmically spaced above. 
  
  int nbelowrexo=nrpts/2;
  
  Real logmax = log(rmaxatm-rMars);
  Real logmin = log(rexo-rMars);
  Real logspace = (logmax-logmin)/((Real) nrpts-nbelowrexo);

  Real linspace = (rexo-rminatm)/((Real) nbelowrexo-1);

  rpts.clear();
  for (int i = 0; i<nrpts; i++) {
    if (i<nbelowrexo) 
      rpts.push_back(rminatm+i*linspace);
    else
      rpts.push_back(
		     exp(
			 logmin + (i-nbelowrexo+1)*logspace )
		     + rMars);
  }

}


#endif
