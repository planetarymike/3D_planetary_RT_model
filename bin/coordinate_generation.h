//coordinate_generation.h -- routines to get point and ray coordinates

#ifndef __coordinate_generation_h
#define __coordinate_generation_h

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "constants.h"
#include "atmo_vec.h"
#include "atmosphere.h"
#include "interp.h"
#include "gauss_legendre_quadrature.h" // numerical recipies gauss quadrature points.
// would switch to boost, but gauss weights are templated, not dynamic



using std::vector;
using std::exp;
using std::log;

template <class V>
void gauss_quadrature_points(V &pts,
			     V &wts,
			     double start,
			     double end,
			     int npts) {

  pts.resize(npts);
  wts.resize(npts);
  
  gauleg(start,end,pts,wts);
  
}

template <class V>
void uniform_quadrature_points(V &pts,
			       V &wts,
			       double start,
			       double end,
			       int npts,
			       bool cyclic = false,
			       double offset = 0.0) {

  pts.resize(npts);
  wts.resize(npts);
  
  int ndivisions;
  if (cyclic) 
    ndivisions = npts+1;
  else 
    ndivisions = npts;

  double step = ( end - start )/double(ndivisions-1);

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
				  double rminatm,
				  double rexo,
				  double rmaxatm) {
  // gets radial points that are split, with half linearly spaced below
  //   the exobase and half logarithmically spaced above. 
  
  int nbelowrexo=nrpts/2;
  
  double logmax = log(rmaxatm-rMars);
  double logmin = log(rexo-rMars);
  double logspace = (logmax-logmin)/((double) nrpts-nbelowrexo);

  double linspace = (rexo-rminatm)/((double) nbelowrexo-1);

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
