//interp_1d.h -- base class for interpolation routines

#ifndef __Interp_H_
#define __Interp_H_

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "Real.hpp"

using std::pow;
using std::vector;

struct Base_interp
// Abstract base class used by all interpolation routines in this
// chapter. Only interp is called by the user.
{
  int n, mmm, dj;
  mutable int jsav, cor;
  Real *xx, *yy;

  Base_interp() { }
  
  Base_interp(vector<Real> &x, vector<Real> &y, int m)
  // constructor. set up for interpolating on a table of x's and y's
  // of length m. Normally called by derived class, not user.
    : n(x.size()), mmm(m), jsav(0), cor(0) {
    if (n==0) {
      xx = NULL;
      yy = NULL;
    } else {
      xx = new Real[n];
      yy = new Real[n];
    }
    for (int i=0;i<n;i++) {
      xx[i]=x[i];
      yy[i]=y[i];
    }
    
    dj = std::min(1, (int) pow((Real) n, 0.25) );
  }

  Base_interp(const Base_interp &B) {
    n=B.n;
    mmm=B.mmm;
    jsav=B.jsav;
    cor=B.cor;
    dj=B.dj;
    if (n==0) {
      xx = NULL;
      yy = NULL;
    } else {
      xx = new Real[n];
      yy = new Real[n];
    }
    for (int i=0; i<n; i++) {
      xx[i]=B.xx[i];
      yy[i]=B.yy[i];
    }
  }

  Base_interp& operator= (const Base_interp &B) {
    n=B.n;
    mmm=B.mmm;
    jsav=B.jsav;
    cor=B.cor;
    dj=B.dj;
    if (n==0) {
      xx = NULL;
      yy = NULL;
    } else {
      xx = new Real[n];
      yy = new Real[n];
    }
    for (int i=0; i<n; i++) {
      xx[i]=B.xx[i];
      yy[i]=B.yy[i];
    }
    return *this;
  }

  ~Base_interp() {
    if (n>0) {
      delete [] xx;
      delete [] yy;
    }
  }


  Real operator()(const Real x) const 
  // Given a value x, return an interpolated value, using data
  // pointed to by xx and yy.
  {
    int jlo = cor ? hunt(x) : locate (x);
    return rawinterp(jlo, x);
  }

  int index(const Real x) const
  // Given a value x, return an interpolated value, using data
  // pointed to by xx and yy.
  {
    int jlo = cor ? hunt(x) : locate(x);
    return jlo;
  }

  int locate(const Real x) const;
  int hunt(const Real x) const;

  virtual Real rawinterp(const int jlo, const Real x) const = 0;
  // derived classes provide this as the actual interpolation method used.

};

struct Linear_interp : Base_interp
// piecewise linear interpolation object. Construct with x and y
// vectors, then call interp for interpolated values
{
  Linear_interp() : Base_interp() {}

  Linear_interp(vector<Real> &xv, vector<Real> &yv) : Base_interp(xv,yv,2) {}
  
  Real rawinterp(const int j, const Real x) const {
    if(xx[j] == xx[j+1]) return yy[j]; // table defective, but recover
    else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
  }
};

struct Bilinear_interp
// object for bilinear interpolation on a matrix. Construct with a
// vector of x1 values, a vector of x2 values, and a matrix of
// tabulated function values yij. Then call interp for interpolated
// values.
{
  int m, n;
  MatrixX y;
  Linear_interp x1terp, x2terp;

  Bilinear_interp() { }
  
  Bilinear_interp(vector<Real> &x1v, vector<Real> &x2v, MatrixX &ym)
    : m(x1v.size()), n(x2v.size()), x1terp(x1v,x1v), x2terp(x2v,x2v) { y = ym; }
    // we need dummy 1 dim interp objects for their locate and hunt functions

  Real interp(const Real x1p, const Real x2p) const {
    int i, j;
    Real yy, t, u;
    //find the grid square:
    i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
    j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);

    //interpolate:
    t = (x1p-x1terp.xx[i])/(x1terp.xx[i+1]-x1terp.xx[i]);
    u = (x2p-x2terp.xx[j])/(x2terp.xx[j+1]-x2terp.xx[j]);
    yy =  (1.-t)*(1.-u) * y(i  ,j  )
	 +    t *(1.-u) * y(i+1,j  )
	 +(1.-t)*    u  * y(i  ,j+1)
	 +    t *    u  * y(i+1,j+1);

    return yy;
  }

};

#endif
