//interp_1d.h -- base class for interpolation routines

#ifndef __Interp_H_
#define __Interp_H_

#include <cmath>
#include <vector>
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

  Base_interp() : n(0) { }
  Base_interp(const vector<Real> &x, const vector<Real> &y, const int m);
  Base_interp(const Base_interp &B);
  Base_interp& operator= (const Base_interp &B);
  ~Base_interp();

  Real operator()(const Real x) const;

  int index(const Real x) const;

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
  Linear_interp(const vector<Real> &xv, const vector<Real> &yv);
  
  Real rawinterp(const int j, const Real x) const;
};

struct Bilinear_interp
// object for bilinear interpolation on a matrix. Construct with a
// vector of x1 values, a vector of x2 values, and a matrix of
// tabulated function values yij. Then call interp for interpolated
// values.
{
  int m, n;
  MatrixX y;
  // we need dummy 1 dim interp objects for their locate and hunt functions
  Linear_interp x1terp, x2terp;

  Bilinear_interp() { }
  Bilinear_interp(const vector<Real> &x1v, const vector<Real> &x2v, const MatrixX &ym);

  Real interp(const Real x1p, const Real x2p) const;
};

#endif
