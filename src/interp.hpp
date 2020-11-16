//interp_1d.h -- base class for interpolation routines

#ifndef __Interp_H_
#define __Interp_H_

#include <cmath>
#include <vector>
#include <Eigen/Dense>

using std::pow;
using std::vector;

template <typename Real>
struct Base_interp
// Abstract base class used by all interpolation routines in this
// chapter. Only interp is called by the user.
{
  int n, mmm, dj;
  mutable int jsav, cor;
  Real *xx, *yy;

  Base_interp() : n(0) { }
  Base_interp(const vector<Real> &x, const vector<Real> &y, const int m)
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

  Base_interp(const Base_interp<Real> &B) {
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
  
  Base_interp<Real>& operator= (const Base_interp<Real> &B) {
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
  
  int locate(const Real x) const
  // given a value x, return a value j such that x is (insofar as
  // possible) centered in the subrange xx[j..j+mmm-1], where xx is the
  // stored pointer. The values in xx must be monotonic, either
  // increasing or decreasing. The returned value is not less than 0,
  // nor greater than n-1).
  {
    int ju, jm, jl;
    //  std::cout << "n = " << n << ", mmm = " << mmm << "\n";
    if (n < 2 || mmm < 2 || mmm > n) {
      //    std::cout << "n = " << n << ", mmm = " << mmm << std::endl;
      assert(false && "locate size error");
    }
    
    bool ascnd = (xx[n-1] >= xx[0]); // true if ascending order of
    // table, false otherwise
    jl = 0;//lower and upper size limits
    ju = n-1;
    while (ju-jl > 1) {  // until the appropriate range is found
      jm = (ju+jl) >> 1; // find the midpoint
      if ((x >= xx[jm]) == ascnd)
	jl = jm; // and replace either the lower limit
      else
	ju = jm; // or the upper limit, as appropriate
    }
    cor = abs(jl-jsav) > dj ? 0 : 1; // decide whether to use hunt() or
    // locate() next time
    jsav = jl;
    return std::max(0,std::min(n-mmm,jl-((mmm-2)>>1)));
  }

  int hunt(const Real x) const
  // given a value x, return a value j such that x is (insofar as
  // possible) centered in the subrange xx[j..j+mmm-1], where xx is the
  // stored pointer. The values in xx must be monotonic, either
  // increasing or decreasing. The returned value is not less than 0,
  // nor greater than n-1).
  {
    int jl = jsav, jm, ju, inc = 1;
    //  std::cout << "n = " << n << ", mmm = " << mmm << std::endl;
    if (n < 2 || mmm < 2 || mmm > n) {
      //    std::cout << "n = " << n << ", mmm = " << mmm << std::endl;
      assert(false && "hunt size error");
    }
    bool ascnd = (xx[n-1] >= xx[0]); // does table ascend?
    if (jl < 0 || jl > n-1) {
      // input guess not useful. go directly to bisection
      jl = 0;
      ju = n-1;
    } else {
      if ((x >= xx[jl]) == ascnd) { // hunt up:
	for(;;) {
	  ju = jl + inc;
	  if (ju >= n-1) { ju = n-1; break; } // off end of table
	  else if ((x < xx[ju]) == ascnd) break; // found bracket
	  else {
	    jl = ju;
	    inc += inc;
	  }
	}
      } else { // hunt down:
	ju = jl;
	for (;;) {
	  jl = jl - inc;
	  if (jl <= 0) {jl = 0; break;} // off end of table
	  else if ((x >= xx[jl]) == ascnd) break; // found bracket
	  else {
	    ju = jl;
	    inc += inc;
	  }
	}
      }
    }
    while (ju - jl > 1) { // hunt is done, so begin final bisection
      jm = (ju+jl) >> 1;
      if ((x >= xx[jm]) == ascnd)
	jl = jm;
      else
	ju = jm;
    }
    cor = abs(jl - jsav) > dj ? 0 : 1; // decide whether to use hunt or
    // locate next time
    jsav = jl;
    return std::max(0,std::min(n-mmm,jl-((mmm-2)>>1)));
  }
  
  virtual Real rawinterp(const int jlo, const Real x) const = 0;
  // derived classes provide this as the actual interpolation method used.
};

template <typename Real>
struct Linear_interp : Base_interp<Real>
// piecewise linear interpolation object. Construct with x and y
// vectors, then call interp for interpolated values
{

  Linear_interp() : Base_interp<Real>() {}

  Linear_interp(const vector<Real> &xv, const vector<Real> &yv)
    : Base_interp<Real>(xv,yv,2) {}


  Real rawinterp(const int j, const Real x) const {
    Real* xptr = this->xx;
    Real* yptr = this->yy;
    
    if(xptr[j] == xptr[j+1]) return yptr[j]; // table defective, but recover
    else return yptr[j] + ((x-xptr[j])/(xptr[j+1]-xptr[j]))*(yptr[j+1]-yptr[j]);
  }
};

template <typename Real>
struct Bilinear_interp
// object for bilinear interpolation on a matrix. Construct with a
// vector of x1 values, a vector of x2 values, and a matrix of
// tabulated function values yij. Then call interp for interpolated
// values.
{
  int m, n;
  // we need dummy 1 dim interp objects for their locate and hunt functions
  Linear_interp<Real> x1terp, x2terp;
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> y;

  Bilinear_interp() { }

  Bilinear_interp(const vector<Real> &x1v,
		  const vector<Real> &x2v,
		  const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> &ym)
    : m(x1v.size()), n(x2v.size()), x1terp(x1v,x1v), x2terp(x2v,x2v), y(ym) {}
  
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
