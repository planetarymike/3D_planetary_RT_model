#include "interp.hpp"

Base_interp::Base_interp(const vector<Real> &x, const vector<Real> &y, const int m)
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

Base_interp::Base_interp(const Base_interp &B) {
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

Base_interp& Base_interp::operator= (const Base_interp &B) {
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

Base_interp::~Base_interp() {
  if (n>0) {
    delete [] xx;
    delete [] yy;
  }
}

Real Base_interp::operator()(const Real x) const 
// Given a value x, return an interpolated value, using data
// pointed to by xx and yy.
{
  int jlo = cor ? hunt(x) : locate (x);
  return rawinterp(jlo, x);
}

int Base_interp::index(const Real x) const
// Given a value x, return an interpolated value, using data
// pointed to by xx and yy.
{
  int jlo = cor ? hunt(x) : locate(x);
  return jlo;
}




int Base_interp::locate(const Real x) const
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


int Base_interp::hunt(const Real x) const
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
  if (jl < 0 || jl > n-1) { // input guess not useful. go directly to
			    // bisection
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

Linear_interp::Linear_interp(const vector<Real> &xv, const vector<Real> &yv)
  : Base_interp(xv,yv,2) {}



Real Linear_interp::rawinterp(const int j, const Real x) const {
  if(xx[j] == xx[j+1]) return yy[j]; // table defective, but recover
  else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
}

Bilinear_interp::Bilinear_interp(const vector<Real> &x1v,
				 const vector<Real> &x2v,
				 const MatrixX &ym)
  : m(x1v.size()), n(x2v.size()), x1terp(x1v,x1v), x2terp(x2v,x2v) { y = ym; }

Real Bilinear_interp::interp(const Real x1p, const Real x2p) const {
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

