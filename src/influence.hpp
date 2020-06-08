//influence.hpp -- functions to compute influence from one region of
// the atmosphere to another, including Holstein T and G functions

#ifndef __INFLUENCE_H
#define __INFLUENCE_H

#include <cmath>
#include <vector>
using std::vector;

#include "Real.hpp"
#include "cuda_compatibility.hpp"
#include "interp.hpp"

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using boost::math::interpolators::cardinal_cubic_b_spline;


struct influence {
public:
  virtual Real Tintint(const Real /*tau*/) const { return 0; };
  vector<Real> Tintint(vector<Real> tau) const;

  virtual Real Tint(const Real /*tau*/) const { return 0; };
  vector<Real> Tint(vector<Real> tau) const;

  virtual Real T(const Real /*tau*/) const { return 0; };
  vector<Real> T(vector<Real> tau) const;
};










struct holstein_integrand {
  int order;
  Real tau;

  holstein_integrand();

  holstein_integrand(int orderr, Real tauu);
  
  Real operator()(Real x) const;
};











struct holstein_exact {
  constexpr static Real integration_constant = M_2_SQRTPI;

  Real result(Real tau, int order = 1) const;

  Real operator()(Real tau) const;  
};




struct holstein_T_integral_exact {

  Real operator()(Real tau) const;
};













struct holstein_approx : influence {
  holstein_exact exact;
  holstein_T_integral_exact Tint_exact;
  
  Real taumin;
  Real taumax;
  static const int ntau = 200;

  Real logtaumin;
  Real logtaustep;
  Real invlogtaustep;

  Real logtau[ntau];
  Real logG[ntau];
  Real logT[ntau];
  //Real logTincrement[ntau-1];
  Real logTint[ntau];
  //Real logTintincrement[ntau-1];

  cardinal_cubic_b_spline<Real> loglogGspline;
  cardinal_cubic_b_spline<Real> loglogTspline;
  cardinal_cubic_b_spline<Real> loglogTintspline;

  holstein_approx(Real tauminn=1e-6,
		  Real taumaxx=1e6);

  inline int get_lerp_loc(Real logtau_target) const {
    int itau = (int) ((logtau_target-logtaumin)*invlogtaustep);
    if (logtau_target<logtau[itau])
      itau--;
    if (logtau[itau+1]<=logtau_target)
      itau++;
    assert(0<=itau && itau<ntau && "itau must be in range");
    assert(logtau[itau]<=logtau_target && logtau_target < logtau[itau+1] && "target tau must be in this increment");
    return itau;
  }

    

  Real Tint(const Real tau) const;
  CUDA_CALLABLE_MEMBER
  Real Tint_lerp(const Real tau) const;
  
  Real T(const Real tau) const;
  CUDA_CALLABLE_MEMBER
  Real T_lerp(const Real tau) const;
  
  Real G(const Real tau) const;  
};







#endif
