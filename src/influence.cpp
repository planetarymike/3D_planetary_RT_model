//influence.cpp -- functions to compute influence from one region of
// the atmosphere to another, including Holstein T and G functions

#include "influence.hpp"
using std::vector;

#include <cassert>

#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>

using std::exp;
using std::log;
using boost::math::quadrature::exp_sinh;
using boost::math::quadrature::tanh_sinh;

vector<Real> influence::Tintint(vector<Real> tau) const {
  vector<Real> retval;
  for(auto&& t: tau)
    retval.push_back(Tintint(t));
  return retval;
}

vector<Real> influence::Tint(vector<Real> tau) const {
  vector<Real> retval;
  for(auto&& t: tau)
    retval.push_back(Tint(t));
  return retval;
}

vector<Real> influence::T(vector<Real> tau) const {
  vector<Real> retval;
  for(auto&& t: tau)
    retval.push_back(T(t));
  return retval;
}




holstein_integrand::holstein_integrand() { }

holstein_integrand::holstein_integrand(int orderr, Real tauu) :
  order(orderr), tau(tauu) { }
  
Real holstein_integrand::operator()(Real x) const {
  Real phi = exp(-x*x);
    
  assert((order == 1 || order == 2));

  Real retval;
  if (order == 1)
    retval = phi*exp(-tau*phi);
  if (order == 2)
    retval = phi*phi*exp(-tau*phi);

  return retval;
}    










Real holstein_exact::result(Real tau, int order) const {
  holstein_integrand integrand(order,tau);
  exp_sinh<Real> integrator;

  Real result = integrator.integrate(integrand);

  return integration_constant*result;
}

Real holstein_exact::operator()(Real tau) const {
  return this->result(tau, 1);
}




Real holstein_T_integral_exact::operator()(Real tau) const {
  holstein_exact exact;
  tanh_sinh<Real> integrator;
  
  //make faster somehow?? startup takes ~0.5s
  Real result = integrator.integrate(exact, 0.0, tau, ABS);
  
  return result;
}







holstein_approx::holstein_approx(Real tauminn, Real taumaxx) :
  taumin(tauminn), taumax(taumaxx)
{
  
  logtaumin = log(taumin);
  logtaustep = ( log(taumax) - logtaumin ) / ( ntau - 1 );
  invlogtaustep = 1.0/logtaustep;
    
  for (int itau = 0; itau < ntau; itau++) {
    logtau[itau]  = logtaumin + itau * logtaustep;
    logG[itau]    =     log( exact.result(       exp( logtau[itau] ), 2));
    logT[itau]    =     log( exact.result(       exp( logtau[itau] ), 1));
    logTint[itau] =     log(   Tint_exact(       exp( logtau[itau] )   ));
  }

  loglogGspline = cardinal_cubic_b_spline<Real>(logG,
						ntau,
						logtaumin,
						logtaustep);
  loglogTspline = cardinal_cubic_b_spline<Real>(logT,
						ntau,
						logtaumin,
						logtaustep);
  loglogTintspline = cardinal_cubic_b_spline<Real>(logTint,
						   ntau,
						   logtaumin,
						   logtaustep);
}

// CUDA_CALLABLE_MEMBER
// int holstein_approx::get_lerp_loc(Real logtau_target) const {
//   int itau = (int) ((logtau_target-logtaumin)*invlogtaustep);
//   if (logtau_target<logtau[itau])
//     itau--;
//   if (logtau[itau+1]<=logtau_target)
//     itau++;
//   assert(0<=itau && itau<ntau && "itau must be in range");
//   assert(logtau[itau]<=logtau_target && logtau_target < logtau[itau+1] && "target tau must be in this increment");
//   return itau;
// }
  

Real holstein_approx::Tint(const Real tau) const {
  assert(tau<taumax && "tau must be in the simulated range.");

  if (tau < taumin) {
    return 0.0;
  }  else {
    return exp(loglogTintspline(log(tau)));
  }
}

CUDA_CALLABLE_MEMBER
Real holstein_approx::Tint_lerp(const Real tau) const {
  assert(tau<taumax && "tau must be in the simulated range.");
    
  if (tau < taumin) {
    return 0.0;
  }  else {
    Real logtau_target = log(tau);
    int itau = get_lerp_loc(logtau_target);
    Real lerp = logTint[itau]+(logTint[itau+1]-logTint[itau])*(logtau_target-logtau[itau])/(logtau[itau+1]-logtau[itau]);
    return exp(lerp);
  }
}

  
Real holstein_approx::T(const Real tau) const {
  assert(tau<taumax && "tau must be in the simulated range.");

  if (tau < taumin) {
    return 1.0;
  } else {
    return exp(loglogTspline(log(tau)));
  }
}

CUDA_CALLABLE_MEMBER
Real holstein_approx::T_lerp(const Real tau) const {
  assert(tau<taumax && "tau must be in the simulated range.");
    
  if (tau < taumin) {
    return 1.0;
  }  else {
    Real logtau_target = log(tau);
    int itau = get_lerp_loc(logtau_target);
    Real lerp = logT[itau]+(logT[itau+1]-logT[itau])*(logtau_target-logtau[itau])/(logtau[itau+1]-logtau[itau]);
      
    return exp(lerp);
  }
}


  
Real holstein_approx::G(const Real tau) const {
  assert(tau<taumax && "tau must be in the simulated range.");

  if (tau < taumin) {
    return M_SQRT1_2; // 1/sqrt(2)
  } else {
    return exp(loglogGspline(log(tau)));
  }
}
