//influence.h

// functions to compute influence from one region of the atmosphere to
// another, including Holstein T and G functions

#ifndef __INFLUENCE_H
#define __INFLUENCE_H

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <interp.h>
#include <fstream>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include "cuda_compatibility.h"
//#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
//problem with boost redefines this ^^^ ???

using std::exp;
using std::log;
using std::vector;
using std::string;
using boost::math::quadrature::exp_sinh;
using boost::math::quadrature::tanh_sinh;
using boost::math::interpolators::cardinal_cubic_b_spline;

struct influence {\
  virtual Real Tintint(const Real /*tau*/) const { return 0; };
  vector<Real> Tintint(vector<Real> tau) const {
    vector<Real> retval;
    for(auto&& t: tau)
      retval.push_back(Tintint(t));
    return retval;
  }

  virtual Real Tintabs(const Real /*tau*/, const Real /*abs*/) const { return 0; };

  virtual Real Tint(const Real /*tau*/) const { return 0; };
  vector<Real> Tint(vector<Real> tau) const {
    vector<Real> retval;
    for(auto&& t: tau)
      retval.push_back(Tint(t));
    return retval;
  }

  virtual Real T(const Real /*tau*/) const { return 0; };
  vector<Real> T(vector<Real> tau) const {
    vector<Real> retval;
    for(auto&& t: tau)
      retval.push_back(T(t));
    return retval;
  }
};

struct holstein_integrand {
  int order;
  Real tau;
  Real abs;

  holstein_integrand() { }

  holstein_integrand(int orderr, Real tauu, Real abss) :
    order(orderr), tau(tauu), abs(abss) { }
  
  Real operator()(Real x) const {
    Real phi = exp(-x*x);
    
    assert((order == 1 || order == 2));

    Real retval;
    if (order == 1)
      retval = phi*exp(-tau*phi-tau*abs);
    if (order == 2)
      retval = phi*phi*exp(-tau*phi-tau*abs);

    return retval;
  }    

};

struct holstein_exact {
  Real integration_constant;
  Real abs;

  holstein_exact(Real abss = 0.0)
    : abs(abss)
  {
    integration_constant = M_2_SQRTPI; // 2/sqrt(pi)
  }

  Real result(Real tau, int order = 1) const {
    holstein_integrand integrand(order,tau,abs);
    exp_sinh<Real> integrator;

    Real result = integrator.integrate(integrand);

    return integration_constant*result;
  }

  Real operator()(Real tau) const {
    return this->result(tau, 1);
  }
  
};

//make faster somehow?? startup takes ~0.5s
//also, still not sure how to prove this integral is convergent
struct holstein_T_integral_exact {

  holstein_T_integral_exact() { }

  Real operator()(Real tau, Real abs = 0.0) const {
    holstein_exact exact(abs);
    tanh_sinh<Real> integrator;

    Real result = integrator.integrate(exact, 0.0, tau, ABS);

    return result;
  }

};

// struct holstein_Tint_integral_exact {

//   holstein_Tint_integral_exact() { }

//   Real operator()(Real tau) const {
//     holstein_T_integral_exact exact;
//     tanh_sinh<Real> integrator;
    
//     Real result = integrator.integrate(exact, 0.0, tau, ABS);

//     return result;
//   }

// };


// struct holstein_T_integral_with_absorption_interp {
//   vector<Real> logtau_vec, tau_vec;
//   vector<Real> logabs_vec, abs_vec;
//   Bilinear_interp interp;

//   holstein_T_integral_with_absorption_interp()
//   {
//     const int ntaupts=100;
//     const Real logtaumin = 1e-3;
//     const Real logtaumax = 1e3;
//     const Real logtaustep = (logtaumax-logtaumin)/(ntaupts-1);
//     for (int itau = 0; itau < ntaupts; itau++) {
//       logtau_vec.push_back(logtaumin+itau*logtaustep);
//       tau_vec.push_back(exp(logtau_vec[itau]));
//     }

//     const int nabspts=40;
//     const Real logabsmin = 1e-5;
//     const Real logabsmax = 1e2;
//     const Real logabsstep = (logabsmax-logabsmin)/(nabspts-1);
//     for (int iabs = 0; iabs < nabspts; iabs++) {
//       logabs_vec.push_back(logabsmin+iabs*logabsstep);
//       abs_vec.push_back(exp(logabs_vec[iabs]));
//     }

//     MatrixX logHol;
//     logHol.resize(ntaupts,nabspts);

//     Real tau,abs;
// #pragma omp parallel for private(tau,abs) shared(logHol) default(none)
//     for (int itau = 0; itau < ntaupts; itau++) {
//       tau = tau_vec[itau];
	
//       for (int iabs = 0; iabs < nabspts; iabs++) {
// 	abs = abs_vec[iabs];

// 	holstein_exact exact(abs);
// 	tanh_sinh<Real> integrator;
	
// 	Real result = integrator.integrate(exact, 0.0, tau, ABS);
// 	logHol(itau,iabs) = result;
//       }
//     }
    
//     interp = Bilinear_interp(logtau_vec, logabs_vec, logHol);
//   }
  
//   Real operator()(const Real tau, const Real abs) {
//     assert(abs < abs_vec.back() && "absorption must be less than max absorption.");
//     assert(tau < tau_vec.back() && "optical depth must be less than max tau");

//     Real logtau=log(tau);
//     Real logabs=log(abs);
    
//     if (logabs < logabs_vec[0])
//       logabs = logabs_vec[0];
//     if (logtau < logtau_vec[0]) 
//       logtau = logtau_vec[0];
    
//     return exp(interp.interp(logtau,logabs));
//   }
// };

struct holstein_approx : influence {
  holstein_exact exact;
  holstein_T_integral_exact Tint_exact;
  //  holstein_Tint_integral_exact Tintint_exact;
  //  holstein_T_integral_with_absorption_interp Tintabs_interp;
  
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
		  Real taumaxx=1e6) :
    taumin(tauminn), taumax(taumaxx) {
    
    logtaumin = log(taumin);
    logtaustep = ( log(taumax) - logtaumin ) / ( ntau - 1 );
    invlogtaustep = 1.0/logtaustep;
    
    for (int itau = 0; itau < ntau; itau++) {
      logtau[itau]  = logtaumin + itau * logtaustep;
      logG[itau]    =     log( exact.result(       exp( logtau[itau] ), 2));
      logT[itau]    =     log( exact.result(       exp( logtau[itau] ), 1));
      logTint[itau] =     log(   Tint_exact(       exp( logtau[itau] )   ));
    }

    // for (int itau = 0; itau < ntau-1; itau++) {
    //   logTincrement[itau] = (logT[itau+1]-logT[itau])/(logtau[itau+1]-logtau[itau]);
    //   logTintincrement[itau] = (logTint[itau+1]-logTint[itau])/(logtau[itau+1]-logtau[itau]);
    // }
    
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

  Real Tint(const Real tau) const {
    assert(tau<taumax && "tau must be in the simulated range.");

    if (tau < taumin) {
      return 0.0;
    }  else {
      return exp(loglogTintspline(log(tau)));
    }
  }

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
  
  CUDA_CALLABLE_MEMBER
  Real Tint_lerp(const Real tau) const {
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

  
  Real T(const Real tau) const {
    assert(tau<taumax && "tau must be in the simulated range.");

    if (tau < taumin) {
      return 1.0;
    } else {
      return exp(loglogTspline(log(tau)));
    }
  }

  CUDA_CALLABLE_MEMBER
  Real T_lerp(const Real tau) const {
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


  
  Real G(const Real tau) const {
    assert(tau<taumax && "tau must be in the simulated range.");

    if (tau < taumin) {
      return M_SQRT1_2; // 1/sqrt(2)
    } else {
      return exp(loglogGspline(log(tau)));
    }
  }

  
};







#endif
