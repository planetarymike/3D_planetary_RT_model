//influence.h

// functions to compute influence from one region of the atmosphere to
// another, including Holstein T and G functions

#ifndef __INFLUENCE_H
#define __INFLUENCE_H

#include <cmath>
#include <vector>
#include <cassert>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
//#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
//problem with boost redefines this ^^^ ???

using std::exp;
using std::log;
using std::vector;
using boost::math::quadrature::exp_sinh;
using boost::math::quadrature::tanh_sinh;
using boost::math::interpolators::cardinal_cubic_b_spline;

struct influence {
  virtual double Tint(const double /*tau*/) { return 0; };
  virtual double T(const double /*tau*/) { return 0; };
};

struct holstein_integrand {
  int order;
  double tau;

  holstein_integrand() { }

  holstein_integrand(int orderr, double tauu) :
    order(orderr), tau(tauu) { }
  
  double operator()(double x) const {
    double phi = exp(-x*x);
    
    assert((order == 1 || order == 2));

    double retval;
    if (order == 1)
      retval = phi*exp(-tau*phi);
    if (order == 2)
      retval = phi*phi*exp(-tau*phi);

    return retval;
  }    

};

struct holstein_exact {
  double integration_constant;

  holstein_exact() {
    integration_constant = M_2_SQRTPI; // 2/sqrt(pi)
  }

  double result(double tau, int order = 1) const {
    holstein_integrand integrand;
    exp_sinh<double> integrator;

    integrand.tau = tau;
    integrand.order = order;
    
    double result = integrator.integrate(integrand);

    return integration_constant*result;
  }

  double operator()(double tau) const {
    return this->result(tau, 1);
  }
  
};

//make faster somehow?? startup takes ~0.5s
//also, still not sure how to prove this integral is convergent
struct holstein_T_integral_exact {
  holstein_exact exact;
  tanh_sinh<double> integrator;

  holstein_T_integral_exact() { }

  double operator()(double tau) {
    
    double result = integrator.integrate(exact, 0.0, tau);

    return result;
  }

};

struct holstein_approx : influence {
  holstein_exact exact;
  holstein_T_integral_exact Tint_exact;
  
  double taumin;
  double taumax;
  double ntau;

  vector<double> logtau;
  vector<double> logG;
  vector<double> logT;
  vector<double> logTint;

  cardinal_cubic_b_spline<double> loglogGspline;
  cardinal_cubic_b_spline<double> loglogTspline;
  cardinal_cubic_b_spline<double> loglogTintspline;

  holstein_approx(double tauminn=1e-6,
		  double taumaxx=1e6,
		  int ntauu=100) :
    taumin(tauminn), taumax(taumaxx), ntau(ntauu) {
    
    double logtaumin = log(taumin);
    double logtaustep = ( log(taumax) - logtaumin ) / ( ntau - 1 );

    for (int itau = 0; itau < ntau; itau++) {
      logtau.push_back(  logtaumin + itau * logtaustep);
      logG.push_back(    log( exact.result(       exp( logtau[itau] ), 2)));
      logT.push_back(    log( exact.result(       exp( logtau[itau] ), 1)));
      logTint.push_back( log(   Tint_exact(       exp( logtau[itau] )   )));
    }

    loglogGspline = cardinal_cubic_b_spline<double>(logG.begin(),
						    logG.end(),
						    logtaumin,
						    logtaustep);
    loglogTspline = cardinal_cubic_b_spline<double>(logT.begin(),
						    logT.end(),
						    logtaumin,
						    logtaustep);
    loglogTintspline = cardinal_cubic_b_spline<double>(logTint.begin(),
    						       logTint.end(),
    						       logtaumin,
    						       logtaustep);

    
  }



  double Tint(const double tau) {
    if (tau < taumin) {
      return 0.0;
    } else if (tau > taumax) {
      return exp(loglogTintspline(log(taumax)));;
    } else {
      return exp(loglogTintspline(log(tau)));
    }
  }

  double T(const double tau) {
    if (tau < taumin) {
      return 1.0;
    } else if (tau > taumax) {
      return 0.0;
    } else {
      return exp(loglogTspline(log(tau)));
    }
  }

  double G(const double tau) {
    if (tau < taumin) {
      return M_SQRT1_2; // 1/sqrt(2)
    } else if (tau > taumax) {
      return 0.0;
    } else {
      return exp(loglogGspline(log(tau)));
    }
  }

  
};







#endif
