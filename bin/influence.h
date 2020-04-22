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
//#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
//problem with boost redefines this ^^^ ???

using std::exp;
using std::log;
using std::vector;
using std::string;
using boost::math::quadrature::exp_sinh;
using boost::math::quadrature::tanh_sinh;
using boost::math::interpolators::cardinal_cubic_b_spline;

struct influence {
  virtual double Tintint(const double /*tau*/) { return 0; };
  vector<double> Tintint(vector<double> tau) {
    vector<double> retval;
    for(auto&& t: tau)
      retval.push_back(Tintint(t));
    return retval;
  }

  virtual double Tintabs(const double /*tau*/, const double /*abs*/) { return 0; };

  virtual double Tint(const double /*tau*/) { return 0; };
  vector<double> Tint(vector<double> tau) {
    vector<double> retval;
    for(auto&& t: tau)
      retval.push_back(Tint(t));
    return retval;
  }

  virtual double T(const double /*tau*/) { return 0; };
  vector<double> T(vector<double> tau) {
    vector<double> retval;
    for(auto&& t: tau)
      retval.push_back(T(t));
    return retval;
  }
};

struct holstein_integrand {
  int order;
  double tau;
  double abs;

  holstein_integrand() { }

  holstein_integrand(int orderr, double tauu, double abss) :
    order(orderr), tau(tauu), abs(abss) { }
  
  double operator()(double x) const {
    double phi = exp(-x*x);
    
    assert((order == 1 || order == 2));

    double retval;
    if (order == 1)
      retval = phi*exp(-tau*phi-tau*abs);
    if (order == 2)
      retval = phi*phi*exp(-tau*phi-tau*abs);

    return retval;
  }    

};

struct holstein_exact {
  double integration_constant;
  double abs;

  holstein_exact(double abss = 0.0)
    : abs(abss)
  {
    integration_constant = M_2_SQRTPI; // 2/sqrt(pi)
  }

  double result(double tau, int order = 1) const {
    holstein_integrand integrand(order,tau,abs);
    exp_sinh<double> integrator;

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

  holstein_T_integral_exact() { }

  double operator()(double tau, double abs = 0.0) const {
    holstein_exact exact(abs);
    tanh_sinh<double> integrator;

    double result = integrator.integrate(exact, 0.0, tau, 1e-6);

    return result;
  }

};

// struct holstein_Tint_integral_exact {

//   holstein_Tint_integral_exact() { }

//   double operator()(double tau) const {
//     holstein_T_integral_exact exact;
//     tanh_sinh<double> integrator;
    
//     double result = integrator.integrate(exact, 0.0, tau, 1e-3);

//     return result;
//   }

// };


// struct holstein_T_integral_with_absorption_interp {
//   vector<double> logtau_vec, tau_vec;
//   vector<double> logabs_vec, abs_vec;
//   Bilinear_interp interp;

//   holstein_T_integral_with_absorption_interp()
//   {
//     const int ntaupts=100;
//     const double logtaumin = 1e-3;
//     const double logtaumax = 1e3;
//     const double logtaustep = (logtaumax-logtaumin)/(ntaupts-1);
//     for (int itau = 0; itau < ntaupts; itau++) {
//       logtau_vec.push_back(logtaumin+itau*logtaustep);
//       tau_vec.push_back(exp(logtau_vec[itau]));
//     }

//     const int nabspts=40;
//     const double logabsmin = 1e-5;
//     const double logabsmax = 1e2;
//     const double logabsstep = (logabsmax-logabsmin)/(nabspts-1);
//     for (int iabs = 0; iabs < nabspts; iabs++) {
//       logabs_vec.push_back(logabsmin+iabs*logabsstep);
//       abs_vec.push_back(exp(logabs_vec[iabs]));
//     }

//     MatrixXd logHol;
//     logHol.resize(ntaupts,nabspts);

//     double tau,abs;
// #pragma omp parallel for private(tau,abs) shared(logHol) default(none)
//     for (int itau = 0; itau < ntaupts; itau++) {
//       tau = tau_vec[itau];
	
//       for (int iabs = 0; iabs < nabspts; iabs++) {
// 	abs = abs_vec[iabs];

// 	holstein_exact exact(abs);
// 	tanh_sinh<double> integrator;
	
// 	double result = integrator.integrate(exact, 0.0, tau, 1e-6);
// 	logHol(itau,iabs) = result;
//       }
//     }
    
//     interp = Bilinear_interp(logtau_vec, logabs_vec, logHol);
//   }
  
//   double operator()(const double tau, const double abs) {
//     assert(abs < abs_vec.back() && "absorption must be less than max absorption.");
//     assert(tau < tau_vec.back() && "optical depth must be less than max tau");

//     double logtau=log(tau);
//     double logabs=log(abs);
    
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
  
  double taumin;
  double taumax;
  double ntau;

  vector<double> logtau;
  vector<double> logG;
  vector<double> logT;
  vector<double> logTint;
  //  vector<double> logTintint;

  cardinal_cubic_b_spline<double> loglogGspline;
  cardinal_cubic_b_spline<double> loglogTspline;
  cardinal_cubic_b_spline<double> loglogTintspline;
  //  cardinal_cubic_b_spline<double> loglogTintintspline;

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
      //      logTintint.push_back( log(   Tintint_exact(       exp( logtau[itau] )   )));
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
    // loglogTintintspline = cardinal_cubic_b_spline<double>(logTint.begin(),
    // 							  logTint.end(),
    // 							  logtaumin,
    // 							  logtaustep);

    
  }



  // double Tintint(const double tau) const {
  //   assert(tau<taumax && "tau must be in the simulated range.");

  //   if (tau < taumin) {
  //     return 0.0;
  //   } else {
  //     return exp(loglogTintspline(log(tau)));
  //   }
  // }

  // double Tintabs(const double tau, const double abs) const {
  //   return Tintabs_interp(tau,abs);
  // }

  double Tint(const double tau) const {
    assert(tau<taumax && "tau must be in the simulated range.");

    if (tau < taumin) {
      return 0.0;
    }  else {
      return exp(loglogTintspline(log(tau)));
    }
  }

  double T(const double tau) const {
    assert(tau<taumax && "tau must be in the simulated range.");

    if (tau < taumin) {
      return 1.0;
    } else {
      return exp(loglogTspline(log(tau)));
    }
  }

  double G(const double tau) const {
    assert(tau<taumax && "tau must be in the simulated range.");

    if (tau < taumin) {
      return M_SQRT1_2; // 1/sqrt(2)
    } else {
      return exp(loglogGspline(log(tau)));
    }
  }

  
};







#endif
