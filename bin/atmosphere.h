//atmosphere.h -- routines to compute CO2 and hydrogen number density.

#ifndef __ATMOSPHERE_H_
#define __ATMOSPHERE_H_

#include "constants.h"
#include "atmo_vec.h"
#include "interp.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>
#include <type_traits>
#include <boost/numeric/odeint/integrate/integrate.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/tools/roots.hpp>

using std::exp;
using std::log;
using boost::math::gamma_p;
using std::vector;
using boost::numeric::odeint::runge_kutta4;
using boost::numeric::odeint::integrate;
using boost::numeric::odeint::integrate_const;
using boost::math::interpolators::cardinal_cubic_b_spline;


//temperature

//generic temperature class
struct temperature {
  double T_exo;
  double T;
  double Tprime;
  
  virtual void get(const double &r) { };
};

struct krasnopolsky_temperature : public temperature {
  // computes the analytic thermospheric temperature as given by Krasnopolsky (2002)
  
  double T_tropo;
  double r_tropo;
  double shape_parameter;

  krasnopolsky_temperature() { }
  
  krasnopolsky_temperature(double T_exoo,
			   double T_tropoo = 125.0,
			   double r_tropoo = rMars + 90e5,
			   double shape_parameterr = 11.4) {
    T_exo = T_exoo; 
    T_tropo = T_tropoo; 
    r_tropo = r_tropoo; 
    shape_parameter = shape_parameterr;
  }

  void get(const double &r) {
    const double rdiff = r - r_tropo;
    T = T_exo - (T_exo - T_tropo)*exp(-rdiff*rdiff*1e-10/(shape_parameter*T_exo));
    Tprime = ( T_exo - T ) * ( 2*rdiff*1e-10 / (shape_parameter*T_exo) );
  }

};

//diffusion coefficients
struct diffusion_coefs {
  double DH; // diffusion coefficient of hydrogen through CO2
  double KK; // eddy diffusion coefficient
  
  void get(const double &r, const double &T, const double &Texo, const double &nCO2) {
    const double DH0 = 8.4e17;// cm^2 s^-1
    const double s = 0.6;
    DH = std::pow(T,s) * DH0/nCO2;
    KK = 1.2e12 * std::sqrt(Texo/nCO2); 
  }
  
};

// define the differential equation used to solve for the CO2 and H
// number densities in the themosphere
struct thermosphere_diffeq {
  diffusion_coefs diff;
  temperature &temp;
  double &H_escape_flux;
  double &rexo;

  thermosphere_diffeq(temperature &tempp, double &H_escape_fluxx, double &rexoo)
    : temp(tempp), H_escape_flux(H_escape_fluxx), rexo(rexoo) { }

  // x[0] = log(nCO2), x[1] = log(nH)
  void operator()( const vector<double> &x , vector<double> &dxdr , const double &r ) {
    temp.get(r);
    diff.get(r, temp.T, temp.T_exo, exp(x[0]) );

    double Hninv = G*mMars*mCO2/(kB*temp.T*r*r)+temp.Tprime/temp.T;
    double alpha = -0.25;
    double HHinv = G*mMars*mH/(kB*temp.T*r*r)+(1+alpha)*temp.Tprime/temp.T;
    
    dxdr[0] = -Hninv;
    dxdr[1] = -1./(diff.DH+diff.KK)*( H_escape_flux * ( rexo*rexo /r /r ) / exp(x[1])
				      + (diff.DH*HHinv + diff.KK*Hninv) );
    
  }

};

struct push_back_quantities
{
  vector< vector<double>* > vec_ptrs;

  template<bool...> struct bool_pack{};
  template<class... Ts>
  using conjunction = std::is_same<bool_pack<true,Ts::value...>, bool_pack<Ts::value..., true>>;
  template<typename... Ts>
  using AllVecs = typename std::enable_if<conjunction<std::is_convertible<Ts, vector<double>*>...>::value>::type;
  template<typename... Ts, typename = AllVecs<Ts...>>
  push_back_quantities(Ts... v0) {
    vector<vector<double>*> vecs = { v0... };
    for (auto v : vecs) {
      vec_ptrs.push_back(v);
    }
  }

  void operator()( const vector<double> &x , double r )
  {
    for(unsigned int i=0; i<x.size(); i++)
      (*(vec_ptrs[i])).push_back(x[i]);
    (*(vec_ptrs.back())).push_back( r );
  }
};


struct chamberlain_exosphere {
  double &rexo;
  double &Texo;
  double &nHexo;
  double lambdac;
  double effusion_velocity;
  double H_escape_flux;

  chamberlain_exosphere(double &rexoo, double &Texoo, double &nHexoo)
    : rexo(rexoo), Texo(Texoo), nHexo(nHexoo)
  {
    lambdac = G*mMars*mH/(kB*Texo*rexo);//chamberlain lambda @ rexo
    effusion_velocity = 0.5 * sqrt( 2.0*kB*Texo / (mH*pi) ) * (1.0 + lambdac) * exp(-lambdac); // effusion velocity
    H_escape_flux = nHexo * effusion_velocity;
  }

  double operator()(double r) {
    // computes hydrogen number density as a function of altitude above
    // the exobase, assuming a chamberlain exosphere w/o satellite particles

    if ((r-rexo)/rexo<1e-6) 
      //if we are very close to the exobase due to a rounding error,
      //fix the issue
      r = rexo;
    
    const double lambda = G*mMars*mH/(kB*Texo*r);//chamberlain lambda
    const double psione = lambda*lambda/(lambda+lambdac);

    //gamma_p = complementary normalized incomplete gamma function
    const double g1 = gamma_p(1.5, lambda);
    const double g1c = gamma_p(1.5, lambdac);
    const double g2 = gamma_p(1.5, lambda - psione);
  
    // calculate the fraction of the exobase density at this altitude:
    // this looks different than the formula in C&H b/c the incomplete
    // gamma function in the gsl is normalized; to get gamma(3/2,x) from
    // C&H we have to multiply g1 by Gamma(3/2). This enables us to pull
    // the complete Gamma function outside the parens.
    double frac = ( 1.0 + g1 - sqrt( 1.0 - lambda*lambda/(lambdac*lambdac) ) * exp(-psione) * ( 1.0 + g2 ) );  
    // normalize to the fraction at the exobase:
    double norm = (1.0 + g1c);
    frac /= norm;

    // multiply by the exobase density and return
    return nHexo*frac*exp(lambda-lambdac);
  }
};

struct atmosphere {
  double rmin;// cm, exobase altitude
  double rexo;// cm, minimum altitude in model atmosphere
  double rmax;// cm, max altitude in model atmosphere

  virtual double nH(const atmo_point pt) { return 0; };
  virtual vector<double> nH(const vector<atmo_point> pts) {
    vector<double> ret;
    ret.resize(pts.size());
    for(unsigned int i=0;i<pts.size();i++)
      ret[i] = nH(pts[i]);
    return ret;
  }
  virtual double nH(const double r) {
    //return subsolar densities if not overridden
    atmo_point p;
    p.rtp(r,0.,0.);
    return nH(p);
  }

  virtual double r_from_nH(const double nH) { return 0; }

  virtual double nCO2(const atmo_point pt) { return 0; };
  virtual vector<double> nCO2(const vector<atmo_point> pts) {
    vector<double> ret;
    ret.resize(pts.size());
    for(unsigned int i=0;i<pts.size();i++)
      ret[i] = nCO2(pts[i]);
    return ret;
  }
  virtual double nCO2(const double r) {
    //return subsolar densities if not overridden
    atmo_point p;
    p.rtp(r,0.,0.);
    return nCO2(p);
  }

  atmosphere(double rminn, double rexoo, double rmaxx)
    : rmin(rminn), rexo(rexoo), rmax(rmaxx) { }

  double s_null(const atmo_point pt) {
    return 0.0;
  }

  //really ought to refactor so cross section info is stored in a
  //totally seperate object
  virtual double sH_lya(const atmo_point pt) {
    return 0.0;
  }
  virtual double sH_lya(const double r) {
    //return subsolar densities if not overridden
    atmo_point p;
    p.rtp(r,0.,0.);
    return sH_lya(p);
  }
  
  virtual double sCO2_lya(const atmo_point pt) {
    return 0.0;
  }
  virtual double sCO2_lya(const double r) {
    //return subsolar densities if not overridden
    atmo_point p;
    p.rtp(r,0.,0.);
    return sCO2_lya(p);
  }

};


//now a structure to hold the atmosphere and interpolate when
//necessary
struct chamb_diff_1d : atmosphere {
  double nHexo;   // cm-3, H density at exobase
  double nCO2exo; // cm-3, CO2 density at exobase

  double rmindiffusion; // cm, minimum altitude to solve the diffusion equation 
  double nHrmindiffusion; // nH at this altitude
  double nCO2rmindiffusion;
  
  temperature &temp;
  chamberlain_exosphere exosphere;
  thermosphere_diffeq diffeq;

  //thermosphere interpolation object
  int n_thermosphere_steps;
  double thermosphere_step_r;
  vector<double> lognCO2thermosphere;
  vector<double> lognHthermosphere;
  vector<double> r_thermosphere;
  cardinal_cubic_b_spline<double> lognCO2_thermosphere_spline;
  Linear_interp invlognCO2_thermosphere;
  cardinal_cubic_b_spline<double> lognH_thermosphere_spline;
  Linear_interp invlognH_thermosphere;
  
  //exosphere interpolation
  int n_exosphere_steps;
  double exosphere_step_logr;
  vector<double> lognHexosphere;
  vector<double> logr_exosphere;
  cardinal_cubic_b_spline<double> lognH_exosphere_spline;
  Linear_interp invlognH_exosphere;

  chamb_diff_1d(double nHexoo, // a good number is 10^5-6
		double nCO2exoo, //a good number is 10^9 (?)
		temperature &tempp)
    : chamb_diff_1d(/*          rmin = */rMars + 80e5,
		    /*          rexo = */rMars + 200e5,
		    /*          rmax = */rMars + 50000e5,
		    /* rmindiffusion = */rMars + 120e5,
		    nHexoo,
		    nCO2exoo,
		    tempp)   { }

  chamb_diff_1d(double rminn,
		double rexoo,
		double rmaxx,
		double rmindiffusionn,
		double nHexoo, // a good number is 10^5-6
		double nCO2exoo, //a good number is 10^9 (?)
		temperature &tempp)
    : atmosphere(rminn,rexoo,rmaxx),
      nHexo(nHexoo),
      nCO2exo(nCO2exoo),
      rmindiffusion(rmindiffusionn),
      temp(tempp), 
      exosphere(rexo, temp.T_exo, nHexo),
      diffeq(tempp, exosphere.H_escape_flux, rexo)      
  {

    //integrate the differential equation to get the species densities
    //in the thermosphere
    vector<double> nexo(2);
    nexo[0] = log(nCO2exo);
    nexo[1] = log(nHexo);

    //use a constant stepper for easy interpolation
    runge_kutta4< vector<double> > stepper;
    n_thermosphere_steps = 20;
    thermosphere_step_r = -(rexo-rmin)/(n_thermosphere_steps-1.);
    integrate_const( stepper , diffeq ,
		     nexo , rexo , rmin , thermosphere_step_r,
		     push_back_quantities( &lognCO2thermosphere,
					   &lognHthermosphere,
					   &r_thermosphere ) );
    //interpolate the densities in the thermosphere
    lognCO2_thermosphere_spline = cardinal_cubic_b_spline<double>(lognCO2thermosphere.rbegin(),
								  lognCO2thermosphere.rend(),
								  rmin,
								  -thermosphere_step_r);
    invlognCO2_thermosphere = Linear_interp(lognCO2thermosphere,r_thermosphere);
    lognH_thermosphere_spline = cardinal_cubic_b_spline<double>(lognHthermosphere.rbegin(),
								lognHthermosphere.rend(),
								rmin,
								-thermosphere_step_r);
    invlognH_thermosphere = Linear_interp(lognHthermosphere,r_thermosphere);

    nHrmindiffusion = nH(rmindiffusion);
    nCO2rmindiffusion = nCO2(rmindiffusion);
    
    //now get the interpolation points in the exosphere
    n_exosphere_steps = 20;
    exosphere_step_logr = (log(rmax) - log(rexo))/(n_exosphere_steps - 1.);
    for (int iexo = 0; iexo < n_exosphere_steps; iexo++) {
      logr_exosphere.push_back( log(rexo) + iexo * exosphere_step_logr );
      lognHexosphere.push_back( log( exosphere( exp( logr_exosphere[iexo] ) ) ) );
    }
    lognH_exosphere_spline = cardinal_cubic_b_spline<double>(lognHexosphere.begin(),
							     lognHexosphere.end(),
							     log(rexo),
							     exosphere_step_logr);
    invlognH_exosphere = Linear_interp(lognHexosphere,logr_exosphere);

  }


  double nCO2(const double &r) {
    if (r>rexo)
      return 0.0;
    else {
      assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
      return exp(lognCO2_thermosphere_spline(r));
    }
  }

  double nCO2(const atmo_point pt) {
    return nCO2(pt.r);
  }


  double nH(const double &r) {
    if (r>=rexo)
      return exp(lognH_exosphere_spline(log(r)));
    else {
      if (r>=rmindiffusion)
	return exp(lognH_thermosphere_spline(r));
      else {
	assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
	return nHrmindiffusion/nCO2rmindiffusion*nCO2(r);
      }
    }
  }

  double nH(const atmo_point pt) {
    return nH(pt.r);
  }

  double sH_lya(const double r) {
    temp.get(r);
    return lyman_alpha_line_center_cross_section_coef/sqrt(temp.T);
  }    
  double sH_lya(const atmo_point pt) {
    return sH_lya(pt.r);
  }

  double sCO2_lya(const double r) {
    return CO2_lyman_alpha_absorption_cross_section;
  }
  double sCO2_lya(const atmo_point pt) {
    return sCO2_lya(pt.r);
  }


  double r_from_nH(double nHtarget) {
    if (nHtarget==nHexo) {
      return rexo;
    } else if (nHtarget<nHexo) {
      return exp(invlognH_exosphere(log(nHtarget)));
    } else if (nHtarget>nHrmindiffusion) {
      return invlognCO2_thermosphere(log(nHtarget*nCO2rmindiffusion/nHrmindiffusion));
    } else {
      return invlognH_thermosphere(log(nHtarget));
    }
  }
  

  double nCO2_exact(const double &r) {
    if (r>rexo)
      return 0.0;
    else {
      assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
      vector<double> lognCO2thermosphere_tmp;
      vector<double> lognHthermosphere_tmp;
      vector<double> r_thermosphere_tmp;

      vector<double> nexo(2);
      nexo[0] = log(nCO2exo);
      nexo[1] = log(nHexo);
      integrate( diffeq , nexo , rexo , r , thermosphere_step_r,
		 push_back_quantities( &lognCO2thermosphere_tmp,
				       &lognHthermosphere_tmp,
				       &r_thermosphere_tmp ));
      return exp(lognCO2thermosphere_tmp.back());
    }
  }


  double nH_exact(const double &r) {
    if (r>=rexo)
      return exosphere(r);
    else {
      assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
      if (r>=rmindiffusion) {
	vector<double> lognCO2thermosphere_tmp;
	vector<double> lognHthermosphere_tmp;
	vector<double> r_thermosphere_tmp;
	
	vector<double> nexo(2);
	nexo[0] = log(nCO2exo);
	nexo[1] = log(nHexo);
	integrate( diffeq , nexo , rexo , r , thermosphere_step_r,
		   push_back_quantities( &lognCO2thermosphere_tmp,
					 &lognHthermosphere_tmp, 
					 &r_thermosphere_tmp ));
	return  exp(lognHthermosphere_tmp.back());
      } else 
	return nH_exact(rmindiffusion)/nCO2_exact(rmindiffusion)*nCO2_exact(r);
    }
  }



  
};



#endif
