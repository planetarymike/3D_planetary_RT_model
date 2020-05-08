//atmosphere.h -- routines to compute CO2 and hydrogen number density.

#ifndef __ATMOSPHERE_H_
#define __ATMOSPHERE_H_

#include "constants.h"
#include "atmo_vec.h"
#include "interp.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <type_traits>
#include <Eigen/Dense>
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
  Real T_exo;
  Real T;
  Real Tprime;
  
  virtual void get(const Real &r) { };
};

struct krasnopolsky_temperature : public temperature {
  // computes the analytic thermospheric temperature as given by Krasnopolsky (2002)
  
  Real T_tropo;
  Real r_tropo;
  Real shape_parameter;

  krasnopolsky_temperature(Real T_exoo = 200.0,
			   Real T_tropoo = 125.0,
			   Real r_tropoo = rMars + 90e5,
			   Real shape_parameterr = 11.4) {
    T_exo = T_exoo; 
    T_tropo = T_tropoo; 
    r_tropo = r_tropoo; 
    shape_parameter = shape_parameterr;
  }

  void get(const Real &r) {
    const Real rdiff = r - r_tropo;
    T = T_exo - (T_exo - T_tropo)*exp(-rdiff*rdiff*1e-10/(shape_parameter*T_exo));
    Tprime = ( T_exo - T ) * ( 2*rdiff*1e-10 / (shape_parameter*T_exo) );
  }

};

//diffusion coefficients
struct diffusion_coefs {
  Real DH; // diffusion coefficient of hydrogen through CO2
  Real KK; // eddy diffusion coefficient
  
  void get(const Real &r, const Real &T, const Real &Texo, const Real &nCO2) {
    const Real DH0 = 8.4e17;// cm^2 s^-1
    const Real s = 0.6;
    DH = std::pow(T,s) * DH0/nCO2;
    KK = 1.2e12 * std::sqrt(Texo/nCO2); 
  }
  
};

// define the differential equation used to solve for the CO2 and H
// number densities in the themosphere
struct thermosphere_diffeq {
  diffusion_coefs diff;
  temperature *temp;
  Real H_escape_flux;
  Real rexo;

  thermosphere_diffeq(temperature &tempp, Real &H_escape_fluxx, Real &rexoo)
    : temp(&tempp), H_escape_flux(H_escape_fluxx), rexo(rexoo) { }

  // x[0] = log(nCO2), x[1] = log(nH)
  void operator()( const vector<Real> &x , vector<Real> &dxdr , const Real &r ) {
    temp->get(r);
    diff.get(r, temp->T, temp->T_exo, exp(x[0]) );

    Real Hninv = G*mMars*mCO2/(kB*(temp->T)*r*r)+(temp->Tprime)/(temp->T);
    Real alpha = -0.25;
    Real HHinv = G*mMars*mH/(kB*temp->T*r*r)+(1+alpha)*(temp->Tprime)/(temp->T);
    
    dxdr[0] = -Hninv;
    dxdr[1] = -1./(diff.DH+diff.KK)*( H_escape_flux * ( rexo*rexo /r /r ) / exp(x[1])
				      + (diff.DH*HHinv + diff.KK*Hninv) );
    
  }

};

struct push_back_quantities
{
  vector< vector<Real>* > vec_ptrs;

  template<bool...> struct bool_pack{};
  template<class... Ts>
  using conjunction = std::is_same<bool_pack<true,Ts::value...>, bool_pack<Ts::value..., true>>;
  template<typename... Ts>
  using AllVecs = typename std::enable_if<conjunction<std::is_convertible<Ts, vector<Real>*>...>::value>::type;
  template<typename... Ts, typename = AllVecs<Ts...>>
  push_back_quantities(Ts... v0) {
    vector<vector<Real>*> vecs = { v0... };
    for (auto&& v : vecs) {
      vec_ptrs.push_back(v);
    }
  }

  void operator()( const vector<Real> &x , Real r )
  {
    for(unsigned int i=0; i<x.size(); i++)
      (*(vec_ptrs[i])).push_back(x[i]);
    (*(vec_ptrs.back())).push_back( r );
  }
};


struct chamberlain_exosphere {
  Real rexo;
  Real Texo;
  Real nHexo;
  Real lambdac;
  Real effusion_velocity;
  Real H_escape_flux;

  chamberlain_exosphere(Real &rexoo, Real &Texoo, Real &nHexoo)
    : rexo(rexoo), Texo(Texoo), nHexo(nHexoo)
  {
    lambdac = G*mMars*mH/(kB*Texo*rexo);//chamberlain lambda @ rexo
    effusion_velocity = 0.5 * sqrt( 2.0*kB*Texo / (mH*pi) ) * (1.0 + lambdac) * exp(-lambdac); // effusion velocity
    H_escape_flux = nHexo * effusion_velocity;
  }

  Real nH(Real r) {
    // computes hydrogen number density as a function of altitude above
    // the exobase, assuming a chamberlain exosphere w/o satellite particles

    if ((r-rexo)/rexo<ABS) 
      //if we are very close to the exobase due to a rounding error,
      //fix the issue
      r = rexo;
    
    const Real lambda = G*mMars*mH/(kB*Texo*r);//chamberlain lambda
    const Real psione = lambda*lambda/(lambda+lambdac);

    //gamma_p = complementary normalized incomplete gamma function
    const Real g1 = gamma_p(1.5, lambda);
    const Real g1c = gamma_p(1.5, lambdac);
    const Real g2 = gamma_p(1.5, lambda - psione);
  
    // calculate the fraction of the exobase density at this altitude:
    // this looks different than the formula in C&H b/c the incomplete
    // gamma function in the gsl is normalized; to get gamma(3/2,x) from
    // C&H we have to multiply g1 by Gamma(3/2). This enables us to pull
    // the complete Gamma function outside the parens.
    Real frac = ( 1.0 + g1 - sqrt( 1.0 - lambda*lambda/(lambdac*lambdac) ) * exp(-psione) * ( 1.0 + g2 ) );  
    // normalize to the fraction at the exobase:
    Real norm = (1.0 + g1c);
    frac /= norm;

    // multiply by the exobase density and return
    return nHexo*frac*exp(lambda-lambdac);
  }
  Real operator ()(Real r) {
    return this->nH(r);
  }


  template <typename T>
  struct nHfinder {
    T *parent;
    Real nHtarget;

    nHfinder(T *parentt, Real &nHtargett)
      : parent(parentt), nHtarget(nHtargett)
    { }

    Real operator()(Real r) {
      return parent->nH(r) - nHtarget;
    }
  };

  Real r(Real &nHtarget) {  //find r corresponding to a given nH
    assert(nHtarget<nHexo && "exosphere nH must be less than nHexo");


    using boost::math::tools::bracket_and_solve_root;
    using boost::math::tools::eps_tolerance;

    
    Real guess = kB*Texo/G/mMars/mH*log(nHexo/nHtarget)+rexo;
    Real factor = 2;

    const boost::uintmax_t maxit = 20;
    boost::uintmax_t it = maxit;      
    bool is_rising = false;
    int get_digits = 8;
    eps_tolerance<Real> tol(get_digits);

    nHfinder<chamberlain_exosphere> find(this,nHtarget);
    std::pair<Real, Real> r = bracket_and_solve_root(find,
							 guess, factor, is_rising, tol, it);


    assert(it < maxit && "we must find the value we want in less than the maximum number of iterations.");
    // if (it >= maxit)
    //   {
    // 	std::cout << "Unable to locate solution in " << maxit << " iterations:"
    // 	  " Current best guess is between " << r.first << " and " << r.second << std::endl;
    //   }
    // else
    //   {
    // 	std::cout << "Converged after " << it << " (from maximum of " << maxit << " iterations)." << std::endl;
    // 	std::cout << "Converged to " << r.first + (r.second - r.first)/2 << std::endl;
    //   }
    

    return r.first + (r.second - r.first)/2;
  }


  
};

struct atmosphere {
  Real rmin;// cm, exobase altitude
  Real rexo;// cm, minimum altitude in model atmosphere
  Real rmax;// cm, max altitude in model atmosphere

  virtual Real nH(const atmo_point pt) { return 0; };
  virtual vector<Real> nH(const vector<atmo_point> pts) {
    vector<Real> ret;
    ret.resize(pts.size());
    for(unsigned int i=0;i<pts.size();i++)
      ret[i] = nH(pts[i]);
    return ret;
  }
  virtual Real nH(const Real r) {
    //return subsolar densities if not overridden
    atmo_point p;
    p.rtp(r,0.,0.);
    return nH(p);
  }

  virtual Real r_from_nH(const Real nH) { return 0; }

  virtual Real nCO2(const atmo_point pt) { return 0; };
  virtual vector<Real> nCO2(const vector<atmo_point> pts) {
    vector<Real> ret;
    ret.resize(pts.size());
    for(unsigned int i=0;i<pts.size();i++)
      ret[i] = nCO2(pts[i]);
    return ret;
  }
  virtual Real nCO2(const Real r) {
    //return subsolar densities if not overridden
    atmo_point p;
    p.rtp(r,0.,0.);
    return nCO2(p);
  }

  atmosphere(Real rminn, Real rexoo, Real rmaxx)
    : rmin(rminn), rexo(rexoo), rmax(rmaxx) { }

  Real s_null(const atmo_point pt) {
    return 0.0;
  }

  //really ought to refactor so cross section info is stored in a
  //totally seperate object
  virtual Real sH_lya(const atmo_point pt) {
    return 0.0;
  }
  virtual Real sH_lya(const Real r) {
    //return subsolar densities if not overridden
    atmo_point p;
    p.rtp(r,0.,0.);
    return sH_lya(p);
  }
  
  virtual Real sCO2_lya(const atmo_point pt) {
    return 0.0;
  }
  virtual Real sCO2_lya(const Real r) {
    //return subsolar densities if not overridden
    atmo_point p;
    p.rtp(r,0.,0.);
    return sCO2_lya(p);
  }

};


//now a structure to hold the atmosphere and interpolate when
//necessary
struct chamb_diff_1d : atmosphere {
  Real nHexo;   // cm-3, H density at exobase
  Real nCO2exo; // cm-3, CO2 density at exobase

  Real rmindiffusion; // cm, minimum altitude to solve the diffusion equation 
  Real nHrmindiffusion; // nH at this altitude
  Real nCO2rmindiffusion;
  
  temperature *temp;
  chamberlain_exosphere exosphere;
  thermosphere_diffeq diffeq;

  //thermosphere interpolation object
  int n_thermosphere_steps;
  Real thermosphere_step_r;
  vector<Real> lognCO2thermosphere;
  vector<Real> lognHthermosphere;
  vector<Real> r_thermosphere;
  cardinal_cubic_b_spline<Real> lognCO2_thermosphere_spline;
  Linear_interp invlognCO2_thermosphere;
  cardinal_cubic_b_spline<Real> lognH_thermosphere_spline;
  Linear_interp invlognH_thermosphere;
  
  //exosphere interpolation
  int n_exosphere_steps;
  Real exosphere_step_logr;
  vector<Real> lognHexosphere;
  vector<Real> logr_exosphere;
  cardinal_cubic_b_spline<Real> lognH_exosphere_spline;
  Linear_interp invlognH_exosphere;

  chamb_diff_1d(Real nHexoo, // a good number is 10^5-6
		Real nCO2exoo, //a good number is 10^9 (?)
		temperature &tempp)
    : chamb_diff_1d(/*          rmin = */rMars + 80e5,
		    /*          rexo = */rMars + 200e5,
		    /*         nHmin = */10,
		    /* rmindiffusion = */rMars + 120e5,
		    nHexoo,
		    nCO2exoo,
		    tempp)   { }

  chamb_diff_1d(Real rminn,
		Real rexoo,
		Real nHmin,
		Real rmindiffusionn,
		Real nHexoo, // a good number is 10^5-6
		Real nCO2exoo, //a good number is 10^9 (?)
		temperature &tempp)
    : atmosphere(rminn,rexoo,-1),
      nHexo(nHexoo),
      nCO2exo(nCO2exoo),
      rmindiffusion(rmindiffusionn),
      temp(&tempp), 
      exosphere(rexo, temp->T_exo, nHexo),
      diffeq(tempp, exosphere.H_escape_flux, rexo)      
  {

    //set the max altitude by finding the density at which the exosphere = nHmin
    rmax = exosphere.r(nHmin);
    
    //integrate the differential equation to get the species densities
    //in the thermosphere
    vector<Real> nexo(2);
    nexo[0] = log(nCO2exo);
    nexo[1] = log(nHexo);

    //use a constant stepper for easy interpolation
    runge_kutta4< vector<Real> > stepper;
    n_thermosphere_steps = 20;
    thermosphere_step_r = -(rexo-rmin)/(n_thermosphere_steps-1.);
    integrate_const( stepper , diffeq ,
		     nexo , rexo , rmin , thermosphere_step_r,
		     push_back_quantities( &lognCO2thermosphere,
					   &lognHthermosphere,
					   &r_thermosphere ) );
    //interpolate the densities in the thermosphere
    lognCO2_thermosphere_spline = cardinal_cubic_b_spline<Real>(lognCO2thermosphere.rbegin(),
								  lognCO2thermosphere.rend(),
								  rmin,
								  -thermosphere_step_r);
    invlognCO2_thermosphere = Linear_interp(lognCO2thermosphere,r_thermosphere);
    lognH_thermosphere_spline = cardinal_cubic_b_spline<Real>(lognHthermosphere.rbegin(),
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
    lognH_exosphere_spline = cardinal_cubic_b_spline<Real>(lognHexosphere.begin(),
							     lognHexosphere.end(),
							     log(rexo),
							     exosphere_step_logr);
    invlognH_exosphere = Linear_interp(lognHexosphere,logr_exosphere);

  }


  Real nCO2(const Real &r) {
    if (r>rexo)
      return 0.0;
    else {
      assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
      return exp(lognCO2_thermosphere_spline(r));
    }
  }

  Real nCO2(const atmo_point pt) {
    return nCO2(pt.r);
  }


  Real nH(const Real &r) {
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

  Real nH(const atmo_point pt) {
    return nH(pt.r);
  }

  Real sH_lya(const Real r) {
    temp->get(r);
    return lyman_alpha_line_center_cross_section_coef/sqrt(temp->T);
  }    
  Real sH_lya(const atmo_point pt) {
    return sH_lya(pt.r);
  }

  Real sCO2_lya(const Real r) {
    return CO2_lyman_alpha_absorption_cross_section;
  }
  Real sCO2_lya(const atmo_point pt) {
    return sCO2_lya(pt.r);
  }

  Real sH_lyb(const Real r) {
    temp->get(r);
    return lyman_beta_line_center_cross_section_coef/sqrt(temp->T);
  }    
  Real sH_lyb(const atmo_point pt) {
    return sH_lyb(pt.r);
  }

  Real sCO2_lyb(const Real r) {
    return CO2_lyman_beta_absorption_cross_section;
  }
  Real sCO2_lyb(const atmo_point pt) {
    return sCO2_lyb(pt.r);
  }


  Real r_from_nH(Real nHtarget) {
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
  

  Real nCO2_exact(const Real &r) {
    if (r>rexo)
      return 0.0;
    else {
      assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
      vector<Real> lognCO2thermosphere_tmp;
      vector<Real> lognHthermosphere_tmp;
      vector<Real> r_thermosphere_tmp;

      vector<Real> nexo(2);
      nexo[0] = log(nCO2exo);
      nexo[1] = log(nHexo);
      integrate( diffeq , nexo , rexo , r , thermosphere_step_r,
		 push_back_quantities( &lognCO2thermosphere_tmp,
				       &lognHthermosphere_tmp,
				       &r_thermosphere_tmp ));
      return exp(lognCO2thermosphere_tmp.back());
    }
  }


  Real nH_exact(const Real &r) {
    if (r>=rexo)
      return exosphere(r);
    else {
      assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
      if (r>=rmindiffusion) {
	vector<Real> lognCO2thermosphere_tmp;
	vector<Real> lognHthermosphere_tmp;
	vector<Real> r_thermosphere_tmp;
	
	vector<Real> nexo(2);
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


  void write_vector(std::ofstream &file, std::string preamble, vector<Real> &data) {
    VectorX write_out = Eigen::Map<VectorX>(data.data(),
					    data.size());

    file << preamble << write_out.transpose() << "\n";
  }
  
  void save(std::string fname) {

    std::ofstream file(fname.c_str());
    if (file.is_open())
      {

	file << "chaffin atmosphere for:\n"
	     << "rexo = "    << rexo << " cm,\n"
	     << "Texo = "    << temp->T_exo << " K,\n"
	     << "nHexo = "   << nHexo << " cm-3,\n"
	     << "nCO2exo = " << nCO2exo << " cm-3,\n\n";

	file << "Thermosphere is defined by solution to Krasnopolsky (2002) differential equation:\n";
	file << "Thermosphere interpolation is log-linear:\n";
	write_vector(file, "r [cm] = ", r_thermosphere);
	write_vector(file, "log(nH) [cm-3] = ", lognHthermosphere);
	write_vector(file, "log(nCO2) [cm-3] = ", lognCO2thermosphere);
	file << std::endl;

	file << "Exosphere is spherically symmetric Chamberlain (nCO2 assumed zero):\n";
	file << "Exosphere interpolation is log-log:\n";
	write_vector(file, "logr [cm] = ", logr_exosphere);
	write_vector(file, "lognH [cm-3] = ", lognHexosphere);

	file.close();
      }
  }
  
};



#endif
