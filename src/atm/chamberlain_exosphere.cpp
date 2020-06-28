#include "chamberlain_exosphere.hpp"

#include <cmath>
using std::exp;
using std::sqrt;
using std::pow;

#include <boost/math/special_functions/gamma.hpp>

chamberlain_exosphere::chamberlain_exosphere() { }
chamberlain_exosphere::chamberlain_exosphere(const Real &rexoo,
					     const Real &Texoo,
					     const Real &nHexoo)
  : rexo(rexoo), Texo(Texoo), nHexo(nHexoo)
{
  lambdac = G*mMars*mH/(kB*Texo*rexo);//chamberlain lambda @ rexo
  effusion_velocity = 0.5 * sqrt( 2.0*kB*Texo / (mH*pi) ) * (1.0 + lambdac) * exp(-lambdac); // effusion velocity
  H_escape_flux = nHexo * effusion_velocity;
}

Real chamberlain_exosphere::nH(const Real rr) const {
  // computes hydrogen number density as a function of altitude above
  // the exobase, assuming a chamberlain exosphere w/o satellite particles
  Real r = rr;
  
  using boost::math::gamma_p;


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
  Real frac = ( 1.0 + g1 - sqrt( 1.0 - lambda*lambda/(lambdac*lambdac) ) * exp(-psione) * ( 1.0 + g2 ) );  
  // this looks different than the formula in C&H b/c the incomplete
  // gamma function in the gsl is normalized; to get gamma(3/2,x) from
  // C&H we have to multiply g1 by Gamma(3/2). This enables us to pull
  // the complete Gamma function outside the parens.

  // normalize to the fraction at the exobase:
  Real norm = (1.0 + g1c);
  frac /= norm;

  // multiply by the exobase density and return
  return nHexo*frac*exp(lambda-lambdac);
}
Real chamberlain_exosphere::operator()(const Real r) const {
  return this->nH(r);
}

//find r corresponding to a given nH
Real chamberlain_exosphere::r(const Real &nHtarget) const {  
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



Temp_converter::Temp_converter(Real rexoo)
  : rexo(rexoo)
{
  for (int iT = 0;iT<nT;iT++) {
    T_list[iT]   = Tmin  + iT*Tstep; 
    lc_list[iT]  = lc_from_T_exact(T_list[iT]);
    eff_list[iT] = eff_from_T_exact(T_list[iT]);
  }
  eff_spline = cardinal_cubic_b_spline<Real>(eff_list,
					     nT,
					     Tmin,
					     Tstep);

  lc_spline = cardinal_cubic_b_spline<Real>(lc_list,
					    nT,
					    Tmin,
					    Tstep);

  vector<Real> Tvec(T_list,T_list+nT);
  vector<Real> lcvec(lc_list,lc_list+nT);
  vector<Real> effvec(eff_list,eff_list+nT);
  inv_eff = Linear_interp(effvec,Tvec);
  inv_lc = Linear_interp(lcvec,Tvec);
}

Real Temp_converter::lc_from_T_exact(const Real T) const {
  return G*mMars*mH/(kB*T*rexo);
}
Real Temp_converter::eff_from_T_exact(const Real T) const {
  Real lambdac = lc_from_T_exact(T);
  return 0.5 * sqrt( 2.0*kB*T / (mH*pi) ) * (1.0 + lambdac) * exp(-lambdac);
}

Real Temp_converter::eff_from_T(const Real T) const {
  return eff_spline(T);
}
Real Temp_converter::T_from_eff(const Real eff) const {
  return inv_eff(eff);
}

Real Temp_converter::lc_from_T(const Real T) const {
  return lc_spline(T);
}
Real Temp_converter::T_from_lc(const Real lc) const {
  return inv_lc(lc);
}
