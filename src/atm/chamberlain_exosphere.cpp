#include "chamberlain_exosphere.hpp"

#include <cmath>
using std::exp;
using std::sqrt;
using std::pow;

#include <boost/math/special_functions/gamma.hpp>

chamberlain_exosphere::chamberlain_exosphere() { }
chamberlain_exosphere::chamberlain_exosphere(const doubReal rexoo,
					     const doubReal Texoo,
					     const doubReal nexoo,
					     const doubReal m_speciess)
  : rexo(rexoo), Texo(Texoo), nexo(nexoo), m_species(m_speciess)
{
  lambdac = G*mMars*m_species/(kB*Texo*rexo);//chamberlain lambda @ rexo
  effusion_velocity = 0.5 * sqrt( 2.0*kB*Texo / (m_species*pi) ) * (1.0 + lambdac) * exp(-lambdac); // effusion velocity
  escape_flux = nexo * effusion_velocity;
}

doubReal chamberlain_exosphere::n(const doubReal rr) const {
  // computes hydrogen number density as a function of altitude above
  // the exobase, assuming a chamberlain exosphere w/o satellite particles
  doubReal r = rr;
  
  using boost::math::gamma_p;


  if ((r-rexo)/rexo<ATMEPS) 
    //if we are very close to the exobase due to a rounding error,
    //fix the issue
    r = rexo;
    
  const doubReal lambda = G*mMars*m_species/(kB*Texo*r);//chamberlain lambda
  const doubReal psione = lambda*lambda/(lambda+lambdac);

  //gamma_p = complementary normalized incomplete gamma function
  const doubReal g1 = gamma_p(1.5, lambda);
  const doubReal g1c = gamma_p(1.5, lambdac);
  const doubReal g2 = gamma_p(1.5, lambda - psione);
  
  // calculate the fraction of the exobase density at this altitude:
  doubReal frac = ( 1.0 + g1 - sqrt( 1.0 - lambda*lambda/(lambdac*lambdac) ) * exp(-psione) * ( 1.0 + g2 ) );  
  // this looks different than the formula in C&H b/c the incomplete
  // gamma function in the gsl is normalized; to get gamma(3/2,x) from
  // C&H we have to multiply g1 by Gamma(3/2). This enables us to pull
  // the complete Gamma function outside the parens.

  // normalize to the fraction at the exobase:
  doubReal norm = (1.0 + g1c);
  frac /= norm;

  // multiply by the exobase density and return
  doubReal retval = nexo*frac*exp(lambda-lambdac);
  assert(retval > 0 && "n must be positive");
  
  return retval;
}
doubReal chamberlain_exosphere::operator()(const doubReal r) const {
  return this->n(r);
}

//find r corresponding to a given n
doubReal chamberlain_exosphere::r(const doubReal &ntarget) const {  
  assert(ntarget<nexo && "exosphere n must be less than nexo");


  using boost::math::tools::bracket_and_solve_root;
  using boost::math::tools::eps_tolerance;

    
  doubReal guess = kB*Texo/G/mMars/m_species*log(nexo/ntarget)+rexo;
  doubReal factor = 2;

  const boost::uintmax_t maxit = 40;
  boost::uintmax_t it = maxit;      
  bool is_rising = false;
  int get_digits = 8;
  eps_tolerance<doubReal> tol(get_digits);

  nfinder<chamberlain_exosphere> find(this, ntarget);
  std::pair<doubReal, doubReal> r = bracket_and_solve_root(find,
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



Temp_converter::Temp_converter(const doubReal rexoo/* = rexo_typical*/, const doubReal m_speciess/* = mH*/)
  : rexo(rexoo), m_species(m_speciess)
{
  for (int iT = 0;iT<nT;iT++) {
    T_list[iT]   = Tmin  + iT*Tstep; 
    lc_list[iT]  = lc_from_T_exact(T_list[iT]);
    eff_list[iT] = eff_from_T_exact(T_list[iT]);
  }
  eff_spline = cardinal_cubic_b_spline<doubReal>(eff_list,
					       nT,
					       Tmin,
					       Tstep);
  
  lc_spline = cardinal_cubic_b_spline<doubReal>(lc_list,
					      nT,
					      Tmin,
					      Tstep);

  vector<doubReal> Tvec(T_list,T_list+nT);
  vector<doubReal> lcvec(lc_list,lc_list+nT);
  vector<doubReal> effvec(eff_list,eff_list+nT);
  inv_eff = Linear_interp<doubReal>(effvec,Tvec);
  inv_lc = Linear_interp<doubReal>(lcvec,Tvec);
}

doubReal Temp_converter::lc_from_T_exact(const doubReal T) const {
  return G*mMars*m_species/(kB*T*rexo);
}
doubReal Temp_converter::eff_from_T_exact(const doubReal T) const {
  doubReal lambdac = lc_from_T_exact(T);
  return 0.5 * sqrt( 2.0*kB*T / (m_species*pi) ) * (1.0 + lambdac) * exp(-lambdac);
}

doubReal Temp_converter::eff_from_T(const doubReal T) const {
  return eff_spline(T);
}
doubReal Temp_converter::T_from_eff(const doubReal eff) const {
  return inv_eff(eff);
}

doubReal Temp_converter::lc_from_T(const doubReal T) const {
  return lc_spline(T);
}
doubReal Temp_converter::T_from_lc(const doubReal lc) const {
  return inv_lc(lc);
}
