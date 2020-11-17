#include "thermosphere_exosphere.hpp"
using std::vector;

#include <boost/numeric/odeint/integrate/integrate.hpp>
using boost::numeric::odeint::integrate;

#include <boost/numeric/odeint/integrate/integrate_const.hpp>
using boost::numeric::odeint::integrate_const;

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
using boost::numeric::odeint::runge_kutta4;

#include <boost/math/quadrature/tanh_sinh.hpp>
using boost::math::quadrature::tanh_sinh;

#include <fstream>

thermosphere_exosphere::thermosphere_exosphere(Real nHexoo, // a good number is 10^5-6
			     Real nCO2exoo, //a good number is 10^9 (?)
			     temperature &tempp)
  : thermosphere_exosphere(/*          rmin = */rMars + 80e5,
			   /*          rexo = */rexo_typical,
			   /*         nHmin = */10,
			   /* rmindiffusion = */rMars + 80e5,
			   nHexoo,
			   nCO2exoo,
			   tempp)   { }

thermosphere_exosphere::thermosphere_exosphere(Real rminn,
					       Real rexoo,
					       Real rmaxx_or_nHmin,
					       Real rmindiffusionn,
					       Real nHexoo, // a good number is 10^5-6
					       Real nCO2rmin_or_nCO2exoo, //a good number is 10^9 (?)
					       temperature &tempp,
					       const int method)
  : atmosphere(rminn, rexoo, -1), // -1 is immediately overwritten
    exosphere(rexoo, tempp.T_exo, nHexoo),// overwritten in setup()
    diffeq(tempp, exosphere.H_escape_flux, rexoo) //overwritten in setup()
{
  if (method == method_nHmin_nCO2exo)
    this->setup_nHmin_nCO2exo(rminn,
			      rexoo,
			      rmaxx_or_nHmin,
			      rmindiffusionn,
			      nHexoo,
			      nCO2rmin_or_nCO2exoo,
			      tempp);
  else if (method == method_rmax_nCO2rmin)
    this->setup_rmax_nCO2rmin(rminn,
			      rexoo,
			      rmaxx_or_nHmin,
			      rmindiffusionn,
			      nHexoo,
			      nCO2rmin_or_nCO2exoo,
			      tempp);
}

void thermosphere_exosphere::setup_nHmin_nCO2exo(Real rminn,
						 Real rexoo,
						 Real nHmin,
						 Real rmindiffusionn,
						 Real nHexoo, // a good number is 10^5-6
						 Real nCO2exoo, //a good number is 10^9 (?)
						 temperature &tempp)
{
  //set the max altitude by finding the density at which the exosphere = nHmin
  Real rmaxx = exosphere.r(nHmin);
  
  setup_rmax_nCO2exo(rminn,
		     rexoo,
		     rmaxx,
		     rmindiffusionn,
		     nHexoo, // a good number is 10^5-6
		     nCO2exoo, //a good number is 10^9 (?)
		     tempp);
}

void thermosphere_exosphere::setup_rmax_nCO2rmin(Real rminn,
						 Real rexoo,
						 Real rmaxx,
						 Real rmindiffusionn,
						 Real nHexoo, // a good number is 10^5-6
						 Real nCO2rmin, //a good number is 2.6e13 (Chaufray 2008)
						 temperature &tempp)
{
  //integrate from rmin upward to rexo to get nCO2exo
  vector<Real> n80km(2);
  n80km[0] = log(nCO2rmin);
  n80km[1] = log(1);//H density is arbitrary, not used later on

  //use a constant stepper for easy interpolation
  Real initial_step = (rexo-rmin)/200.;
  integrate( diffeq , n80km, rmin , rexo , initial_step);

  Real nCO2exoo = exp(n80km[0]);
  
  setup_rmax_nCO2exo(rminn,
		     rexoo,
		     rmaxx,
		     rmindiffusionn,
		     nHexoo, // a good number is 10^5-6
		     nCO2exoo, //a good number is 10^9 (?)
		     tempp);
}

void thermosphere_exosphere::setup_rmax_nCO2exo(Real rminn,
						Real rexoo,
						Real rmaxx,
						Real rmindiffusionn,
						Real nHexoo, // a good number is 10^5-6
						Real nCO2exoo, //a good number is 10^9 (?)
						temperature &tempp)
{
  rmin = rminn;
  rexo = rexoo;
  rmax = rmaxx;
    
  nHexo = nHexoo;
  nCO2exo = nCO2exoo;
  rmindiffusion = rmindiffusionn;

  temp = &tempp;

  exosphere = chamberlain_exosphere(rexoo, temp->T_exo, nHexo);
  diffeq = thermosphere_diffeq(tempp, exosphere.H_escape_flux, rexo);

  //integrate the differential equation to get the species densities
  //in the thermosphere
  vector<Real> nexo(2);
  nexo[0] = log(nCO2exo);
  nexo[1] = log(nHexo);

  //use a constant stepper for easy interpolation
  runge_kutta4< vector<Real> > stepper;
  thermosphere_step_r = -(rexo-rmin)/(n_thermosphere_steps-1.);
  integrate_const( stepper , diffeq ,
		   nexo , rexo , rmin , thermosphere_step_r,
		   push_back_quantities( &lognCO2thermosphere,
					 &lognHthermosphere,
					 &r_thermosphere ) );

#ifndef DNDEBUG
  //check the thermosphere values for any negatives
  for (unsigned int i=0; i < lognCO2thermosphere.size(); i++) {
    assert(exp(lognCO2thermosphere[i]) > 0 && "densities must be positive.");
    assert(exp(lognHthermosphere[i]) > 0 && "densities must be positive.");
    assert(r_thermosphere[i] > 0 && "radii must be positive.");
  }
#endif
  
  //interpolate the densities in the thermosphere
  // lognCO2_thermosphere_spline = cardinal_cubic_b_spline<Real>(lognCO2thermosphere.rbegin(),
  // 							      lognCO2thermosphere.rend(),
  // 							      rmin,
  // 							      -thermosphere_step_r);
  lognCO2_thermosphere_spline = Linear_interp<Real>(r_thermosphere,lognCO2thermosphere);
  invlognCO2_thermosphere = Linear_interp<Real>(lognCO2thermosphere,r_thermosphere);
  // lognH_thermosphere_spline = cardinal_cubic_b_spline<Real>(lognHthermosphere.rbegin(),
  // 							    lognHthermosphere.rend(),
  // 							    rmin,
  // 							    -thermosphere_step_r);
  lognH_thermosphere_spline = Linear_interp<Real>(r_thermosphere,lognHthermosphere);
  invlognH_thermosphere = Linear_interp<Real>(lognHthermosphere,r_thermosphere);

  nHrmindiffusion = nH(rmindiffusion);
  nCO2rmindiffusion = nCO2(rmindiffusion);
    
  //now get the interpolation points in the exosphere
  exosphere_step_logr = (log(rmax) - log(rexo))/(n_exosphere_steps - 1.);
  for (int iexo = 0; iexo < n_exosphere_steps; iexo++) {
    logr_exosphere.push_back( log(rexo) + iexo * exosphere_step_logr );
    assert(logr_exosphere.back() > 0 && "radii must be positive");
    lognHexosphere.push_back( log( exosphere( exp( logr_exosphere[iexo] ) ) ) );
    assert(exp(lognHexosphere.back()) > 0 && "densities must be positive");
  }
  // lognH_exosphere_spline = cardinal_cubic_b_spline<Real>(lognHexosphere.begin(),
  // 							 lognHexosphere.end(),
  // 							 log(rexo),
  // 							 exosphere_step_logr);
  lognH_exosphere_spline = Linear_interp<Real>(logr_exosphere,lognHexosphere);
  invlognH_exosphere = Linear_interp<Real>(lognHexosphere,logr_exosphere);

  init=true;
}


Real thermosphere_exosphere::nCO2(const Real &r) const {
  if (r>rexo)
    return 0.0;
  else {
    assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
    return exp(lognCO2_thermosphere_spline(r));
  }
}
Real thermosphere_exosphere::n_absorber(const Real &r) const {
  return nCO2(r);
}

Real thermosphere_exosphere::nH(const Real &r) const {
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
Real thermosphere_exosphere::n_species(const Real &r) const {
  return nH(r);
}

Real thermosphere_exosphere::Temp(const Real &r) const {
  return temp->T(r);
}

Real thermosphere_exosphere::r_from_nH(const Real &nHtarget) const {
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
Real thermosphere_exosphere::r_from_n_species(const Real &n_species) const {
  return r_from_nH(n_species);
}

Real thermosphere_exosphere::nCO2_exact(const Real &r) const {
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

Real thermosphere_exosphere::nH_exact(const Real &r) const {
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



void thermosphere_exosphere::write_vector(std::ofstream &file, const std::string &preamble,
				 const vector<Real> &data) const {
  VectorX write_out = Eigen::Map<const VectorX>(data.data(),
						data.size());

  file << preamble << write_out.transpose() << "\n";
}
  

void thermosphere_exosphere::save(std::string fname) const {

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
