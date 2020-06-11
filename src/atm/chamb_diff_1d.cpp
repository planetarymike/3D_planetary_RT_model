#include "chamb_diff_1d.hpp"

#include <boost/numeric/odeint/integrate/integrate.hpp>
using boost::numeric::odeint::integrate;

#include <boost/numeric/odeint/integrate/integrate_const.hpp>
using boost::numeric::odeint::integrate_const;

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
using boost::numeric::odeint::runge_kutta4;

#include <boost/math/quadrature/tanh_sinh.hpp>
using boost::math::quadrature::tanh_sinh;

#include <fstream>

chamb_diff_1d::chamb_diff_1d(Real nHexoo, // a good number is 10^5-6
			     Real nCO2exoo, //a good number is 10^9 (?)
			     temperature &tempp)
  : chamb_diff_1d(/*          rmin = */rMars + 80e5,
		  /*          rexo = */rexo_typical,
		  /*         nHmin = */10,
		  /* rmindiffusion = */rMars + 120e5,
		  nHexoo,
		  nCO2exoo,
		  tempp)   { }

chamb_diff_1d::chamb_diff_1d(Real rminn,
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


  //now we need to integrate the relevant quantites so we can compute averages
  int n_int_steps = 1000;
  Real log_r_int_step = (log(rmax-rMars) - log(rmin-rMars))/(n_int_steps - 1.);
  log_r_int.push_back(log((rmax-rMars)*(1-ABS)));
  nH_int.push_back(0);
  nCO2_int.push_back(0);
  Tint.push_back(0);
  for (int i_int=1; i_int<n_int_steps; i_int++) {
    log_r_int.push_back(log(rmax-rMars)-i_int*log_r_int_step);
    if (exp(log_r_int.back()) < rmin-rMars)
      log_r_int.back() = log((rmin-rMars)*(1+ABS));
    
    Real r0 = exp(log_r_int[i_int-1])+rMars;
    Real r1 = exp(log_r_int[i_int])+rMars;
    Real dr = r0-r1;

    nH_int.push_back( (nH(r1) + nH(r0))/2.0 * dr + nH_int.back() );
    nCO2_int.push_back( (nCO2(r1) + nCO2(r0))/2.0 * dr + nCO2_int.back() );

    temp->get(r0);
    Real T0 = temp->T;
    temp->get(r1);
    Real T1 = temp->T;

    Tint.push_back( (T1 + T0)/2.0 * dr + Tint.back() );
  }
  nH_int_spline = cardinal_cubic_b_spline<Real>(nH_int.rbegin(),
						nH_int.rend(),
						log_r_int.back(),
						log_r_int_step);
  nCO2_int_spline = cardinal_cubic_b_spline<Real>(nCO2_int.rbegin(),
						  nCO2_int.rend(),
						  log_r_int.back(),
						  log_r_int_step);
  Tint_spline = Linear_interp(log_r_int, Tint);
  // Tint_spline = cardinal_cubic_b_spline<Real>(Tint.rbegin(),
  // 					      Tint.rend(),
  // 					      log_r_int.back(),
  // 					      log_r_int_step);
}


Real chamb_diff_1d::nCO2(const Real &r) {
  if (r>rexo)
    return 0.0;
  else {
    assert(r>=rmin && "r must be above the lower boundary of the atmosphere.");
    return exp(lognCO2_thermosphere_spline(r));
  }
}
Real chamb_diff_1d::nCO2(const atmo_point pt) {
  return nCO2(pt.r);
}
void chamb_diff_1d::nCO2(const atmo_voxel vox, Real &ret_avg, Real &ret_pt) {
  ret_avg =  (nCO2_int_spline(log(vox.rbounds[0]-rMars))
	      -nCO2_int_spline(log(vox.rbounds[1]-rMars)))/(vox.rbounds[1]
							    -vox.rbounds[0]);
  ret_pt = nCO2(vox.pt.r);
}


Real chamb_diff_1d::nH(const Real &r) {
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
Real chamb_diff_1d::nH(const atmo_point pt) {
  return nH(pt.r);
}
void chamb_diff_1d::nH(const atmo_voxel vox, Real &ret_avg, Real &ret_pt) {
  ret_avg =  (nH_int_spline(log(vox.rbounds[0]-rMars))
	      -nH_int_spline(log(vox.rbounds[1]-rMars)))/(vox.rbounds[1]
							  -vox.rbounds[0]);
  ret_pt = nH(vox.pt.r);
}



Real chamb_diff_1d::sH_lya(const Real r) {
  temp->get(r);
  return lyman_alpha_line_center_cross_section_coef/sqrt(temp->T);
}    
Real chamb_diff_1d::sH_lya(const atmo_point pt) {
  return sH_lya(pt.r);
}
void chamb_diff_1d::sH_lya(const atmo_voxel vox, Real &ret_avg, Real &ret_pt) {
  Real Tavg = (Tint_spline(log(vox.rbounds[0]-rMars))
	       -Tint_spline(log(vox.rbounds[1]-rMars)))/(vox.rbounds[1]
							 -vox.rbounds[0]);
  ret_avg =  lyman_alpha_line_center_cross_section_coef/sqrt(Tavg);
  ret_pt = sH_lya(vox.pt.r);
}

Real chamb_diff_1d::sCO2_lya(const Real r) {
  return CO2_lyman_alpha_absorption_cross_section;
}
Real chamb_diff_1d::sCO2_lya(const atmo_point pt) {
  return sCO2_lya(pt.r);
}
void chamb_diff_1d::sCO2_lya(const atmo_voxel vox, Real &ret_avg, Real &ret_pt) {
  ret_avg = ret_pt = CO2_lyman_alpha_absorption_cross_section;
}

Real chamb_diff_1d::sH_lyb(const Real r) {
  temp->get(r);
  return lyman_beta_line_center_cross_section_coef/sqrt(temp->T);
}    
Real chamb_diff_1d::sH_lyb(const atmo_point pt) {
  return sH_lyb(pt.r);
}
void chamb_diff_1d::sH_lyb(const atmo_voxel vox, Real &ret_avg, Real &ret_pt) {
  Real Tavg = (Tint_spline(log(vox.rbounds[0]-rMars))
	       -Tint_spline(log(vox.rbounds[1]-rMars)))/(vox.rbounds[1]
							 -vox.rbounds[0]);
  ret_avg = lyman_beta_line_center_cross_section_coef/sqrt(Tavg);
  ret_pt = sH_lyb(vox.pt.r);
}

Real chamb_diff_1d::sCO2_lyb(const Real r) {
  return CO2_lyman_beta_absorption_cross_section;
}
Real chamb_diff_1d::sCO2_lyb(const atmo_point pt) {
  return sCO2_lyb(pt.r);
}
void chamb_diff_1d::sCO2_lyb(const atmo_voxel vox, Real &ret_avg, Real &ret_pt) {
  ret_avg = ret_pt = CO2_lyman_beta_absorption_cross_section;
}




Real chamb_diff_1d::r_from_nH(Real nHtarget) {
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
  




Real chamb_diff_1d::nCO2_exact(const Real &r) {
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

Real chamb_diff_1d::nH_exact(const Real &r) {
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






void chamb_diff_1d::write_vector(std::ofstream &file, std::string preamble, vector<Real> &data) {
  VectorX write_out = Eigen::Map<VectorX>(data.data(),
					  data.size());

  file << preamble << write_out.transpose() << "\n";
}
  



void chamb_diff_1d::save(std::string fname) {

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
