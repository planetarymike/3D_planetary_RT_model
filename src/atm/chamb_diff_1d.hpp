#ifndef __CHAMB_DIFF_1D
#define __CHAMB_DIFF_1D

#include "Real.hpp"
#include "constants.hpp"
#include "atmosphere_base.hpp"
using std::vector;

#include "push_back.hpp"
#include "temperature.hpp"
#include "chamberlain_exosphere.hpp"
#include "thermosphere_diffeq.hpp"
#include "interp.hpp"

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using boost::math::interpolators::cardinal_cubic_b_spline;

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

  //integrated quantities to get averages
  vector<Real> log_r_int;
  vector<Real> nH_int;
  cardinal_cubic_b_spline<Real> nH_int_spline;
  vector<Real> nCO2_int;
  cardinal_cubic_b_spline<Real> nCO2_int_spline;
  vector<Real> Tint;
  Linear_interp Tint_spline;
  //  cardinal_cubic_b_spline<Real> Tint_spline;

  chamb_diff_1d(Real nHexoo, // a good number is 10^5-6
		Real nCO2exoo, //a good number is 10^9 (?)
		temperature &tempp);

  chamb_diff_1d(Real rminn,
		Real rexoo,
		Real nHmin,
		Real rmindiffusionn,
		Real nHexoo, // a good number is 10^5-6
		Real nCO2exoo, //a good number is 10^9 (?)
		temperature &tempp);


  Real nCO2(const Real &r);
  Real nCO2(const atmo_point pt);
  void nCO2(const atmo_voxel vox, Real &ret_avg, Real &ret_pt);

  Real nH(const Real &r);
  Real nH(const atmo_point pt);
  void nH(const atmo_voxel vox, Real &ret_avg, Real &ret_pt);

  Real sH_lya(const Real r);
  Real sH_lya(const atmo_point pt);
  void sH_lya(const atmo_voxel vox, Real &ret_avg, Real &ret_pt);

  Real sCO2_lya(const Real r);
  Real sCO2_lya(const atmo_point pt);
  void sCO2_lya(const atmo_voxel vox, Real &ret_avg, Real &ret_pt);

  Real sH_lyb(const Real r);
  Real sH_lyb(const atmo_point pt);
  void sH_lyb(const atmo_voxel vox, Real &ret_avg, Real &ret_pt);

  Real sCO2_lyb(const Real r);
  Real sCO2_lyb(const atmo_point pt);
  void sCO2_lyb(const atmo_voxel vox, Real &ret_avg, Real &ret_pt);

  Real r_from_nH(Real nHtarget);  

  Real nCO2_exact(const Real &r);
  Real nH_exact(const Real &r);

  void write_vector(std::ofstream &file, std::string preamble, vector<Real> &data);  
  void save(std::string fname);  
};

#endif
