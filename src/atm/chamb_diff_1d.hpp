#ifndef __CHAMB_DIFF_1D
#define __CHAMB_DIFF_1D

#include "Real.hpp"
#include "constants.hpp"
#include "atmosphere_base.hpp"
#include <vector>
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
  bool temp_dependent_sH;
  Real constant_temp_sH;

  bool no_CO2_absorption;

  chamberlain_exosphere exosphere;
  thermosphere_diffeq diffeq;

  //thermosphere interpolation object
  static const int n_thermosphere_steps = 40;
  Real thermosphere_step_r;
  vector<Real> lognCO2thermosphere;
  vector<Real> lognHthermosphere;
  vector<Real> r_thermosphere;
  cardinal_cubic_b_spline<Real> lognCO2_thermosphere_spline;
  Linear_interp invlognCO2_thermosphere;
  cardinal_cubic_b_spline<Real> lognH_thermosphere_spline;
  Linear_interp invlognH_thermosphere;
  
  //exosphere interpolation
  static const int n_exosphere_steps = 40;
  Real exosphere_step_logr;
  vector<Real> lognHexosphere;
  vector<Real> logr_exosphere;
  cardinal_cubic_b_spline<Real> lognH_exosphere_spline;
  Linear_interp invlognH_exosphere;

  //integrated quantities to get averages
  static const int n_int_steps = 100;
  static const int r_int_scale = 1e8;
  vector<Real> log_r_int;
  vector<Real> nH_int;
  cardinal_cubic_b_spline<Real> nH_int_spline;
  vector<Real> nH_int_spherical;
  cardinal_cubic_b_spline<Real> nH_int_spline_spherical;
  vector<Real> nCO2_int;
  cardinal_cubic_b_spline<Real> nCO2_int_spline;
  vector<Real> nCO2_int_spherical;
  cardinal_cubic_b_spline<Real> nCO2_int_spline_spherical;
  vector<Real> Tint;
  Linear_interp Tint_spline;
  vector<Real> Tint_spherical;
  Linear_interp Tint_spline_spherical;
  //  cardinal_cubic_b_spline<Real> Tint_spline;

  bool spherical;//whether to compute averages in spherical geometry
		 //or not
  
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


  Real ravg(const Real &r0, const Real &r1,
	    const Real &q0, const Real &q1) const;
  
  Real nCO2(const Real &r) const;
  Real nCO2avg(const Real &r0, const Real &r1) const;
  Real nCO2(const atmo_point &pt) const; 
  void nCO2(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  Real nH(const Real &r) const;
  Real nHavg(const Real &r0, const Real &r1) const;
  Real nH(const atmo_point &pt) const;
  void nH(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
  Real n_species(const Real &r) const;
  
  Real Tavg(const Real &r0, const Real &r1) const;
  void H_Temp(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  Real sH_lya(const Real &r) const;
  Real sH_lya(const atmo_point &pt) const;
  void sH_lya(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  Real sCO2_lya(const Real &r) const;
  Real sCO2_lya(const atmo_point &pt) const;
  void sCO2_lya(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  Real sH_lyb(const Real &r) const;
  Real sH_lyb(const atmo_point &pt) const;
  void sH_lyb(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  Real sCO2_lyb(const Real &r) const;
  Real sCO2_lyb(const atmo_point &pt) const;
  void sCO2_lyb(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  Real r_from_n_species(const Real &n_species) const;
  Real r_from_nH(const Real &nHtarget) const;  

  Real nCO2_exact(const Real &r) const;
  Real nH_exact(const Real &r) const;

  void write_vector(std::ofstream &file, const std::string &preamble,
		    const vector<Real> &data) const;  
  void save(std::string fname) const;  
};

#endif
