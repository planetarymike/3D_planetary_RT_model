//thermosphere_exosphere.hpp --- routine to combine thermosphere + exosphere simulations

#ifndef __THERM_EXO
#define __THERM_EXO

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


//this class combines a diffusive thermosphere and an exosphere with
//an assumed exobase density and temperature profile

struct thermosphere_exosphere : virtual public atmosphere {
  Real nHexo;   // cm-3, H density at exobase
  Real nCO2exo; // cm-3, CO2 density at exobase

  Real rmindiffusion; // cm, minimum altitude to solve the diffusion equation 
  Real nHrmindiffusion; // nH at this altitude
  Real nCO2rmindiffusion;
  
  temperature *temp;

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

  thermosphere_exosphere(Real nHexoo, // a good number is 10^5-6
			 Real nCO2exoo, //a good number is 10^9 (?)
			 temperature &tempp);
  
  thermosphere_exosphere(Real rminn,
			 Real rexoo,
			 Real nHmin,
			 Real rmindiffusionn,
			 Real nHexoo, // a good number is 10^5-6
			 Real nCO2exoo, //a good number is 10^9 (?)
			 temperature &tempp);


  Real nCO2(const Real &r) const;
  Real n_absorber(const Real &r) const;

  Real nH(const Real &r) const;
  Real n_species(const Real &r) const;
  
  Real Temp(const Real &r) const;

  Real r_from_n_species(const Real &n_species) const;
  Real r_from_nH(const Real &nHtarget) const;  

  Real nCO2_exact(const Real &r) const;
  Real nH_exact(const Real &r) const;

  void write_vector(std::ofstream &file, const std::string &preamble,
		    const vector<Real> &data) const;  
  void save(std::string fname) const;  
};

#endif
