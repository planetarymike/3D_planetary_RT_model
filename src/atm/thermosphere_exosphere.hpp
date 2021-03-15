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
  static const int method_nHmin_nCO2exo = 0;
  static const int method_rmax_nCO2rmin = 1;
  
  Real nHexo;   // cm-3, H density at exobase
  Real nCO2exo; // cm-3, CO2 density at exobase

  Real rmindiffusion; // cm, minimum altitude to solve the diffusion equation 
  Real nHrmindiffusion; // nH at this altitude
  Real nCO2rmindiffusion;
  
  temperature *temp;

  chamberlain_exosphere exosphere;
  thermosphere_diffeq diffeq;

  //thermosphere interpolation object
  static const int n_thermosphere_steps = 100;
  Real thermosphere_step_r;
  vector<Real> lognCO2thermosphere;
  vector<Real> lognHthermosphere;
  vector<Real> r_thermosphere;
  //cardinal_cubic_b_spline<Real> lognCO2_thermosphere_spline;
  Linear_interp<Real> lognCO2_thermosphere_spline;
  Linear_interp<Real> invlognCO2_thermosphere;
  //  cardinal_cubic_b_spline<Real> lognH_thermosphere_spline;
  Linear_interp<Real> lognH_thermosphere_spline;
  Linear_interp<Real> invlognH_thermosphere;
  
  //exosphere interpolation
  static const int n_exosphere_steps = 100;
  Real exosphere_step_logr;
  vector<Real> lognHexosphere;
  vector<Real> logr_exosphere;
  //cardinal_cubic_b_spline<Real> lognH_exosphere_spline;
  Linear_interp<Real> lognH_exosphere_spline;
  Linear_interp<Real> invlognH_exosphere;

  thermosphere_exosphere(Real nHexoo, // a good number is 10^5-6
			 Real nCO2exoo, //a good number is 10^9 (?)
			 temperature &tempp);
  
  thermosphere_exosphere(Real rminn,
			 Real rexoo,
			 Real rmaxx_or_nHmin,
			 Real rmindiffusionn,
			 Real nHexoo, // a good number is 10^5-6
			 Real nCO2rmin_or_nCO2exoo, //a good number is 10^9 (?)
			 temperature &tempp,
			 const int method = method_nHmin_nCO2exo);

  void setup_nHmin_nCO2exo(Real nHexoo, // a good number is 10^5-6
			   Real nCO2exoo, //a good number is 10^9 (?)
			   temperature &tempp);
  
  void setup_nHmin_nCO2exo(Real rminn,
			   Real rexoo,
			   Real nHmin,
			   Real rmindiffusionn,
			   Real nHexoo, // a good number is 10^5-6
			   Real nCO2exoo, //a good number is 10^9 (?)
			   temperature &tempp);

  void setup_rmax_nCO2rmin(Real nHexoo, // a good number is 10^5-6
			   Real nCO2rmin, //a good number is 10^9 (?)
			   temperature &tempp);

  void setup_rmax_nCO2rmin(Real rminn,
			   Real rexoo,
			   Real rmaxx,
			   Real rmindiffusionn,
			   Real nHexoo, // a good number is 10^5-6
			   Real nCO2rmin, //a good number is 10^9 (?)
			   temperature &tempp);

  void setup_rmax_nCO2exo(Real nHexoo, // a good number is 10^5-6
			  Real nCO2exoo, //a good number is 10^9 (?)
			  temperature &tempp);

  void setup_rmax_nCO2exo(Real rminn,
			  Real rexoo,
			  Real nHmin,
			  Real rmindiffusionn,
			  Real nHexoo, // a good number is 10^5-6
			  Real nCO2exoo, //a good number is 10^9 (?)
			  temperature &tempp);
  
  Real nCO2(const Real &r) const;
  Real n_absorber(const Real &r) const override;

  Real nH(const Real &r) const;
  Real n_species(const Real &r) const override;
  
  Real Temp(const Real &r) const override;

  Real r_from_n_species(const Real &n_species) const override;
  Real r_from_nH(const Real &nHtarget) const;  

  Real nCO2_exact(const Real &r) const;
  Real nH_exact(const Real &r) const;

  void write_vector(std::ofstream &file, const std::string &preamble,
		    const vector<Real> &data) const;  
  void save(std::string fname) const;  
};

#endif
