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
  
  double nHexo;   // cm-3, H density at exobase
  double nCO2exo; // cm-3, CO2 density at exobase

  double rmindiffusion; // cm, minimum altitude to solve the diffusion equation 
  double nHrmindiffusion; // nH at this altitude
  double nCO2rmindiffusion;
  
  temperature *temp;

  chamberlain_exosphere exosphere;
  thermosphere_diffeq diffeq;

  //thermosphere interpolation object
  static const int n_thermosphere_steps = 100;
  double thermosphere_step_r;
  vector<double> lognCO2thermosphere;
  vector<double> lognHthermosphere;
  vector<double> r_thermosphere;
  //cardinal_cubic_b_spline<double> lognCO2_thermosphere_spline;
  Linear_interp<double> lognCO2_thermosphere_spline;
  Linear_interp<double> invlognCO2_thermosphere;
  //  cardinal_cubic_b_spline<double> lognH_thermosphere_spline;
  Linear_interp<double> lognH_thermosphere_spline;
  Linear_interp<double> invlognH_thermosphere;
  
  //exosphere interpolation
  static const int n_exosphere_steps = 100;
  double exosphere_step_logr;
  vector<double> lognHexosphere;
  vector<double> logr_exosphere;
  //cardinal_cubic_b_spline<double> lognH_exosphere_spline;
  Linear_interp<double> lognH_exosphere_spline;
  Linear_interp<double> invlognH_exosphere;

  thermosphere_exosphere(double nHexoo, // a good number is 10^5-6
			 double nCO2exoo, //a good number is 10^9 (?)
			 temperature &tempp);
  
  thermosphere_exosphere(double rminn,
			 double rexoo,
			 double rmaxx_or_nHmin,
			 double rmindiffusionn,
			 double nHexoo, // a good number is 10^5-6
			 double nCO2rmin_or_nCO2exoo, //a good number is 10^9 (?)
			 temperature &tempp,
			 const int method = method_nHmin_nCO2exo);

  void setup_nHmin_nCO2exo(double nHexoo, // a good number is 10^5-6
			   double nCO2exoo, //a good number is 10^9 (?)
			   temperature &tempp);
  
  void setup_nHmin_nCO2exo(double rminn,
			   double rexoo,
			   double nHmin,
			   double rmindiffusionn,
			   double nHexoo, // a good number is 10^5-6
			   double nCO2exoo, //a good number is 10^9 (?)
			   temperature &tempp);

  void setup_rmax_nCO2rmin(double nHexoo, // a good number is 10^5-6
			   double nCO2rmin, //a good number is 10^9 (?)
			   temperature &tempp);

  void setup_rmax_nCO2rmin(double rminn,
			   double rexoo,
			   double rmaxx,
			   double rmindiffusionn,
			   double nHexoo, // a good number is 10^5-6
			   double nCO2rmin, //a good number is 10^9 (?)
			   temperature &tempp);

  void setup_rmax_nCO2exo(double nHexoo, // a good number is 10^5-6
			  double nCO2exoo, //a good number is 10^9 (?)
			  temperature &tempp);

  void setup_rmax_nCO2exo(double rminn,
			  double rexoo,
			  double nHmin,
			  double rmindiffusionn,
			  double nHexoo, // a good number is 10^5-6
			  double nCO2exoo, //a good number is 10^9 (?)
			  temperature &tempp);
  
  double nCO2(const double &r) const;
  double n_absorber(const double &r) const override;

  double nH(const double &r) const;
  double n_species(const double &r) const override;
  
  double Temp(const double &r) const override;

  double r_from_n_species(const double &n_species) const override;
  double r_from_nH(const double &nHtarget) const;  

  double nCO2_exact(const double &r) const;
  double nH_exact(const double &r) const;

  void write_vector(std::ofstream &file, const std::string &preamble,
		    const vector<double> &data) const;  
  void save(std::string fname) const;  
};

#endif
