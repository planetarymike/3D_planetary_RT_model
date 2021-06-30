//thermosphere_exosphere.hpp --- routine to combine thermosphere + exosphere simulations

#ifndef __THERM_EXO
#define __THERM_EXO

#include "Real.hpp"
#include "constants.hpp"
#include "atmosphere_base.hpp"
#include <vector>
using std::vector;

//#include "push_back.hpp"
#include "temperature.hpp"
#include "chamberlain_exosphere.hpp"
#include "species_density_parameters.hpp"
#include "interp.hpp"

// #include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
// using boost::math::interpolators::cardinal_cubic_b_spline;


// this class combines a diffusive thermosphere and an exosphere with
// an assumed exobase density and temperature profile

struct thermosphere_exosphere : virtual public atmosphere {
  static const int method_nspmin_nCO2exo = 0;
  static const int method_rmax_nCO2rmin = 1;
  
  double n_species_exo;   // cm-3, species density at exobase
  double nCO2exo; // cm-3, CO2 density at exobase

  double rmindiffusion; // cm, minimum altitude to solve the diffusion equation 
  double n_species_rmindiffusion; // n_species at this altitude
  double nCO2rmindiffusion;
  
  temperature *temp;

  chamberlain_exosphere exosphere;
  species_density_parameters *species_thermosphere;

  //thermosphere interpolation object
  static const int n_thermosphere_steps = 200;
  double thermosphere_step_r;
  vector<double> log_nCO2_thermosphere;
  vector<double> log_n_species_thermosphere;
  vector<double> r_thermosphere;
  Linear_interp<double> log_nCO2_thermosphere_spline;
  Linear_interp<double> invlog_nCO2_thermosphere;
  Linear_interp<double> log_n_species_thermosphere_spline;
  Linear_interp<double> invlog_n_species_thermosphere;
  
  //exosphere interpolation
  static const int n_exosphere_steps = 100;
  double exosphere_step_logr;
  vector<double> log_n_species_exosphere;
  vector<double> log_r_exosphere;
  Linear_interp<double> log_n_species_exosphere_spline;
  Linear_interp<double> invlog_n_species_exosphere;

  thermosphere_exosphere(double n_species_exoo, // for H, a good number is 10^5-6
			 double nCO2exoo, // a good number is ~10^9
			 temperature *tempp,
			 species_density_parameters *species_thermospheree);
  
  thermosphere_exosphere(double rminn,
			 double rexoo,
			 double rmaxx_or_nspmin,
			 double rmindiffusionn,
			 double n_species_exoo, // for H, a good number is 10^5-6
			 double nCO2rmin_or_nCO2exoo, //a good number is ~10^9
			 temperature *tempp,
			 species_density_parameters *species_thermospheree,
			 const int method = method_nspmin_nCO2exo);

  void setup_nspmin_nCO2exo(double n_species_exoo, // for H, a good number is 10^5-6
			   double nCO2_exoo, // a good number is ~10^9
			   temperature *tempp);
  
  void setup_nspmin_nCO2exo(double rminn,
			   double rexoo,
			   double n_species_min,
			   double rmindiffusionn,
			   double n_species_exoo, // for H, a good number is 10^5-6
			   double nCO2_exoo, // a good number is ~10^9
			   temperature *tempp);

  void setup_rmax_nCO2rmin(double n_species_exoo, // for H, a good number is 10^5-6
			   double nCO2rmin, // a good number is ~10^9
			   temperature *tempp);

  void setup_rmax_nCO2rmin(double rminn,
			   double rexoo,
			   double rmaxx,
			   double rmindiffusionn,
			   double n_species_exoo, // for H, a good number is 10^5-6
			   double nCO2rmin, // a good number is ~10^9
			   temperature *tempp);

  void setup_rmax_nCO2exo(double n_species_exoo, // for H, a good number is 10^5-6
			  double nCO2_exoo, // a good number is ~10^9
			  temperature *tempp);

  void setup_rmax_nCO2exo(double rminn,
			  double rexoo,
			  double n_species_min,
			  double rmindiffusionn,
			  double n_species_exoo, // for H, a good number is 10^5-6
			  double nCO2_exoo, // a good number is ~10^9
			  temperature *tempp);
  
  double nCO2(const double &r) const;
  double n_absorber(const double &r) const override;

  double n_species(const double &r) const override;
  
  double Temp(const double &r) const override;

  double r_from_n_species(const double &n_species) const override;

  double nCO2_exact(const double &r) const;
  double n_species_exact(const double &r) const;

  void write_vector(std::ofstream &file, const std::string &preamble,
		    const vector<double> &data) const;  
  void save(std::string fname) const;  
};

#endif
