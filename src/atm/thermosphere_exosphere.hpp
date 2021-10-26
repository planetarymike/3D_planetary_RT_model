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
  
  doubReal n_species_exo;   // cm-3, species density at exobase
  doubReal nCO2exo; // cm-3, CO2 density at exobase

  doubReal rmindiffusion; // cm, minimum altitude to solve the diffusion equation 
  doubReal n_species_rmindiffusion; // n_species at this altitude
  doubReal nCO2rmindiffusion;
  
  temperature *temp;

  chamberlain_exosphere exosphere;
  species_density_parameters *species_thermosphere;
  chamberlain_exosphere CO2_exosphere;


  //thermosphere interpolation object
  static const int n_thermosphere_steps = 200;
  doubReal thermosphere_step_r;
  vector<doubReal> log_nCO2_thermosphere;
  vector<doubReal> log_n_species_thermosphere;
  vector<doubReal> r_thermosphere;
  Linear_interp<doubReal> log_nCO2_thermosphere_spline;
  Linear_interp<doubReal> invlog_nCO2_thermosphere;
  Linear_interp<doubReal> log_n_species_thermosphere_spline;
  Linear_interp<doubReal> invlog_n_species_thermosphere;
  
  //exosphere interpolation
  static const int n_exosphere_steps = 100;
  doubReal exosphere_step_logr;
  vector<doubReal> log_n_species_exosphere;
  vector<doubReal> log_r_exosphere;
  Linear_interp<doubReal> log_n_species_exosphere_spline;
  Linear_interp<doubReal> invlog_n_species_exosphere;

  doubReal CO2_exo_zero_level = rexo;// + 500e5;
  doubReal exosphere_step_CO2_logr;
  vector<doubReal> log_n_CO2_exosphere;
  vector<doubReal> log_r_CO2_exosphere;
  Linear_interp<doubReal> log_n_CO2_exosphere_spline;
  Linear_interp<doubReal> invlog_n_CO2_exosphere;

  thermosphere_exosphere(doubReal n_species_exoo, // for H, a good number is 10^5-6
			 doubReal nCO2exoo, // a good number is ~10^9
			 temperature *tempp,
			 species_density_parameters *species_thermospheree);
  
  thermosphere_exosphere(doubReal rminn,
			 doubReal rexoo,
			 doubReal rmaxx_or_nspmin,
			 doubReal rmindiffusionn,
			 doubReal n_species_exoo, // for H, a good number is 10^5-6
			 doubReal nCO2rmin_or_nCO2exoo, //a good number is ~10^9
			 temperature *tempp,
			 species_density_parameters *species_thermospheree,
			 const int method = method_nspmin_nCO2exo);

  void setup_nspmin_nCO2exo(doubReal n_species_exoo, // for H, a good number is 10^5-6
			   doubReal nCO2_exoo, // a good number is ~10^9
			   temperature *tempp);
  
  void setup_nspmin_nCO2exo(doubReal rminn,
			   doubReal rexoo,
			   doubReal n_species_min,
			   doubReal rmindiffusionn,
			   doubReal n_species_exoo, // for H, a good number is 10^5-6
			   doubReal nCO2_exoo, // a good number is ~10^9
			   temperature *tempp);

  void setup_rmax_nCO2rmin(doubReal n_species_exoo, // for H, a good number is 10^5-6
			   doubReal nCO2rmin, // a good number is ~10^9
			   temperature *tempp);

  void setup_rmax_nCO2rmin(doubReal rminn,
			   doubReal rexoo,
			   doubReal rmaxx,
			   doubReal rmindiffusionn,
			   doubReal n_species_exoo, // for H, a good number is 10^5-6
			   doubReal nCO2rmin, // a good number is ~10^9
			   temperature *tempp);

  void setup_rmax_nCO2exo(doubReal n_species_exoo, // for H, a good number is 10^5-6
			  doubReal nCO2_exoo, // a good number is ~10^9
			  temperature *tempp);

  void setup_rmax_nCO2exo(doubReal rminn,
			  doubReal rexoo,
			  doubReal n_species_min,
			  doubReal rmindiffusionn,
			  doubReal n_species_exoo, // for H, a good number is 10^5-6
			  doubReal nCO2_exoo, // a good number is ~10^9
			  temperature *tempp);
  
  doubReal nCO2(const doubReal &r) const;
  doubReal n_absorber(const doubReal &r) const override;

  doubReal n_species(const doubReal &r) const override;
  
  doubReal Temp(const doubReal &r) const override;

  doubReal r_from_n_species(const doubReal &n_species) const override;

  doubReal nCO2_exact(const doubReal &r) const;
  doubReal n_species_exact(const doubReal &r) const;

  void write_vector(std::ofstream &file, const std::string &preamble,
		    const vector<doubReal> &data) const;  
  void save(std::string fname) const;  
};

#endif
