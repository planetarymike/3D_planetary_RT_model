//tabular_atmosphere.hpp --- use splines to interpolate a table of atmosphere parameters

#include <string>
#include <vector>
using std::vector;

#include "Real.hpp"
#include "constants.hpp"
#include "atmosphere_base.hpp"
#include "interp.hpp"

#include "chamberlain_exosphere.hpp"

struct tabular_atmosphere : virtual public atmosphere {
protected:
  Linear_interp<double> log_n_species_spline;
  Linear_interp<double> inv_log_n_species_spline;
  Linear_interp<double> Temp_spline;
  Linear_interp<double> log_n_absorber_spline;

  bool compute_exosphere;
  chamberlain_exosphere exosphere; 
  
public:
  tabular_atmosphere(double rminn, double rexoo, double rmaxx, bool compute_exospheree = false);

  void load_log_species_density(const vector<double> &alt, const vector<double> &log_n_species);
  void load_log_absorber_density(const vector<double> &alt, const vector<double> &log_n_absorber);
  void load_temperature(const vector<double> &alt, const vector<double> &temp);
  void check_init();
  void init_exosphere();

  double r_from_n_species(const double &n_species_target) const override;

  double n_species(const double &r) const override;

  double Temp(const double &r) const override;

  double n_absorber(const double &r) const override;
};
