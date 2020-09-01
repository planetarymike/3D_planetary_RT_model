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
  Linear_interp log_n_species_spline;
  Linear_interp inv_log_n_species_spline;
  Linear_interp Temp_spline;
  Linear_interp log_n_absorber_spline;

  bool compute_exosphere;
  chamberlain_exosphere exosphere; 
  
public:
  tabular_atmosphere(Real rminn, Real rexoo, Real rmaxx, bool compute_exospheree = false);

  void load_log_species_density(const vector<Real> &alt, const vector<Real> &log_n_species);
  void load_log_absorber_density(const vector<Real> &alt, const vector<Real> &log_n_absorber);
  void load_temperature(const vector<Real> &alt, const vector<Real> &temp);
  void check_init();
  void init_exosphere();

  Real r_from_n_species(const Real &n_species_target) const;

  Real n_species(const Real &r) const;

  Real Temp(const Real &r) const;

  Real n_absorber(const Real &r) const;
};
