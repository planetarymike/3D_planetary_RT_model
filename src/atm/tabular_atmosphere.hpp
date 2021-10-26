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
  Linear_interp<doubReal> log_n_species_spline;
  Linear_interp<doubReal> inv_log_n_species_spline;
  Linear_interp<doubReal> Temp_spline;
  Linear_interp<doubReal> log_n_absorber_spline;

  bool compute_exosphere;
  doubReal m_species;
  chamberlain_exosphere exosphere; 
  
public:
  tabular_atmosphere(const doubReal rminn, const doubReal rexoo, const doubReal rmaxx,
		     const bool compute_exospheree = false,
		     const doubReal m_species = mH);

  void load_log_species_density(const vector<doubReal> &alt, const vector<doubReal> &log_n_species);
  void load_log_absorber_density(const vector<doubReal> &alt, const vector<doubReal> &log_n_absorber);
  void load_temperature(const vector<doubReal> &alt, const vector<doubReal> &temp);
  void check_init();
  void init_exosphere();

  doubReal r_from_n_species(const doubReal &n_species_target) const override;

  doubReal n_species(const doubReal &r) const override;

  doubReal Temp(const doubReal &r) const override;

  doubReal n_absorber(const doubReal &r) const override;
};
