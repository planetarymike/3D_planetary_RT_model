#ifndef __TABULAR_1D
#define __TABULAR_1D

#include "hydrogen_cross_sections.hpp"
#include "tabular_atmosphere.hpp"
#include "atmosphere_average_1d.hpp"
using std::vector;

struct tabular_1d: public tabular_atmosphere,
		   public atmosphere_average_1d,
		   public H_cross_sections
{
  tabular_1d();
  tabular_1d(double rminn, double rexoo, double rmaxx, bool compute_exospheree = false);

  //wrapper functions for the stuff in tabular_atmosphere to ensure
  //atmosphere_average_1d is always accurate
  void load_log_species_density(const vector<double> &alt, const vector<double> &log_n_species);
  void load_log_absorber_density(const vector<double> &alt, const vector<double> &log_n_absorber);
  void load_temperature(const vector<double> &alt, const vector<double> &temp);

  using tabular_atmosphere::n_species;
  using tabular_atmosphere::n_absorber;
  using tabular_atmosphere::Temp;

  //this function is also used in chamb_diff_1d
  void Temp_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const;
};

#endif
