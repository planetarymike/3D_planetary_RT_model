#include "tabular_1d.hpp"
using std::vector;

tabular_1d::tabular_1d()
    : atmosphere(-1, -1, -1),
      tabular_atmosphere(-1, -1, -1)
{ }

tabular_1d::tabular_1d(double rminn, double rexoo, double rmaxx,
		       bool compute_exospheree/* = false*/)
  : atmosphere(rminn,rexoo,rmaxx),
    tabular_atmosphere(rminn,rexoo,rmaxx,compute_exospheree)
{
  atmosphere_average_1d::setup();
}

//wrapper functions for the stuff in tabular_atmosphere to ensure
//atmosphere_average_1d is always accurate
void tabular_1d::load_log_species_density(const vector<double> &alt,
                                          const vector<double> &log_n_species) {
  tabular_atmosphere::load_log_species_density(alt,log_n_species);
  atmosphere_average_1d::setup();
  
}
void tabular_1d::load_log_absorber_density(const vector<double> &alt,
					   const vector<double> &log_n_absorber) {
  tabular_atmosphere::load_log_absorber_density(alt,log_n_absorber);
  atmosphere_average_1d::setup();

}
void tabular_1d::load_temperature(const vector<double> &alt,
				  const vector<double> &temp) {
  tabular_atmosphere::load_temperature(alt,temp);
  atmosphere_average_1d::setup();
}


void tabular_1d::Temp_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const {
  // this code is also used in tabular_1d.cpp
  if (temp_dependent_sH) {
    atmosphere_average_1d::Temp_voxel_avg(vox, ret_avg, ret_pt);
  } else {
    ret_avg = ret_pt = constant_temp_sH;
  }
}
