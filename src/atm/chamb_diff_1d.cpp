#include "chamb_diff_1d.hpp"
using std::vector;

chamb_diff_1d::chamb_diff_1d(double n_species_exoo,
			     double nCO2_exoo, 
			     temperature *tempp,
			     species_density_parameters *species_thermospheree)
  : chamb_diff_1d(/*          rmin = */rMars + 80e5,
		  /*          rexo = */rexo_typical,
		  /* n_species_min = */10,
		  /* rmindiffusion = */rMars + 80e5,
		  n_species_exoo,
		  nCO2_exoo,
		  tempp,
		  species_thermospheree)   { }

chamb_diff_1d::chamb_diff_1d(double rminn,
			     double rexoo,
			     double rmaxx_or_nspmin,
                             double rmindiffusionn,
                             double n_species_exoo,   
                             double nCO2rmin_or_nCO2exoo, 
                             temperature *tempp,
			     species_density_parameters *species_thermospheree,
			     const int method/* = thermosphere_exosphere::method_nspmin_nCO2exo*/)
  : atmosphere(rminn,rexoo,rmaxx_or_nspmin),
    thermosphere_exosphere(rminn,
			   rexoo,
			   rmaxx_or_nspmin,
			   rmindiffusionn,
			   n_species_exoo,
			   nCO2rmin_or_nCO2exoo,
			   tempp,
			   species_thermospheree,
			   method)
{
  setup();
}

void chamb_diff_1d::setup()
{
  atmosphere_average_1d::setup();
}

void chamb_diff_1d::Temp_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const {
  // this code is also used in tabular_1d.cpp
  if (temp_dependent_sH) {
    atmosphere_average_1d::Temp_voxel_avg(vox, ret_avg, ret_pt);
  } else {
    ret_avg = ret_pt = constant_temp_sH;
  }
}
