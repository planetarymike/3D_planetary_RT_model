#include "chamb_diff_1d.hpp"
using std::vector;

chamb_diff_1d::chamb_diff_1d(double nHexoo, // a good number is 10^5-6
			     double nCO2exoo, //a good number is 10^9 (?)
			     temperature &tempp)
  : chamb_diff_1d(/*          rmin = */rMars + 80e5,
		  /*          rexo = */rexo_typical,
		  /*         nHmin = */10,
		  /* rmindiffusion = */rMars + 80e5,
		  nHexoo,
		  nCO2exoo,
		  tempp,
		  method_nHmin_nCO2exo)   { }

chamb_diff_1d::chamb_diff_1d(double rminn,
			     double rexoo,
			     double rmaxx_or_nHmin,
                             double rmindiffusionn,
                             double nHexoo,   // a good number is 10^5-6
                             double nCO2rmin_or_nCO2exoo, // a good number is 10^9 (?)
                             temperature &tempp,
			     const int method)
  : atmosphere(rminn,rexoo,rmaxx_or_nHmin),
    thermosphere_exosphere(rminn,
			   rexoo,
			   rmaxx_or_nHmin,
			   rmindiffusionn,
			   nHexoo,
			   nCO2rmin_or_nCO2exoo,
			   tempp,
			   method)
{
  setup();
}

void chamb_diff_1d::setup()
{
  atmosphere_average_1d::setup();
}

void chamb_diff_1d::nH(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const {
  //shared with tabular_1d
  atmosphere_average_1d::n_species(vox, ret_avg, ret_pt);
}
void chamb_diff_1d::nCO2(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const {
  //shared with tabular_1d
  atmosphere_average_1d::n_absorber(vox, ret_avg, ret_pt);
}

void chamb_diff_1d::H_Temp(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const {
  //shared with tabular_1d
  if (temp_dependent_sH) {
    atmosphere_average_1d::Temp(vox, ret_avg, ret_pt);
  } else {
    ret_avg = ret_pt = constant_temp_sH;
  }
}
