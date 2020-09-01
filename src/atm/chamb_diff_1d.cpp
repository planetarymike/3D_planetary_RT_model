#include "chamb_diff_1d.hpp"
using std::vector;

chamb_diff_1d::chamb_diff_1d(Real nHexoo, // a good number is 10^5-6
			     Real nCO2exoo, //a good number is 10^9 (?)
			     temperature &tempp)
  : chamb_diff_1d(/*          rmin = */rMars + 80e5,
		  /*          rexo = */rexo_typical,
		  /*         nHmin = */10,
		  /* rmindiffusion = */rMars + 120e5,
		  nHexoo,
		  nCO2exoo,
		  tempp)   { }

chamb_diff_1d::chamb_diff_1d(Real rminn,
			     Real rexoo,
			     Real nHmin,
                             Real rmindiffusionn,
                             Real nHexoo,   // a good number is 10^5-6
                             Real nCO2exoo, // a good number is 10^9 (?)
                             temperature &tempp)
  : atmosphere(rminn,rexoo,-1),
    thermosphere_exosphere(rminn,
			   rexoo,
			   nHmin,
			   rmindiffusionn,
			   nHexoo,
			   nCO2exoo,
			   tempp)
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
