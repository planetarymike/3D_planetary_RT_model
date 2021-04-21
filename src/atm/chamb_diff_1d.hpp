#ifndef __CHAMB_DIFF_1D
#define __CHAMB_DIFF_1D

#include "hydrogen_cross_sections.hpp"
#include "thermosphere_exosphere.hpp"
#include "atmosphere_average_1d.hpp"
using std::vector;

struct chamb_diff_1d : public thermosphere_exosphere,
		       public atmosphere_average_1d,
		       public H_cross_sections
{
  chamb_diff_1d(double nHexoo, // a good number is 10^5-6
		double nCO2exoo, //a good number is 10^9 (?)
		temperature &tempp);

  chamb_diff_1d(double rminn,
		double rexoo,
		double rmaxx_or_nHmin,
		double rmindiffusionn,
		double nHexoo, // a good number is 10^5-6
		double nCO2exo_or_nCO2rmin, //a good number is 10^9 (?)
		temperature &tempp,
		const int method = thermosphere_exosphere::method_nHmin_nCO2exo);

  void setup(); //initialize averages
  
  //these functions are shared with tabular_1d
  void nH(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const;
  void nCO2(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const;
  void H_Temp(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const;
};

#endif
