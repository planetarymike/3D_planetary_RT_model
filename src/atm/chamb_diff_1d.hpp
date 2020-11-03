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
  chamb_diff_1d(Real nHexoo, // a good number is 10^5-6
		Real nCO2exoo, //a good number is 10^9 (?)
		temperature &tempp);

  chamb_diff_1d(Real rminn,
		Real rexoo,
		Real nHmin,
		Real rmindiffusionn,
		Real nHexoo, // a good number is 10^5-6
		Real nCO2exoo, //a good number is 10^9 (?)
		temperature &tempp);

  void setup(Real nHexoo, // a good number is 10^5-6
	     Real nCO2exoo, //a good number is 10^9 (?)
	     temperature &tempp);
  
  void setup(Real rminn,
	     Real rexoo,
	     Real nHmin,
	     Real rmindiffusionn,
	     Real nHexoo, // a good number is 10^5-6
	     Real nCO2exoo, //a good number is 10^9 (?)
	     temperature &tempp);
  
  //these functions are shared with tabular_1d
  void nH(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const;
  void nCO2(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const;
  void H_Temp(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const;
};

#endif
