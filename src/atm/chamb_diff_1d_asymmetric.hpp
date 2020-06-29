//wrapper of chamb_diff_1d to force some simple asymmetry

#ifndef __CHAMB_DIFF_1D_ASYM
#define __CHAMB_DIFF_1D_ASYM

#include "Real.hpp"
#include "chamb_diff_1d.hpp"
using std::vector;

struct chamb_diff_1d_asymmetric : public chamb_diff_1d {

  Real asymmetry;//multiplicative factor, how much more H there is at
		 //midnight than at noon

  //we use a linear model from noon to midnight to create the anisotropy
  Real nslope;
  Real n0;
  
  chamb_diff_1d_asymmetric(Real nHexoo, // a good number is 10^5-6
			   Real nCO2exoo, //a good number is 10^9 (?)
			   temperature &tempp);

  chamb_diff_1d_asymmetric(Real rminn,
			   Real rexoo,
			   Real nHmin,
			   Real rmindiffusionn,
			   Real nHexoo, // a good number is 10^5-6
			   Real nCO2exoo, //a good number is 10^9 (?)
			   temperature &tempp);

  void set_asymmetry(const Real &a);

  Real theta_average_factor(const Real &t0, const Real &t1);

  Real nCO2(const atmo_point &pt);
  void nCO2(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt);

  Real nH(const atmo_point &pt);
  void nH(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt);

};

#endif
