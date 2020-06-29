#include "chamb_diff_1d_asymmetric.hpp"

chamb_diff_1d_asymmetric::chamb_diff_1d_asymmetric(Real nHexoo, // a good number is 10^5-6
						   Real nCO2exoo, //a good number is 10^9 (?)
						   temperature &tempp)
  : chamb_diff_1d(nHexoo,nCO2exoo,tempp), asymmetry(1.0)
{ }

chamb_diff_1d_asymmetric::chamb_diff_1d_asymmetric(Real rminn,
						   Real rexoo,
						   Real nHmin,
						   Real rmindiffusionn,
						   Real nHexoo, // a good number is 10^5-6
						   Real nCO2exoo, //a good number is 10^9 (?)
						   temperature &tempp)
  : chamb_diff_1d(rminn,rexoo,nHmin,rmindiffusionn,nHexoo,nCO2exoo,tempp),
    asymmetry(1.0)
{ }


void chamb_diff_1d_asymmetric::set_asymmetry(const Real &a) {
  asymmetry=a;
  
  n0=2.0/(a+1);
  nslope=2.0/pi*(a-1)/(a+1);
}

Real chamb_diff_1d_asymmetric::theta_average_factor(const Real &t0, const Real &t1) {
  //the t voxels might have values <0 or >pi (so the pt lies on the
  //planet-Sun line) deal with that here
  Real tmin = t0 < 0  ? 0 : t0;
  Real tmax = t1 > pi ? 0 : t1;

  Real cos_tmin = std::cos(tmin);
  Real cos_tmax = std::cos(tmax);
  Real sin_tmin = std::sin(tmin);
  Real sin_tmax = std::sin(tmax);

  return n0 + nslope*(( sin_tmax-sin_tmin + tmin*cos_tmin - tmax*cos_tmax )
		      /(cos_tmin-cos_tmax));
}

Real chamb_diff_1d_asymmetric::nCO2(const atmo_point &pt) {
  //return chamb_diff_1d::nCO2(pt.r)*(pt.t*nslope+n0);

  //assume no CO2 asymmetry
  return chamb_diff_1d::nCO2(pt.r);
}
void chamb_diff_1d_asymmetric::nCO2(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) {
  //Real tint=theta_average_factor(vox.tbounds[0],vox.tbounds[1]);
  
  //ret_avg = tint*chamb_diff_1d::nCO2avg(vox.rbounds[0],vox.rbounds[1]);
  //ret_pt  =      chamb_diff_1d::nCO2(vox.pt.r)*(nslope*vox.pt.t+n0);

  //assume no CO2 asymmetry
  ret_avg = chamb_diff_1d::nCO2avg(vox.rbounds[0],vox.rbounds[1]);
  ret_pt  = chamb_diff_1d::nCO2(vox.pt.r);
}


Real chamb_diff_1d_asymmetric::nH(const atmo_point &pt) {
  return chamb_diff_1d::nH(pt.r)*(pt.t*nslope+n0);
}
void chamb_diff_1d_asymmetric::nH(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) {
  Real tint=theta_average_factor(vox.tbounds[0],vox.tbounds[1]);

  ret_avg = tint*chamb_diff_1d::nHavg(vox.rbounds[0],vox.rbounds[1]);
  ret_pt  =      chamb_diff_1d::nH(vox.pt.r)*(nslope*vox.pt.t+n0);
}

