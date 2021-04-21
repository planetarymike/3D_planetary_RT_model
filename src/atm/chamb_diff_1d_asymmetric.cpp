#include "chamb_diff_1d_asymmetric.hpp"

chamb_diff_1d_asymmetric::chamb_diff_1d_asymmetric(double nHexoo, // a good number is 10^5-6
						   double nCO2exoo, //a good number is 10^9 (?)
						   temperature &tempp)
  : chamb_diff_1d_asymmetric(/*          rmin = */rMars + 80e5,
			     /*          rexo = */rexo_typical,
			     /*         nHmin = */10,
			     /* rmindiffusion = */rMars + 80e5,
			     nHexoo,
			     nCO2exoo,
			     tempp)
{ }

chamb_diff_1d_asymmetric::chamb_diff_1d_asymmetric(double rminn,
						   double rexoo,
						   double nHmin,
						   double rmindiffusionn,
						   double nHexoo, // a good number is 10^5-6
						   double nCO2exoo, //a good number is 10^9 (?)
						   temperature &tempp)
  : atmosphere(rminn, rexoo, -1),
    chamb_diff_1d(rminn,rexoo,nHmin,rmindiffusionn,nHexoo,nCO2exoo,tempp),
    asymmetry(1.0)
{ }


void chamb_diff_1d_asymmetric::set_asymmetry(const double &a) {
  asymmetry=a;
  
  n0=2.0/(a+1);
  nslope=2.0/pi*(a-1)/(a+1);
}

double chamb_diff_1d_asymmetric::theta_average_factor(const double &t0, const double &t1) const {
  //the t voxels might have values <0 or >pi (so the pt lies on the
  //planet-Sun line) deal with that here
  double tmin = t0 < 0  ? 0 : t0;
  double tmax = t1 > pi ? 0 : t1;

  double cos_tmin = std::cos(tmin);
  double cos_tmax = std::cos(tmax);
  double sin_tmin = std::sin(tmin);
  double sin_tmax = std::sin(tmax);

  return n0 + nslope*(( sin_tmax-sin_tmin + tmin*cos_tmin - tmax*cos_tmax )
		      /(cos_tmin-cos_tmax));
}

double chamb_diff_1d_asymmetric::nCO2(const atmo_point &pt) const {
  //return chamb_diff_1d::thermosphere_exospherez::nCO2(pt.r)*(pt.t*nslope+n0);

  //assume no CO2 asymmetry
  return chamb_diff_1d::thermosphere_exosphere::nCO2(pt.r);
}
void chamb_diff_1d_asymmetric::nCO2(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  //double tint=theta_average_factor(vox.tbounds[0],vox.tbounds[1]);
  
  //ret_avg = tint*chamb_diff_1d::atmosphere_average_1d::n_absorber_avg(vox.rbounds[0],vox.rbounds[1]);
  //ret_pt  =      chamb_diff_1d::thermosphere_exosphere::nCO2(vox.pt.r)*(nslope*vox.pt.t+n0);

  //assume no CO2 asymmetry
  ret_avg = chamb_diff_1d::atmosphere_average_1d::n_absorber_avg(vox.rbounds[0],vox.rbounds[1]);
  ret_pt  = chamb_diff_1d::thermosphere_exosphere::nCO2(vox.pt.r);
}


double chamb_diff_1d_asymmetric::nH(const atmo_point &pt) const {
  return chamb_diff_1d::thermosphere_exosphere::nH(pt.r)*(pt.t*nslope+n0);
}
void chamb_diff_1d_asymmetric::nH(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  double tint=theta_average_factor(vox.tbounds[0],vox.tbounds[1]);

  ret_avg = tint*chamb_diff_1d::atmosphere_average_1d::n_species_avg(vox.rbounds[0],vox.rbounds[1]);
  ret_pt  =      chamb_diff_1d::thermosphere_exosphere::nH(vox.pt.r)*(nslope*vox.pt.t+n0);
}
