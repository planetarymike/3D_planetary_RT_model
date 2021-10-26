#include "chamb_diff_1d_asymmetric.hpp"

chamb_diff_1d_asymmetric::chamb_diff_1d_asymmetric(doubReal n_species_exoo, // a good number is 10^5-6
						   doubReal nCO2exoo, 
						   temperature *tempp,
						   species_density_parameters *species_thermospheree)
  : chamb_diff_1d_asymmetric(/*          rmin = */rMars + 80e5,
			     /*          rexo = */rexo_typical,
			     /* n_species_min = */10,
			     /* rmindiffusion = */rMars + 80e5,
			     n_species_exoo,
			     nCO2exoo,
			     tempp,
			     species_thermospheree)
{ }

chamb_diff_1d_asymmetric::chamb_diff_1d_asymmetric(doubReal rminn,
						   doubReal rexoo,
						   doubReal n_species_min,
						   doubReal rmindiffusionn,
						   doubReal n_species_exoo, 
						   doubReal nCO2exoo, //a good number is 10^9 (?)
						   temperature *tempp,
						   species_density_parameters *species_thermospheree)
  : atmosphere(rminn, rexoo, -1),
    chamb_diff_1d(rminn,rexoo,n_species_min,rmindiffusionn,n_species_exoo,nCO2exoo,tempp, species_thermospheree),
    asymmetry(1.0)
{ }


void chamb_diff_1d_asymmetric::set_asymmetry(const doubReal &a) {
  asymmetry=a;
  
  n0=2.0/(a+1);
  nslope=2.0/pi*(a-1)/(a+1);
}

doubReal chamb_diff_1d_asymmetric::theta_average_factor(const doubReal &t0, const doubReal &t1) const {
  //the t voxels might have values <0 or >pi (so the pt lies on the
  //planet-Sun line) deal with that here
  doubReal tmin = t0 < 0  ? 0 : t0;
  doubReal tmax = t1 > pi ? 0 : t1;

  doubReal cos_tmin = std::cos(tmin);
  doubReal cos_tmax = std::cos(tmax);
  doubReal sin_tmin = std::sin(tmin);
  doubReal sin_tmax = std::sin(tmax);

  return n0 + nslope*(( sin_tmax-sin_tmin + tmin*cos_tmin - tmax*cos_tmax )
		      /(cos_tmin-cos_tmax));
}

doubReal chamb_diff_1d_asymmetric::nCO2(const atmo_point &pt) const {
  //return chamb_diff_1d::thermosphere_exospherez::nCO2(pt.r)*(pt.t*nslope+n0);

  //assume no CO2 asymmetry
  return chamb_diff_1d::thermosphere_exosphere::nCO2(pt.r);
}
void chamb_diff_1d_asymmetric::nCO2(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  //doubReal tint=theta_average_factor(vox.tbounds[0],vox.tbounds[1]);
  
  //ret_avg = tint*chamb_diff_1d::atmosphere_average_1d::n_absorber_avg(vox.rbounds[0],vox.rbounds[1]);
  //ret_pt  =      chamb_diff_1d::thermosphere_exosphere::nCO2(vox.pt.r)*(nslope*vox.pt.t+n0);

  //assume no CO2 asymmetry
  ret_avg = chamb_diff_1d::atmosphere_average_1d::n_absorber_avg(vox.rbounds[0],vox.rbounds[1]);
  ret_pt  = chamb_diff_1d::thermosphere_exosphere::nCO2(vox.pt.r);
}


// doubReal chamb_diff_1d_asymmetric::n_species(const atmo_point &pt) const {
//   return chamb_diff_1d::thermosphere_exosphere::n_species(pt.r)*(pt.t*nslope+n0);
// }
void chamb_diff_1d_asymmetric::n_species_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const {
  doubReal tint=theta_average_factor(vox.tbounds[0],vox.tbounds[1]);

  ret_avg = tint*chamb_diff_1d::atmosphere_average_1d::n_species_avg(vox.rbounds[0],vox.rbounds[1]);
  ret_pt  =      chamb_diff_1d::thermosphere_exosphere::n_species(vox.pt.r)*(nslope*vox.pt.t+n0);
}
