//wrapper of chamb_diff_1d to force some simple asymmetry

#ifndef __CHAMB_DIFF_1D_ASYM
#define __CHAMB_DIFF_1D_ASYM

#include "Real.hpp"
#include "chamb_diff_1d.hpp"
using std::vector;

struct chamb_diff_1d_asymmetric : public chamb_diff_1d {

  double asymmetry;//multiplicative factor, how much more H there is at
		 //midnight than at noon

  //we use a linear model from noon to midnight to create the anisotropy
  double nslope;
  double n0;
  
  chamb_diff_1d_asymmetric(double n_species_exoo, 
			   double nCO2exoo, //a good number is 10^9 (?)
			   temperature *tempp,
			   species_density_parameters *species_thermospheree);

  chamb_diff_1d_asymmetric(double rminn,
			   double rexoo,
			   double n_species_min,
			   double rmindiffusionn,
			   double n_species_exoo, 
			   double nCO2exoo, //a good number is 10^9 (?)
			   temperature *tempp,
			   species_density_parameters *species_thermospheree);

  void set_asymmetry(const double &a);

  double theta_average_factor(const double &t0, const double &t1) const;

  double nCO2(const atmo_point &pt) const;
  void nCO2(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  //double n_species(const atmo_point &pt) const;
  void n_species_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
};

#endif
