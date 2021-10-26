//wrapper of chamb_diff_1d to force some simple asymmetry

#ifndef __CHAMB_DIFF_1D_ASYM
#define __CHAMB_DIFF_1D_ASYM

#include "Real.hpp"
#include "chamb_diff_1d.hpp"
using std::vector;

struct chamb_diff_1d_asymmetric : public chamb_diff_1d {

  doubReal asymmetry;//multiplicative factor, how much more H there is at
		 //midnight than at noon

  //we use a linear model from noon to midnight to create the anisotropy
  doubReal nslope;
  doubReal n0;
  
  chamb_diff_1d_asymmetric(doubReal n_species_exoo, 
			   doubReal nCO2exoo, //a good number is 10^9 (?)
			   temperature *tempp,
			   species_density_parameters *species_thermospheree);

  chamb_diff_1d_asymmetric(doubReal rminn,
			   doubReal rexoo,
			   doubReal n_species_min,
			   doubReal rmindiffusionn,
			   doubReal n_species_exoo, 
			   doubReal nCO2exoo, //a good number is 10^9 (?)
			   temperature *tempp,
			   species_density_parameters *species_thermospheree);

  void set_asymmetry(const doubReal &a);

  doubReal theta_average_factor(const doubReal &t0, const doubReal &t1) const;

  doubReal nCO2(const atmo_point &pt) const;
  void nCO2(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  //doubReal n_species(const atmo_point &pt) const;
  void n_species_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
};

#endif
