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
  chamb_diff_1d(doubReal n_species_exoo, // for H, a good number is 10^5-6
		doubReal nCO2_exoo, //a good number is ~10^9
		temperature *tempp,
		species_density_parameters *species_thermospheree);

  chamb_diff_1d(doubReal rminn,
		doubReal rexoo,
		doubReal rmaxx_or_nspmin,
		doubReal rmindiffusionn,
		doubReal n_species_exoo, 
		doubReal nCO2exo_or_nCO2rmin,
		temperature *tempp,
		species_density_parameters *species_thermospheree,
		const int method = thermosphere_exosphere::method_nspmin_nCO2exo);

  void setup(); //initialize averages

  // this function is also used in tabular_1d
  using atmosphere_average_1d::n_species_voxel_avg;
  using atmosphere_average_1d::n_absorber_voxel_avg;
  void Temp_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real & ret_pt) const;
};

#endif
