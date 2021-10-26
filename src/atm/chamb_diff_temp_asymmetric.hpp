//wrapper of chamb_diff_1d to include asymmety in temperature noon-midnight
// and nT^5/2 distribution of H

#ifndef __CHAMB_DIFF_TEMP_ASYM
#define __CHAMB_DIFF_TEMP_ASYM

#include "Real.hpp"
#include "chamb_diff_1d.hpp"
using std::vector;

struct chamb_diff_temp_asymmetric : public atmosphere,
				    public H_cross_sections {
protected:
  //key input parameters
  doubReal navg; //average exobase density
  doubReal T0, T1; //temperature at noon and midnight
  doubReal Tpower; //power in nT^p = const
  doubReal A;  //constant A = nT^p 

  //other input parameters
  doubReal nCO2rmin; //CO2 density at rmin

  //  doubReal n_species_min;// minimum H density to track, sets upper boundary
  doubReal rmindiffusion; // minimum altitude to solve H diffusion equation
		      // (constant below)
  
  //internal 1d atmosphere variables
  static const int n_sza=40;// number of SZA values to calculate
  static constexpr doubReal d_sza=M_PI/(n_sza-1);
  vector<doubReal>            sza_vec;
  krasnopolsky_temperature Temp_sza[n_sza];
  chamb_diff_1d            *atm_sza[n_sza]; //pointers are gross but chamb_diff_1d has an initializer


  //2d integrals over SZA, log radius
  doubReal n_log_r = atmosphere_average_1d::n_int_steps;
  vector<doubReal> log_r;
  Eigen::Matrix<doubReal, Eigen::Dynamic, Eigen::Dynamic> n_species_int;
  Eigen::Matrix<doubReal, Eigen::Dynamic, Eigen::Dynamic> n_absorber_int;
  Eigen::Matrix<doubReal, Eigen::Dynamic, Eigen::Dynamic> Temp_int;

  //setup helper routines
  doubReal T_sza(const doubReal & sza) const; //exobase temperature at this SZA
  doubReal n_species_sza(const doubReal  sza) const; //exobase density at this SZA
  doubReal n_species_sza(const doubReal & sza, doubReal AA) const; //exobase density at SZA with different normalization AA

  doubReal sza_int(const doubReal &f0, const doubReal &f1,
		 const doubReal &t0, const doubReal &t1) const;

  //interpolation object
  Bilinear_interp<doubReal> n_species_int_interp;
  Bilinear_interp<doubReal> n_absorber_int_interp;
  Bilinear_interp<doubReal> Temp_int_interp;

  //computing point values
  void sza_interp(const doubReal &sza, int &i_sza, doubReal &sza_wt) const;

  doubReal r_to_log_r(const doubReal &r) const;
  
  //computing averages
  doubReal avg(const Bilinear_interp<doubReal> &terp,
	     const doubReal &r0, const doubReal &r1,
	     const doubReal &t0, const doubReal &t1) const;

public:

  chamb_diff_temp_asymmetric(species_density_parameters *species_thermospheree,
			     const doubReal n00,
			     const doubReal T00,
			     const doubReal T11);

  chamb_diff_temp_asymmetric(species_density_parameters *species_thermospheree,
			     const doubReal n00,
			     const doubReal T00,
			     const doubReal T11,
			     const doubReal nCO2rminn, //a good number is 2.6e13 (Chaufray2008)
			     const doubReal rexoo,
			     const doubReal rminn,
			     const doubReal rmaxx,
			     const doubReal rmindiffusionn,
			     //extra args for krasnopolsky_temp
			     const doubReal T_tropo = 125.0,
			     const doubReal r_tropo = rMars + 90e5,
			     const doubReal shape_parameter = 11.4,
			     //power for temperature in the expression n*T^p = const.
			     const doubReal Tpowerr = 2.5);

  ~chamb_diff_temp_asymmetric();

  bool spherical = true;

  //  doubReal H_Temp(const atmo_point &pt) const;
  //  void H_Temp(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
  doubReal Temp(const atmo_point &pt) const;
  void Temp_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  doubReal n_absorber(const atmo_point &pt) const; 
  void n_absorber_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
  doubReal nCO2(const atmo_point &pt) const;
  void nCO2_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  doubReal n_species(const atmo_point &pt) const;
  void n_species_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  //stuff required by atmosphere base class
  doubReal n_species(const doubReal &r) const;
  doubReal r_from_n_species(const doubReal &n_species) const;
  doubReal Temp(const doubReal &r) const;
  doubReal n_absorber(const doubReal &r) const;
};
#endif
