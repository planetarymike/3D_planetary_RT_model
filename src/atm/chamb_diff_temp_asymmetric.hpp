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
  double navg; //average exobase density
  double T0, T1; //temperature at noon and midnight
  double Tpower; //power in nT^p = const
  double A;  //constant A = nT^p 

  //other input parameters
  double nCO2rmin; //CO2 density at rmin

  double nHmin;// minimum H density to track, sets upper boundary
  double rmindiffusion; // minimum altitude to solve H diffusion equation
		      // (constant below)
  
  //internal 1d atmosphere variables
  static const int n_sza=40;// number of SZA values to calculate
  static constexpr double d_sza=M_PI/(n_sza-1);
  vector<double>            sza_vec;
  krasnopolsky_temperature Temp_sza[n_sza];
  chamb_diff_1d            *atm_sza[n_sza]; //pointers are gross but chamb_diff_1d has an initializer


  //2d integrals over SZA, log radius
  double n_log_r = atmosphere_average_1d::n_int_steps;
  vector<double> log_r;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> n_species_int;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> n_absorber_int;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Temp_int;

  //setup helper routines
  double T_sza(const double & sza) const; //exobase temperature at this SZA
  double nH_sza(const double  sza) const; //exobase density at this SZA
  double nH_sza(const double & sza, double AA) const; //exobase density at SZA with different normalization AA

  double sza_int(const double &f0, const double &f1,
		 const double &t0, const double &t1) const;

  //interpolation object
  Bilinear_interp<double> n_species_int_interp;
  Bilinear_interp<double> n_absorber_int_interp;
  Bilinear_interp<double> Temp_int_interp;

  //computing point values
  void sza_interp(const double &sza, int &i_sza, double &sza_wt) const;

  double r_to_log_r(const double &r) const;
  
  //computing averages
  double avg(const Bilinear_interp<double> &terp,
	     const double &r0, const double &r1,
	     const double &t0, const double &t1) const;

public:

  chamb_diff_temp_asymmetric(const double n00,
			     const double T00,
			     const double T11);

  chamb_diff_temp_asymmetric(const double n00,
			     const double T00,
			     const double T11,
			     const double nCO2rminn, //a good number is 2.6e13 (Chaufray2008)
			     const double rexoo,
			     const double rminn,
			     const double rmaxx,
			     const double rmindiffusionn,
			     //extra args for krasnopolsky_temp
			     const double T_tropo = 125.0,
			     const double r_tropo = rMars + 90e5,
			     const double shape_parameter = 11.4,
			     //power for temperature in the expression n*T^p = const.
			     const double Tpowerr = 2.5);

  ~chamb_diff_temp_asymmetric();

  bool spherical = true;

  Real H_Temp(const atmo_point &pt) const;
  void H_Temp(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
  Real Temp(const atmo_point &pt) const;
  void Temp(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  Real n_absorber(const atmo_point &pt) const; 
  void n_absorber(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
  Real nCO2(const atmo_point &pt) const;
  void nCO2(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  Real n_species(const atmo_point &pt) const;
  void n_species(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
  Real nH(const atmo_point &pt) const;
  void nH(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  //stuff required by atmosphere base class
  Real n_species(const Real &r) const;
  Real r_from_n_species(const Real &n_species) const;
  Real Temp(const Real &r) const;
  Real n_absorber(const Real &r) const;
};

#endif
