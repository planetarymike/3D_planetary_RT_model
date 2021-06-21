//atmosphere_average_1d.hpp --- routines for computing average densities in a 1D atmosphere
//                       --- inherited by chamb_diff_1d

#ifndef __ATMOSPHERE_1D
#define __ATMOSPHERE_1D

#include "Real.hpp"
#include "constants.hpp"
#include "atmosphere_base.hpp"
#include <vector>
using std::vector;

#include "interp.hpp"

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using boost::math::interpolators::cardinal_cubic_b_spline;

//this class works together with some other atmosphere derived class
//that defines n_species(r), n_absorber(r), and Temp(r). This class
//provides functions to get voxel averages of the quantities provided
//by the other class. This allows multiple methods of defining
//densities to share the same averaging code.

struct atmosphere_average_1d : virtual public atmosphere {
  //integrated quantities to get averages
  static const int n_int_steps = 5000;
  static constexpr double r_int_scale = 1e8;

  //use doubles internally or risk losing precision in computing averages
  vector<double> log_r_int;
  vector<double> n_species_int;
  cardinal_cubic_b_spline<double> n_species_int_spline;
  vector<double> n_species_int_spherical;
  cardinal_cubic_b_spline<double> n_species_int_spline_spherical;
  vector<double> n_absorber_int;
  Linear_interp<double> n_absorber_int_spline;
  vector<double> n_absorber_int_spherical;
  Linear_interp<double> n_absorber_int_spline_spherical;
  vector<double> Tint;
  Linear_interp<double> Tint_spline;
  vector<double> Tint_spherical;
  Linear_interp<double> Tint_spline_spherical;

  bool spherical;//whether to compute averages in spherical geometry

  bool average_init;//whether integration has been done yet to compute averages
  
  atmosphere_average_1d();

  void setup();
  void add_integrated(vector<double> &vec,
		      double q0, double q1,
		      double r0, double r1,
		      bool spherical);

  double ravg(const double &r0, const double &r1,
	    const double &q0, const double &q1) const;
  
  double n_absorber_avg(const double &r0, const double &r1) const;
  using atmosphere::n_absorber;
  //  double n_absorber(const atmo_point &pt) const; 
  void n_absorber_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  double n_species_avg(const double &r0, const double &r1) const;
  using atmosphere::n_species;
  //  double n_species(const atmo_point &pt) const;
  void n_species_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
  
  double Temp_avg(const double &r0, const double &r1) const;
  using atmosphere::Temp;
  void Temp_voxel_avg(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
};

#endif
