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
  static const int n_int_steps = 1000;
  static const int r_int_scale = 1e8;
  vector<Real> log_r_int;
  vector<Real> n_species_int;
  cardinal_cubic_b_spline<Real> n_species_int_spline;
  vector<Real> n_species_int_spherical;
  cardinal_cubic_b_spline<Real> n_species_int_spline_spherical;
  vector<Real> n_absorber_int;
  Linear_interp n_absorber_int_spline;
  vector<Real> n_absorber_int_spherical;
  Linear_interp n_absorber_int_spline_spherical;
  vector<Real> Tint;
  Linear_interp Tint_spline;
  vector<Real> Tint_spherical;
  Linear_interp Tint_spline_spherical;

  bool spherical;//whether to compute averages in spherical geometry

  bool average_init;//whether integration has been done yet to compute averages
  
  atmosphere_average_1d();

  void setup();

  Real ravg(const Real &r0, const Real &r1,
	    const Real &q0, const Real &q1) const;
  
  Real n_absorber_avg(const Real &r0, const Real &r1) const;
  using atmosphere::n_absorber;
  Real n_absorber(const atmo_point &pt) const; 
  void n_absorber(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;

  Real n_species_avg(const Real &r0, const Real &r1) const;
  using atmosphere::n_species;
  Real n_species(const atmo_point &pt) const;
  void n_species(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
  
  Real Temp_avg(const Real &r0, const Real &r1) const;
  using atmosphere::Temp;
  void Temp(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const;
};

#endif
