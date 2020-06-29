//csv_atmosphere.hpp --- read data from CSV files into splines for atmospheric parameters

#include "Real.hpp"
#include "constants.hpp"
#include "atmosphere_base.hpp"
#include <vector>
using std::vector;

#include "interp.hpp"

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using boost::math::interpolators::cardinal_cubic_b_spline;

#include <string>

struct csv_atmosphere : atmosphere {
protected:
  vector<Real> Temp_alt;
  vector<Real> Temp;
  Linear_interp T_spline;
  
  vector<Real> n_alt;
  vector<Real> lognCO2;
  Linear_interp lognCO2_spline;
  vector<Real> lognH;
  Linear_interp lognH_spline;
  Linear_interp inv_lognH_spline;

public:
  csv_atmosphere(Real rminn, Real rexoo, Real rmaxx);

  void load_densities(std::string fname);
  void load_temperature(std::string fname);
  void scale_H(Real scale);
  void scale_CO2(Real scale);
  void scale_densities(Real scale);

  Real get_T(const Real &r);

  Real nCO2(const Real &r);
  Real nCO2(const atmo_point pt);
  void nCO2(const atmo_voxel vox, Real &ret_avg, Real &ret_pt);//not implemented

  Real nH(const Real &r);
  Real nH(const atmo_point pt);
  void nH(const atmo_voxel vox, Real &ret_avg, Real &ret_pt);//not implemented

  Real r_from_nH(Real nHtarget);

  Real sH_lya(const atmo_point pt);//not implemented
  void sH_lya(const atmo_voxel vox, Real &ret_avg, Real &ret_pt);//not implemented
  Real sCO2_lya(const atmo_point pt);//not implemented
  void sCO2_lya(const atmo_voxel vox, Real &ret_avg, Real &ret_pt);//not implemented
};
