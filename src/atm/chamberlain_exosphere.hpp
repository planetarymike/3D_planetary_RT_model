#ifndef __CHAMBERLAIN_H
#define __CHAMBERLAIN_H

#include "Real.hpp"

struct chamberlain_exosphere {
  double rexo;
  double Texo;
  double nHexo;
  double lambdac;
  double effusion_velocity;
  double H_escape_flux;

  chamberlain_exosphere();
  chamberlain_exosphere(const double rexoo, const double Texoo, const double nHexoo);
  
  double nH(const double r) const;
  double operator()(double r) const; //alias for nH

  template <typename T>
  struct nHfinder {
    const T *parent;
    const double nHtarget;

    nHfinder(const T *parentt, const double &nHtargett)
      : parent(parentt), nHtarget(nHtargett)
    { }
      
    double operator()(const double r) const {
      return parent->nH(r) - nHtarget;
    }
  };
  
  //find r corresponding to a given nH
  double r(const double &nHtarget) const;
};

#include "constants.hpp"
#include "interp.hpp"
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using boost::math::interpolators::cardinal_cubic_b_spline;

class Temp_converter {
  //class to convert between exobase temperature and escape fraction
protected:
  const double rexo;
  static constexpr double  Tmin = 100;
  static constexpr double  Tmax = 1200;
  static const int     nT = 1101;
  static constexpr double Tstep = (Tmax-Tmin)/(nT-1);


  //forward
  double T_list[nT];
  double lc_list[nT];
  double eff_list[nT];
  cardinal_cubic_b_spline<double> eff_spline;
  cardinal_cubic_b_spline<double> lc_spline;
  Linear_interp<double> inv_eff;
  Linear_interp<double> inv_lc;
public:
  Temp_converter(double rexoo = rexo_typical);

  double lc_from_T_exact(const double T) const;
  double eff_from_T_exact(const double T) const;

  double eff_from_T(const double T) const; 
  double T_from_eff(const double eff) const;

  double lc_from_T(const double T) const; 
  double T_from_lc(const double eff) const;
};


#endif
