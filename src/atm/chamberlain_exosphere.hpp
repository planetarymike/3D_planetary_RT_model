#ifndef __CHAMBERLAIN_H
#define __CHAMBERLAIN_H

#include "Real.hpp"

struct chamberlain_exosphere {
  Real rexo;
  Real Texo;
  Real nHexo;
  Real lambdac;
  Real effusion_velocity;
  Real H_escape_flux;

  chamberlain_exosphere();
  chamberlain_exosphere(const Real rexoo, const Real Texoo, const Real nHexoo);
  
  Real nH(const Real r) const;
  Real operator()(Real r) const; //alias for nH

  template <typename T>
  struct nHfinder {
    const T *parent;
    const Real nHtarget;

    nHfinder(const T *parentt, const Real &nHtargett)
      : parent(parentt), nHtarget(nHtargett)
    { }
      
    Real operator()(const Real r) const {
      return parent->nH(r) - nHtarget;
    }
  };
  
  //find r corresponding to a given nH
  Real r(const Real &nHtarget) const;
};

#include "constants.hpp"
#include "interp.hpp"
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using boost::math::interpolators::cardinal_cubic_b_spline;

class Temp_converter {
  //class to convert between exobase temperature and escape fraction
protected:
  const Real rexo;
  static constexpr Real  Tmin = 100;
  static constexpr Real  Tmax = 1200;
  static const int     nT = 1101;
  static constexpr Real Tstep = (Tmax-Tmin)/(nT-1);


  //forward
  Real T_list[nT];
  Real lc_list[nT];
  Real eff_list[nT];
  cardinal_cubic_b_spline<Real> eff_spline;
  cardinal_cubic_b_spline<Real> lc_spline;
  Linear_interp<Real> inv_eff;
  Linear_interp<Real> inv_lc;
public:
  Temp_converter(Real rexoo = rexo_typical);

  Real lc_from_T_exact(const Real T) const;
  Real eff_from_T_exact(const Real T) const;

  Real eff_from_T(const Real T) const; 
  Real T_from_eff(const Real eff) const;

  Real lc_from_T(const Real T) const; 
  Real T_from_lc(const Real eff) const;
};


#endif
