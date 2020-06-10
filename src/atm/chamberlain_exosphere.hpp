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
  chamberlain_exosphere(const Real &rexoo, const Real &Texoo, const Real &nHexoo);
  
  Real nH(Real r);
  Real operator()(Real r); //alias for nH

  template <typename T>
  struct nHfinder {
    T *parent;
    Real nHtarget;

    nHfinder(T *parentt, Real &nHtargett)
      : parent(parentt), nHtarget(nHtargett)
    { }
      
    Real operator()(Real r) {
      return parent->nH(r) - nHtarget;
    }
  };
  
  //find r corresponding to a given nH
  Real r(Real &nHtarget);
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
  Linear_interp inv_eff;
  Linear_interp inv_lc;
public:
  Temp_converter(Real rexoo = rexo_typical);

  Real lc_from_T_exact(Real T) const;
  Real eff_from_T_exact(Real T) const;

  Real eff_from_T(Real T) const; 
  Real T_from_eff(Real eff);

  Real lc_from_T(Real T) const; 
  Real T_from_lc(Real eff);
};


#endif
