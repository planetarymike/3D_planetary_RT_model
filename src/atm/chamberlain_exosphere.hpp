#ifndef __CHAMBERLAIN_H
#define __CHAMBERLAIN_H

#include "Real.hpp"
#include "constants.hpp"

struct chamberlain_exosphere {
  doubReal rexo;
  doubReal Texo;
  doubReal nexo;
  doubReal m_species;
  doubReal lambdac;
  doubReal effusion_velocity;
  doubReal escape_flux;

  chamberlain_exosphere();
  chamberlain_exosphere(const doubReal rexoo, const doubReal Texoo, const doubReal nexoo, const doubReal m_speciess);
  
  doubReal n(const doubReal r) const;
  doubReal operator()(doubReal r) const; //alias for n

  template <typename T>
  struct nfinder {
    const T *parent;
    const doubReal ntarget;

    nfinder(const T *parentt, const doubReal &ntargett)
      : parent(parentt), ntarget(ntargett)
    { }
      
    doubReal operator()(const doubReal r) const {
      return parent->n(r) - ntarget;
    }
  };
  
  //find r corresponding to a given n
  doubReal r(const doubReal &ntarget) const;
};

#include "constants.hpp"
#include "interp.hpp"
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using boost::math::interpolators::cardinal_cubic_b_spline;

class Temp_converter {
  //class to convert between exobase temperature and escape fraction
protected:
  const doubReal rexo;
  const doubReal m_species;
  static constexpr doubReal  Tmin = 100;
  static constexpr doubReal  Tmax = 1200;
  static const        int    nT = 1101;
  static constexpr doubReal Tstep = (Tmax-Tmin)/(nT-1);


  //forward
  doubReal T_list[nT];
  doubReal lc_list[nT];
  doubReal eff_list[nT];
  cardinal_cubic_b_spline<doubReal> eff_spline;
  cardinal_cubic_b_spline<doubReal> lc_spline;
  Linear_interp<doubReal> inv_eff;
  Linear_interp<doubReal> inv_lc;
public:
  Temp_converter(const doubReal rexoo = rexo_typical, const doubReal m_speciess = mH);

  doubReal lc_from_T_exact(const doubReal T) const;
  doubReal eff_from_T_exact(const doubReal T) const;

  doubReal eff_from_T(const doubReal T) const; 
  doubReal T_from_eff(const doubReal eff) const;

  doubReal lc_from_T(const doubReal T) const; 
  doubReal T_from_lc(const doubReal eff) const;
};


#endif
