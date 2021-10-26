#ifndef __H_XSEC
#define __H_XSEC

#include "Real.hpp"

struct H_cross_sections {
  doubReal H_lya_xsec_coef;
  doubReal H_lyb_xsec_coef;

  doubReal CO2_lya_xsec;
  doubReal CO2_lyb_xsec;
  
  bool temp_dependent_sH;
  doubReal constant_temp_sH;
  
  bool no_CO2_absorption;
  
  H_cross_sections();

  Real sH_lya(const Real &T) const;
  Real sCO2_lya(const Real &T) const;

  Real sH_lyb(const Real &T) const;
  Real sCO2_lyb(const Real &T) const;

  void copy_H_options(const H_cross_sections &copy);
};

#endif
