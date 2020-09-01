#ifndef __H_XSEC
#define __H_XSEC

#include "Real.hpp"

struct H_cross_sections {
  Real H_lya_xsec_coef;
  Real H_lyb_xsec_coef;

  Real CO2_lya_xsec;
  Real CO2_lyb_xsec;
  
  bool temp_dependent_sH;
  Real constant_temp_sH;
  
  bool no_CO2_absorption;
  
  H_cross_sections();

  Real sH_lya(const Real &T) const;
  Real sCO2_lya(const Real &T) const;

  Real sH_lyb(const Real &T) const;
  Real sCO2_lyb(const Real &T) const;

  void copy_H_options(const H_cross_sections &copy);
};

#endif
