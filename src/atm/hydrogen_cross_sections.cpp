#include "hydrogen_cross_sections.hpp"
#include "constants.hpp"

H_cross_sections::H_cross_sections()
    : H_lya_xsec_coef(lyman_alpha_line_center_cross_section_coef),
      H_lyb_xsec_coef(lyman_beta_line_center_cross_section_coef),
      CO2_lya_xsec(CO2_lyman_alpha_absorption_cross_section),
      CO2_lyb_xsec(CO2_lyman_beta_absorption_cross_section),
      temp_dependent_sH(true),
      constant_temp_sH(-1),
      no_CO2_absorption(false)
{ }

Real H_cross_sections::sH_lya(const Real &T) const {
  Real t_sH = temp_dependent_sH ? T : constant_temp_sH;
  return H_lya_xsec_coef/sqrt(t_sH);
} 
Real H_cross_sections::sCO2_lya(__attribute__((unused)) const Real &T) const {
  if (no_CO2_absorption)
    return 0.0;
  else
    return CO2_lya_xsec;
}

Real H_cross_sections::sH_lyb(const Real &T) const {
  Real t_sH = temp_dependent_sH ? T : constant_temp_sH;
  return H_lyb_xsec_coef/sqrt(t_sH);
}
Real H_cross_sections::sCO2_lyb(__attribute__((unused)) const Real &T) const {
  if (no_CO2_absorption)
    return 0.0;
  else
    return CO2_lyb_xsec;
}

void H_cross_sections::copy_H_options(const H_cross_sections &copy) {
  H_lya_xsec_coef=copy.H_lya_xsec_coef;
  H_lyb_xsec_coef=copy.H_lyb_xsec_coef;

  CO2_lya_xsec=copy.CO2_lya_xsec;

  CO2_lyb_xsec=copy.CO2_lyb_xsec;
  
  temp_dependent_sH=copy.temp_dependent_sH;
  constant_temp_sH=copy.constant_temp_sH;
  
  no_CO2_absorption=copy.no_CO2_absorption;

}
