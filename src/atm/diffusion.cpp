#include "diffusion.hpp"

void diffusion_coefs::get(const Real &r, const Real &T, const Real &Texo, const Real &nCO2) {
  const Real DH0 = 8.4e17;// cm^2 s^-1
  const Real s = 0.6;
  DH = std::pow(T,s) * DH0/nCO2;
  KK = 1.2e12 * std::sqrt(Texo/nCO2); 
}
