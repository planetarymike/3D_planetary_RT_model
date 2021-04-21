#include "diffusion.hpp"

void diffusion_coefs::get(const double &T, const double &Texo, const double &nCO2) {
  const double DH0 = 8.4e17;// cm^2 s^-1
  const double s = 0.6;
  DH = std::pow(T,s) * DH0/nCO2;
  KK = 1.2e12 * std::sqrt(Texo/nCO2); 
}
