#include "diffusion.hpp"

diffusion_coefs::diffusion_coefs(const doubReal DH00, const doubReal ss)
  : DH0(DH00), s(ss) { }

void diffusion_coefs::get(const doubReal &T, const doubReal &Texo, const doubReal &nCO2) {
  DH = std::pow(T,s) * DH0/nCO2;
  KK = 1.2e12 * std::sqrt(Texo/nCO2); 
}
