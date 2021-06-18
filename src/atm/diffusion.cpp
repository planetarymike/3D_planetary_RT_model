#include "diffusion.hpp"

diffusion_coefs::diffusion_coefs(const double DH00, const double ss)
  : DH0(DH00), s(ss) { init=true;}

void diffusion_coefs::get(const double &T, const double &Texo, const double &nCO2) {
  assert(init && "Diffusion coefficient parameters must be initialized!");
  DH = std::pow(T,s) * DH0/nCO2;
  KK = 1.2e12 * std::sqrt(Texo/nCO2); 
}

diffusion_coefs_H::diffusion_coefs_H()
  : diffusion_coefs(DH0_hydrogen, s_hydrogen) { }
