
#ifndef __DIFFUSION_H
#define __DIFFUSION_H

#include "Real.hpp"

//diffusion coefficients
struct diffusion_coefs {
  double DH; // diffusion coefficient of hydrogen through CO2
  double KK; // eddy diffusion coefficient

  double DH0;
  double s;
  bool init;

  diffusion_coefs(const double DH00, const double ss);
  
  //sets the above using formulas
  void get(const double &T, const double &Texo, const double &nCO2);  
};

struct diffusion_coefs_H : diffusion_coefs {
  // for hydrogen
  static constexpr double DH0_hydrogen = 8.4e17;// cm^2 s^-1
  static constexpr double s_hydrogen = 0.6;

  diffusion_coefs_H();
};

#endif
