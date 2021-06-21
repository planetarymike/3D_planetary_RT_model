
#ifndef __DIFFUSION_H
#define __DIFFUSION_H

#include "Real.hpp"

//diffusion coefficients
struct diffusion_coefs {
  double DH; // diffusion coefficient of hydrogen through CO2
  double KK; // eddy diffusion coefficient

  const double DH0;
  const double s;

  diffusion_coefs(const double DH00, const double ss);
  
  //sets the above using formulas
  void get(const double &T, const double &Texo, const double &nCO2);  
};

#endif
