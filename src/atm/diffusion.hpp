
#ifndef __DIFFUSION_H
#define __DIFFUSION_H

#include "Real.hpp"

//diffusion coefficients
struct diffusion_coefs {
  doubReal DH; // diffusion coefficient of hydrogen through CO2
  doubReal KK; // eddy diffusion coefficient

  const doubReal DH0;
  const doubReal s;

  diffusion_coefs(const doubReal DH00, const doubReal ss);
  
  //sets the above using formulas
  void get(const doubReal &T, const doubReal &Texo, const doubReal &nCO2);  
};

#endif
