
#ifndef __DIFFUSION_H
#define __DIFFUSION_H

#include "Real.hpp"

//diffusion coefficients
struct diffusion_coefs {
  Real DH; // diffusion coefficient of hydrogen through CO2
  Real KK; // eddy diffusion coefficient
  
  //sets the above using formulas
  void get(const Real &T, const Real &Texo, const Real &nCO2);  
};

#endif
