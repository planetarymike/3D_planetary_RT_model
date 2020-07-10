#ifndef __ATMOSPHERE_BASE_H_
#define __ATMOSPHERE_BASE_H_

#include "Real.hpp"
#include "atmo_vec.hpp"

struct atmosphere {
  Real rmin;// cm, exobase altitude
  Real rexo;// cm, minimum altitude in model atmosphere
  Real rmax;// cm, max altitude in model atmosphere

  //function to return species density at a given radius (either subsolar or average)
  virtual Real n_species(const Real &r) const = 0;

  //function to return radius (at subsolar point, or average) from species density
  virtual Real r_from_n_species(const Real &n_species) const = 0;
  

  atmosphere(Real rminn, Real rexoo, Real rmaxx);
};

#endif
