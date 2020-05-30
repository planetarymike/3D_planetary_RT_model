#ifndef __CHAMBERLAIN_H
#define __CHAMBERLAIN_H

#include "Real.hpp"

struct chamberlain_exosphere {
  Real rexo;
  Real Texo;
  Real nHexo;
  Real lambdac;
  Real effusion_velocity;
  Real H_escape_flux;

  chamberlain_exosphere(Real &rexoo, Real &Texoo, Real &nHexoo);
  
  Real nH(Real r);
  Real operator()(Real r); //alias for nH

  template <typename T>
  struct nHfinder {
    T *parent;
    Real nHtarget;

    nHfinder(T *parentt, Real &nHtargett)
      : parent(parentt), nHtarget(nHtargett)
    { }
      
    Real operator()(Real r) {
      return parent->nH(r) - nHtarget;
    }
  };
  
  //find r corresponding to a given nH
  Real r(Real &nHtarget);
};


#endif
