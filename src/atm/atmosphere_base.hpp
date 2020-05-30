#ifndef __ATMOSPHERE_BASE_H_
#define __ATMOSPHERE_BASE_H_

#include "Real.hpp"
#include "atmo_vec.hpp"
#include <vector>
using std::vector;

struct atmosphere {
  Real rmin;// cm, exobase altitude
  Real rexo;// cm, minimum altitude in model atmosphere
  Real rmax;// cm, max altitude in model atmosphere

  virtual Real nH(const atmo_point pt) = 0;
  virtual vector<Real> nH(const vector<atmo_point> pts);
  virtual Real nH(const Real r);

  virtual Real r_from_nH(const Real nH);

  virtual Real nCO2(const atmo_point pt) = 0;
  virtual vector<Real> nCO2(const vector<atmo_point> pts);
  virtual Real nCO2(const Real r);

  atmosphere(Real rminn, Real rexoo, Real rmaxx);

  Real s_null(const atmo_point pt);

  //really ought to refactor so cross section info is stored in a
  //totally seperate object
  virtual Real sH_lya(const atmo_point pt) = 0;
  virtual Real sH_lya(const Real r);  
  virtual Real sCO2_lya(const atmo_point pt) = 0;
  virtual Real sCO2_lya(const Real r);
};

#endif
