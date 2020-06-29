#include "atmosphere_base.hpp"

atmosphere::atmosphere(Real rminn, Real rexoo, Real rmaxx)
  : rmin(rminn), rexo(rexoo), rmax(rmaxx) { }

Real atmosphere::s_null(__attribute__((unused)) const atmo_point &pt) const {
  return 0.0;
}
