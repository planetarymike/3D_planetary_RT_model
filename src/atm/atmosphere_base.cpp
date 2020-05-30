#include "atmosphere_base.hpp"
using std::vector;

vector<Real> atmosphere::nH(const vector<atmo_point> pts) {
  vector<Real> ret;
  ret.resize(pts.size());
  for(unsigned int i=0;i<pts.size();i++)
    ret[i] = nH(pts[i]);
  return ret;
}
Real atmosphere::nH(const Real r) {
  //return subsolar densities if not overridden
  atmo_point p;
  p.rtp(r,0.,0.);
  return nH(p);
}

Real atmosphere::r_from_nH(const Real nH) { return 0; }

vector<Real> atmosphere::nCO2(const vector<atmo_point> pts) {
  vector<Real> ret;
  ret.resize(pts.size());
  for(unsigned int i=0;i<pts.size();i++)
    ret[i] = nCO2(pts[i]);
  return ret;
}
Real atmosphere::nCO2(const Real r) {
  //return subsolar densities if not overridden
  atmo_point p;
  p.rtp(r,0.,0.);
  return nCO2(p);
}

atmosphere::atmosphere(Real rminn, Real rexoo, Real rmaxx)
  : rmin(rminn), rexo(rexoo), rmax(rmaxx) { }

Real atmosphere::s_null(const atmo_point pt) {
  return 0.0;
}

Real atmosphere::sH_lya(const Real r) {
  //return subsolar densities if not overridden
  atmo_point p;
  p.rtp(r,0.,0.);
  return sH_lya(p);
}

Real atmosphere::sCO2_lya(const Real r) {
  //return subsolar densities if not overridden
  atmo_point p;
  p.rtp(r,0.,0.);
  return sCO2_lya(p);
}
