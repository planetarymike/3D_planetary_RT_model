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

  virtual Real nH(const atmo_point &pt) const = 0;
  virtual void nH(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const = 0;
  virtual vector<Real> nH(const vector<atmo_point> &pts) const;
  virtual Real nH(const Real &r) const;

  virtual Real r_from_nH(const Real &nH) const = 0;

  virtual Real nCO2(const atmo_point &pt) const = 0;
  virtual void nCO2(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const = 0;
  virtual vector<Real> nCO2(const vector<atmo_point> &pts) const;
  virtual Real nCO2(const Real &r) const;

  atmosphere(Real rminn, Real rexoo, Real rmaxx);

  Real s_null(const atmo_point &pt) const ;

  //really ought to refactor so cross section info is stored in a
  //totally seperate object
  virtual void sH_lya(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const = 0;
  virtual Real sH_lya(const atmo_point &pt) const = 0;
  virtual Real sH_lya(const Real &r) const;  
  virtual void sCO2_lya(const atmo_voxel &vox, Real &ret_avg, Real &ret_pt) const = 0;
  virtual Real sCO2_lya(const atmo_point &pt) const = 0;
  virtual Real sCO2_lya(const Real &r) const;
};

#endif
