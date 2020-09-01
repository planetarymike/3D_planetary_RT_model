#ifndef __THERMOSPHERE_DIFFEQ_H
#define __THERMOSPHERE_DIFFEQ_H

#include "temperature.hpp"
#include "diffusion.hpp"

#include <vector>
using std::vector;


// define the differential equation used to solve for the CO2 and H
// number densities in the themosphere
struct thermosphere_diffeq {
  diffusion_coefs diff;
  temperature *temp;
  Real H_escape_flux;
  Real rexo;

  thermosphere_diffeq(temperature &tempp, Real &H_escape_fluxx, Real &rexoo);

  //returns the diffeq derivatives
  void operator()( const vector<Real> &x , vector<Real> &dxdr , const Real &r );

};


#endif
