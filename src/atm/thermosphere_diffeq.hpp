#ifndef __THERMOSPHERE_DIFFEQ_H
#define __THERMOSPHERE_DIFFEQ_H

#include "temperature.hpp"
#include "diffusion.hpp"

#include <vector>
using std::vector;


// define the differential equation used to solve for the CO2 and H
// number densities in the themosphere
struct thermosphere_diffeq {
  diffusion_coefs_H diff;
  temperature *temp;
  double H_escape_flux;
  double rexo;

  thermosphere_diffeq(temperature &tempp, double &H_escape_fluxx, double &rexoo);

  //returns the diffeq derivatives
  void operator()( const vector<double> &x , vector<double> &dxdr , const double &r );

};


#endif
