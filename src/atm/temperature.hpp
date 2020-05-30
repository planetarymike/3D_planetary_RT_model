#ifndef __TEMPERATURE_H_
#define __TEMPERATURE_H_

#include "Real.hpp"
#include "constants.hpp"

//generic temperature class
struct temperature {
  Real T_exo;
  Real T;
  Real Tprime;

  // a function to define T and Tprime
  virtual void get(const Real &r) = 0;
};

struct krasnopolsky_temperature : public temperature {
  // computes the analytic thermospheric temperature as given by Krasnopolsky (2002)
  Real T_tropo;
  Real r_tropo;
  Real shape_parameter;

  krasnopolsky_temperature(Real T_exoo = 200.0,
			   Real T_tropoo = 125.0,
			   Real r_tropoo = rMars + 90e5,
			   Real shape_parameterr = 11.4);

  void get(const Real &r);
};

#endif
