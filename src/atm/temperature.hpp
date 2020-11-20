#ifndef __TEMPERATURE_H_
#define __TEMPERATURE_H_

#include "Real.hpp"
#include "constants.hpp"

//generic temperature class
struct temperature {
protected:
  Real T_internal;
  Real Tprime_internal;
  Real last_r;

  // a function to define T, Tprime, and last_r
  virtual void get(const Real &r) = 0;

public:
  temperature();
  
  Real T_exo;

  Real T(const Real &r);
  Real Tprime(const Real &r);
};

struct krasnopolsky_temperature : public temperature {
  // computes the analytic thermospheric temperature as given by Krasnopolsky (2002)

protected:
  Real T_tropo;
  Real r_tropo;
  Real shape_parameter;

  //implementation of pure virtual from parent
  void get(const Real &r);


public:
  krasnopolsky_temperature(Real T_exoo = 200.0,
			   Real T_tropoo = 125.0,
			   Real r_tropoo = rMars + 90e5,
			   Real shape_parameterr = 11.4,
			   bool shape_parameter_Texo = true);

};

#endif
