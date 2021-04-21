#ifndef __TEMPERATURE_H_
#define __TEMPERATURE_H_

#include "Real.hpp"
#include "constants.hpp"

//generic temperature class
struct temperature {
protected:
  double T_internal;
  double Tprime_internal;
  double last_r;

  // a function to define T, Tprime, and last_r
  virtual void get(const double &r) = 0;

public:
  temperature();
  virtual ~temperature() = default;
  
  double T_exo;

  double T(const double &r);
  double Tprime(const double &r);
};

struct krasnopolsky_temperature : virtual public temperature {
  // computes the analytic thermospheric temperature as given by Krasnopolsky (2002)

protected:
  double T_tropo;
  double r_tropo;
  double shape_parameter;

  //implementation of pure virtual from parent
  void get(const double &r) override;

public:
  krasnopolsky_temperature(double T_exoo = 200.0,
			   double T_tropoo = 125.0,
			   double r_tropoo = rMars + 90e5,
			   double shape_parameterr = 11.4,
			   bool shape_parameter_Texo = true);

};

#endif
