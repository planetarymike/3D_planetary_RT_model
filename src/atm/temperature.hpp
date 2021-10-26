#ifndef __TEMPERATURE_H_
#define __TEMPERATURE_H_

#include "Real.hpp"
#include "constants.hpp"

//generic temperature class
struct temperature {
protected:
  doubReal T_internal;
  doubReal Tprime_internal;
  doubReal last_r;

  // a function to define T, Tprime, and last_r
  virtual void get(const doubReal &r) = 0;

public:
  temperature();
  virtual ~temperature() = default;
  
  doubReal T_exo;

  doubReal T(const doubReal &r);
  doubReal Tprime(const doubReal &r);
};

struct krasnopolsky_temperature : virtual public temperature {
  // computes the analytic thermospheric temperature as given by Krasnopolsky (2002)

protected:
  doubReal T_tropo;
  doubReal r_tropo;
  doubReal shape_parameter;

  //implementation of pure virtual from parent
  void get(const doubReal &r) override;

public:
  krasnopolsky_temperature(doubReal T_exoo = 200.0,
			   doubReal T_tropoo = 125.0,
			   doubReal r_tropoo = rMars + 90e5,
			   doubReal shape_parameterr = 11.4,
			   bool shape_parameter_Texo = true);

};

#endif
