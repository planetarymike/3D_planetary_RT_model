#include "temperature.hpp"
#include <cmath>

using std::exp;

temperature::temperature() : T_internal(0.0),
			     Tprime_internal(0.0),
			     last_r(0.0)
{}

Real temperature::T(const Real &r) {
  if (r != last_r)
    get(r);
  return T_internal;
}

Real temperature::Tprime(const Real &r) {
  if (r != last_r)
    get(r);
  return Tprime_internal;
}

krasnopolsky_temperature::krasnopolsky_temperature(Real T_exoo/* = 200*/,
						   Real T_tropoo/* = 125*/,
						   Real r_tropoo/* = rMars + 90e5 */,
						   Real shape_parameterr/* = 11.4*/,
						   bool shape_parameter_Texo/* = true*/) {
  T_exo = T_exoo; 
  T_tropo = T_tropoo; 
  r_tropo = r_tropoo;
  if (shape_parameter_Texo)
    shape_parameter = shape_parameterr * T_exo;
  else
    shape_parameter = shape_parameterr * shape_parameterr;
}

void krasnopolsky_temperature::get(const Real &r) {
  last_r = r;
  
  const Real rdiff = (r - r_tropo)*1e-5;
  if (rdiff > 0) {
    T_internal      = T_exo - (T_exo - T_tropo)*exp(-rdiff*rdiff/shape_parameter);
    Tprime_internal = ( T_exo - T_internal ) * ( 2*rdiff / shape_parameter ) * 1e-5;
  } else {
    T_internal = T_tropo;
    Tprime_internal = 0;
  }
}
