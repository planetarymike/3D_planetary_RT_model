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

krasnopolsky_temperature::krasnopolsky_temperature(Real T_exoo,
						   Real T_tropoo,
						   Real r_tropoo,
						   Real shape_parameterr,
						   bool shape_parameter_Texo) {
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
  
  const Real rdiff = r - r_tropo;
  if (rdiff > 0) {
    T_internal      = T_exo - (T_exo - T_tropo)*exp(-rdiff*rdiff*1e-10/shape_parameter);
    Tprime_internal = ( T_exo - T_internal ) * ( 2*rdiff*1e-10 / shape_parameter );
  } else {
    T_internal = T_tropo;
    Tprime_internal = 0;
  }
}
