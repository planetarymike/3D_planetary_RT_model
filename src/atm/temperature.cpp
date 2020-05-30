#include "temperature.hpp"
#include <cmath>

using std::exp;

krasnopolsky_temperature::krasnopolsky_temperature(Real T_exoo,
						   Real T_tropoo,
						   Real r_tropoo,
						   Real shape_parameterr) {
  T_exo = T_exoo; 
  T_tropo = T_tropoo; 
  r_tropo = r_tropoo; 
  shape_parameter = shape_parameterr;
}

void krasnopolsky_temperature::get(const Real &r) {
  const Real rdiff = r - r_tropo;
  T = T_exo - (T_exo - T_tropo)*exp(-rdiff*rdiff*1e-10/(shape_parameter*T_exo));
  Tprime = ( T_exo - T ) * ( 2*rdiff*1e-10 / (shape_parameter*T_exo) );
}
