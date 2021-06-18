#include "thermosphere_diffeq.hpp"
using std::vector;

using std::exp;

thermosphere_diffeq::thermosphere_diffeq(temperature &tempp, double &H_escape_fluxx, double &rexoo)
  : temp(&tempp), H_escape_flux(H_escape_fluxx), rexo(rexoo)
{

}

  //returns the diffeq derivatives
void thermosphere_diffeq::operator()( const vector<double> &x , vector<double> &dxdr , const double &r ) {
  // this operator returns the derivatives of log(nCO2) and log(nH)
  // for use with the Boost differential equations library

  //  x[0] = log(nCO2)
  //  x[1] = log(nH)

  // get the temperaure at this location
  const double temp_T      = temp->T(r);
  const double temp_Tprime = temp->Tprime(r);

  // set the diffusion coefficients (stored inside diff)
  diff.get(temp_T, temp->T_exo, exp(x[0]) );

  // get the inverse neutral (CO2) scale height
  double Hninv = G*mMars*mCO2/(kB*temp_T*r*r)+temp_Tprime/temp_T;

  // now get the inverse H scale height
  const double alpha = -0.25; // thermal diffusion coefficient
  double HHinv = G*mMars*mH/(kB*temp_T*r*r)+(1+alpha)*temp_Tprime/temp_T;
  
  // now return the derivative of d(log(nCO2))/dr = - 1 / neutral scale height
  dxdr[0] = -Hninv; 

  // and d(log(nH))/dr = - 1/(D + K)*[ phi_H*(r_exo/r)^2/nH
  //                                   + D*(1 / H scale height) + K*(1 / neutral scale height) ]
  dxdr[1] = -1./(diff.DH+diff.KK)*( H_escape_flux * ( rexo*rexo /r /r ) / exp(x[1])
				    + (diff.DH*HHinv + diff.KK*Hninv) );
  
}
