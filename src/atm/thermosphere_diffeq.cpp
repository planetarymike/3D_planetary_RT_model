#include "thermosphere_diffeq.hpp"
using std::vector;

using std::exp;

thermosphere_diffeq::thermosphere_diffeq(temperature &tempp, Real &H_escape_fluxx, Real &rexoo)
  : temp(&tempp), H_escape_flux(H_escape_fluxx), rexo(rexoo) { }

  //returns the diffeq derivatives
void thermosphere_diffeq::operator()( const vector<Real> &x , vector<Real> &dxdr , const Real &r ) {
  //x[0] = log(nCO2)
  //x[1] = log(nH)

  const Real temp_T      = temp->T(r);
  const Real temp_Tprime = temp->Tprime(r);

  diff.get(temp_T, temp->T_exo, exp(x[0]) );
  
  Real Hninv = G*mMars*mCO2/(kB*temp_T*r*r)+temp_Tprime/temp_T;
  const Real alpha = -0.25;
  Real HHinv = G*mMars*mH/(kB*temp_T*r*r)+(1+alpha)*temp_Tprime/temp_T;
  
  dxdr[0] = -Hninv;
  dxdr[1] = -1./(diff.DH+diff.KK)*( H_escape_flux * ( rexo*rexo /r /r ) / exp(x[1])
				    + (diff.DH*HHinv + diff.KK*Hninv) );
  
}
