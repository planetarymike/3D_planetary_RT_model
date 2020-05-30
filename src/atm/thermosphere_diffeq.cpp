#include "thermosphere_diffeq.hpp"
using std::vector;

using std::exp;

thermosphere_diffeq::thermosphere_diffeq(temperature &tempp, Real &H_escape_fluxx, Real &rexoo)
  : temp(&tempp), H_escape_flux(H_escape_fluxx), rexo(rexoo) { }

  //returns the diffeq derivatives
void thermosphere_diffeq::operator()( const vector<Real> &x , vector<Real> &dxdr , const Real &r ) {
  temp->get(r);
  diff.get(r, temp->T, temp->T_exo, exp(x[0]) );
  
  Real Hninv = G*mMars*mCO2/(kB*(temp->T)*r*r)+(temp->Tprime)/(temp->T);
  Real alpha = -0.25;
  Real HHinv = G*mMars*mH/(kB*temp->T*r*r)+(1+alpha)*(temp->Tprime)/(temp->T);
  
  dxdr[0] = -Hninv;
  dxdr[1] = -1./(diff.DH+diff.KK)*( H_escape_flux * ( rexo*rexo /r /r ) / exp(x[1])
				    + (diff.DH*HHinv + diff.KK*Hninv) );
  
}
