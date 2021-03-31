//iph_model_interface.hpp -- routine to call Quemerais fortran IPH model from C++

#include "Real.hpp"
#include <vector>

using std::vector;

vector<Real> quemerais_iph_model(const Real &g_lya, //Lyman alpha g factor at Mars
				 const std::vector<Real> &marspos, //position of Mars in ecliptic coordinates [AU]
				 const vector<Real> &ra, const vector<Real> &dec);

