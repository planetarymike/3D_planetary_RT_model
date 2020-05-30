//gauss_legendre_quadrature.h -- determine abcissas and weights for gauss-legendre quadrature

#include "gauss_legendre_quadrature.hpp"
#include <cmath>

using std::vector;

void gauleg(const Real x1, const Real x2, vector<Real> &x, vector<Real> &w)
// Given a lower and upper limit, this returns arrays x and w (of length n),
// containing the abcissas and wieghts of the Gauss-Legendre n-point quadrature 
// formula (see NR pg 183)
{
  Real z1, z, xm, xl, pp, p3, p2, p1;
  int n = x.size(); // the number of abcissas and weights
  int m = (n+1)/2; // number of roots to find since roots are symmetric
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);


  for (int i = 0; i < m; i++) // loop to find roots
    {
      z = std::cos(M_PI * (i + 0.75) / (n + 0.5)); // initial guess for root
      // use Newton's method to find the root at the desired relative precision
      do {
	p1 = 1.0;
	p2 = 0.0;
	for (int j = 0; j < n; j++) // use the Legendre polynomial recurrence
	  {                         //   to evaluate the polynomial at z
	    p3 = p2;                //   uses recurrence relation
	    p2 = p1;                // (j+1)*P_(j+1) = (2j+1)*x*P_j - j*P_(j-1)
	    p1 = ((2*j+1)*z*p2 - j*p3)/(j+1);
	  }
	// p1 is now the legendre polynomial evaluated at z.
	// now get the derivative pp, also by a standard recurrence relation
	pp = n * (z*p1 - p2)/(z*z - 1.0);
	z1 = z;
	z = z1 - p1/pp; // Newton's method
      } while (std::abs(z-z1) > STRICTABS);


      x[i] = xm - xl*z; // scale the root to the interval
      x[n-1-i] = xm + xl*z; // and put in its symmetric counterpart
      w[i] = 2.0*xl / ((1.0 - z*z)*pp*pp); // compute the weight (NR pg. 183)
      w[n-1-i] = w[i]; // and its symmetric counterpart
    }
}
