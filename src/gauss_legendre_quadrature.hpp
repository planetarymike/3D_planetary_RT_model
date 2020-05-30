//gauss_legendre_quadrature.h -- determine abcissas and weights for gauss-legendre quadrature

#ifndef __GAUSS_WGTS_H
#define __GAUSS_WGTS_H

#include <vector>
#include "Real.hpp"

void gauleg(const Real x1, const Real x2, std::vector<Real> &x, std::vector<Real> &w);

#endif
