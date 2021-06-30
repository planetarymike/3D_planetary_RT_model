//Real.h -- select type for numeric operations

#ifndef __REAL_H
#define __REAL_H

#include "cuda_compatibility.hpp"
#include <Eigen/Dense>

#ifdef RT_FLOAT

typedef float Real;
#define REAL(N) (N##f)
#define EPS REAL(1e-3)
#define STRICTEPS REAL(1e-5)
#define CONEEPS REAL(1e-2)

#else //default to double

typedef double Real;
#define REAL(N) (N)
#define EPS 1e-6
#define STRICTEPS 1e-10
#define CONEEPS EPS

#endif

#define ATMEPS 1e-10

typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorX;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;
typedef Eigen::Matrix<float, Eigen::Dynamic, 1> VectorXf;


//if EIGEN_ROWMAJOR is defined then row-major ordering is used for the
//influence matrix instead of Eigen's default Column major ordering

//The GPU prefers row-major (by about 20%) for the influence matrix
//calculation because this is done row-by-row
#ifdef __CUDACC__
#define EIGEN_ROWMAJOR
#endif

#ifdef EIGEN_ROWMAJOR
typedef Eigen::Matrix<Real,
		      Eigen::Dynamic, Eigen::Dynamic,
		      Eigen::RowMajor> MatrixX;
#else
typedef Eigen::Matrix<Real,
		      Eigen::Dynamic, Eigen::Dynamic,
		      Eigen::ColMajor> MatrixX;
#endif



// CUDA needs a constexpr sqrt in order to compile O_1026_tracker
// from here: https://stackoverflow.com/questions/8622256/in-c11-is-sqrt-defined-as-constexpr
#include <limits>   

namespace Detail
{
  Real constexpr sqrtNewtonRaphson(Real x, Real curr, Real prev)
  {
    return curr == prev
      ? curr
      : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
  }
}

/*
 * Constexpr version of the square root
 * Return value:
 *   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
 *   - Otherwise, returns NaN
 */
Real constexpr constexpr_sqrt(Real x)
{
  return x >= 0 && x < std::numeric_limits<Real>::infinity()
		       ? Detail::sqrtNewtonRaphson(x, x, 0)
		       : std::numeric_limits<Real>::quiet_NaN();
}

#endif
