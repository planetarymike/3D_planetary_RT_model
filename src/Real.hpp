//Real.h -- select type for numeric operations

#ifndef __REAL_H
#define __REAL_H

#include "cuda_compatibility.hpp"
#include <Eigen/Dense>

#ifdef RT_FLOAT

typedef float Real;
#define ABS 1e-3
#define STRICTABS 1e-5

#else //default to double

typedef double Real;
#define ABS 1e-6
#define STRICTABS 1e-10

#endif

typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorX;

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

#endif
