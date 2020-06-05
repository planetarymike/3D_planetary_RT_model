//Real.h -- select type for numeric operations

#ifndef __REAL_H
#define __REAL_H

#include "cuda_compatibility.hpp"
#include <Eigen/Dense>

#ifdef RT_FLOAT

typedef float Real;
#define CONEABS 3e-2 //gross, but floating point errors compel this
#define ABS 1e-3
#define STRICTABS 1e-5

#else //default to double

typedef double Real;
#define CONEABS 1e-6
#define ABS 1e-6
#define STRICTABS 1e-10

#endif

typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorX;
typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixX;

#endif
