//Real_is_float.h -- select type for numeric operations

#ifndef __REAL_FLOAT_H
#define __REAL_FLOAT_H

#include <Eigen/Dense>

typedef float Real;
typedef Eigen::VectorXf VectorX;
typedef Eigen::MatrixXf MatrixX;
#define ABS 1e-3
#define STRICTABS 1e-5

#endif
