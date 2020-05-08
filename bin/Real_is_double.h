//Real_is_double.h -- select type for numeric operations

#ifndef __REAL_DOUBLE_H
#define __REAL_DOUBLE_H

#include <Eigen/Dense>

typedef double Real;
typedef Eigen::VectorXd VectorX;
typedef Eigen::MatrixXd MatrixX;
#define ABS 1e-6
#define STRICTABS 1e-10
#endif
