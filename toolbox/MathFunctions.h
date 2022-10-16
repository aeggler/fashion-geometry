//
// Created by Anna Maria Eggler on 17.10.22.
//

#ifndef EXAMPLE_MATHFUNCTIONS_H
#define EXAMPLE_MATHFUNCTIONS_H

#include <Eigen/Dense>

typedef double Real;

using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
using Matrix4r = Eigen::Matrix<Real, 4, 4, Eigen::DontAlign>;


class MathFunctions{
public:
    static Real cotTheta(const Vector3r &v, const Vector3r &w);
};
#endif //EXAMPLE_MATHFUNCTIONS_H
