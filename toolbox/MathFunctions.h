//
// Created by Anna Maria Eggler on 17.10.22.
//

#ifndef EXAMPLE_MATHFUNCTIONS_H
#define EXAMPLE_MATHFUNCTIONS_H

#include <Eigen/Dense>

typedef double Real;
using namespace Eigen;
using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
using Matrix4r = Eigen::Matrix<Real, 4, 4, Eigen::DontAlign>;


class MathFunctions{
public:
    static Real cotTheta(const Vector3r &v, const Vector3r &w);
    void Barycentric(const Vector2d& p,const Vector2d& a,const Vector2d& b,const Vector2d& c, Vector3d& ret );

};
#endif //EXAMPLE_MATHFUNCTIONS_H
