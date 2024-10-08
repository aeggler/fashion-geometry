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
    void Barycentric3D(const Vector3d& p,const Vector3d& a,const Vector3d& b,const Vector3d& c, Vector3d& baryP);
    static void setUpRotationMatrix(double angle,Vector3d& axis, Matrix4d& rotationMatrix);

};
bool raySegmentIntersection(const Vector2d& p, const Vector2d& q,const Vector2d & qs, const Vector2d& ray, double rayMaxLength, Vector2d& intersect);

#endif //EXAMPLE_MATHFUNCTIONS_H
