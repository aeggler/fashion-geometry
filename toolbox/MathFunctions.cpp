
#include "MathFunctions.h"
#include <Eigen/Dense>
#include <iostream>
const Real eps = static_cast<Real>(1e-6);

typedef double Real;
using namespace std;
using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
using Matrix4r = Eigen::Matrix<Real, 4, 4, Eigen::DontAlign>;

Real MathFunctions::cotTheta(const Vector3r &v, const Vector3r &w)
{
    const Real cosTheta = v.dot(w);
    const Real sinTheta = (v.cross(w)).norm();
    return  (cosTheta / sinTheta);
}

void MathFunctions::Barycentric(const Vector2d& p,const Vector2d& a,const Vector2d& b,const Vector2d& c, Vector3d& ret)
{
    Vector2d v0 = b - a, v1 = c - a, v2 = p - a;
    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    double denom = d00 * d11 - d01 * d01;
    ret(1) = (d11 * d20 - d01 * d21) / denom;
    ret(2) = (d00 * d21 - d01 * d20) / denom;
    ret(0) = 1.0f - ret(1) - ret(2);
}

void MathFunctions::Barycentric3D(const Vector3d& p,const Vector3d& a,const Vector3d& b,const Vector3d& c, Vector3d& baryP)
{
// // http://gamedev.stackexchange.com/a/23745
    Vector3d v0 = b - a, v1 = c - a, v2 = p - a;

    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    double denom = d00 * d11 - d01 * d01;
    baryP(1) = (d11 * d20 - d01 * d21) / denom;
    baryP(2) = (d00 * d21 - d01 * d20) / denom;
    baryP(0) = 1.0f - baryP(1) - baryP(2);
}