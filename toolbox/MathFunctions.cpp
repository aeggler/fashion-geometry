
#include "MathFunctions.h"
#include <Eigen/Dense>
#include <iostream>
const Real eps = static_cast<Real>(1e-6);

typedef double Real;
using namespace std;
using namespace Eigen;
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

void MathFunctions::setUpRotationMatrix(double angle, Vector3d& axis, Matrix4d& rotationMatrix)
{
    double u= axis(0);
    double v= axis(1);
    double w= axis(2);
    double L = (u*u + v * v + w * w);
    angle = angle * M_PI / 180.0; //converting to radian value
    double u2 = u * u;
    double v2 = v * v;
    double w2 = w * w;

    rotationMatrix(0,0) = (u2 + (v2 + w2) * cos(angle)) / L;
    rotationMatrix(0,1) = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
    rotationMatrix(0,2) = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
    rotationMatrix(0,3) = 0.0;
    rotationMatrix(1,0) = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
    rotationMatrix(1,1) = (v2 + (u2 + w2) * cos(angle)) / L;
    rotationMatrix(1,2) = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
    rotationMatrix(1,3) = 0.0;
    rotationMatrix(2,0) = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
    rotationMatrix(2,1) = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
    rotationMatrix(2,2) = (w2 + (u2 + v2) * cos(angle)) / L;
    rotationMatrix(2,3) = 0.0;
    rotationMatrix(3,0) = 0.0;
    rotationMatrix(3,1) = 0.0;
    rotationMatrix(3,2) = 0.0;
    rotationMatrix(3,3) = 1.0;
    return;
}
double v2cross(const Vector2d& p, const Vector2d& q){
    return ((p(0)*q(1))-(p(1)*q(0)));
}
//https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect/565282#565282
bool raySegmentIntersection(const Vector2d& p, const Vector2d& q,const Vector2d & qs, const Vector2d& ray, double rayMaxLength, Vector2d& intersect){
    bool doIntersect = false;
    Vector2d r = ray.normalized();
    r *= rayMaxLength;
    Vector2d s = qs-q;

    auto rcross_s= v2cross(r,s);
    if(rcross_s == 0){// they are parallel or colinear
        return doIntersect;
    }
    // compute u, t,  and check for 0 ≤ t ≤ 1 and 0 ≤ u ≤ 1 to see if they intersect
    /*
    t = (q − p) × s / (r × s)

    u = (p − q) × r / (s × r) = (q − p) × r / (r × s)
     */

    auto qp = q-p;
    auto t = v2cross(qp, s) / rcross_s;
    auto u = v2cross(qp,r) / rcross_s;
    if(0 <= t && t <= 1 && 0 <= u && u <= 1 ){
        doIntersect = true;
        intersect = p + t*r;
    }

    return doIntersect;


}