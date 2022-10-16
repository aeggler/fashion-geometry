
#include "MathFunctions.h"
#include <Eigen/Dense>
const Real eps = static_cast<Real>(1e-6);

typedef double Real;

using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
using Matrix4r = Eigen::Matrix<Real, 4, 4, Eigen::DontAlign>;

Real MathFunctions::cotTheta(const Vector3r &v, const Vector3r &w)
{
    const Real cosTheta = v.dot(w);
    const Real sinTheta = (v.cross(w)).norm();
    return  (cosTheta / sinTheta);
}