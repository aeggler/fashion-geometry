
#include "PositionBasedDynamics.h"
#include <Eigen/Dense>
#include "MathFunctions.h"
const Real eps = static_cast<Real>(1e-6);

typedef double Real;

using Vector2r = Eigen::Matrix<Real, 2, 1, Eigen::DontAlign>;

using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
using Matrix2r = Eigen::Matrix<Real, 2, 2, Eigen::DontAlign>;
using Matrix3r = Eigen::Matrix<Real, 3, 3, Eigen::DontAlign>;
using Matrix4r = Eigen::Matrix<Real, 4, 4, Eigen::DontAlign>;

// ------------------------------------------------------------------------------------

bool PositionBasedDynamics::init_IsometricBendingConstraint(
        const Vector3r &p0,
        const Vector3r &p1,
        const Vector3r &p2,
        const Vector3r &p3,
        Matrix4r &Q)
{
    // Compute matrix Q for quadratic bending
    const Vector3r *x[4] = { &p2, &p3, &p0, &p1 };

    const Vector3r e0 = *x[1] - *x[0];
    const Vector3r e1 = *x[2] - *x[0];
    const Vector3r e2 = *x[3] - *x[0];
    const Vector3r e3 = *x[2] - *x[1];
    const Vector3r e4 = *x[3] - *x[1];

    const Real c01 = MathFunctions::cotTheta(e0, e1);
    const Real c02 = MathFunctions::cotTheta(e0, e2);
    const Real c03 = MathFunctions::cotTheta(-e0, e3);
    const Real c04 = MathFunctions::cotTheta(-e0, e4);

    const Real A0 = static_cast<Real>(0.5) * (e0.cross(e1)).norm();
    const Real A1 = static_cast<Real>(0.5) * (e0.cross(e2)).norm();

    const Real coef = -3.f / (2.f*(A0 + A1));
    const Real K[4] = { c03 + c04, c01 + c02, -c01 - c03, -c02 - c04 };
    const Real K2[4] = { coef*K[0], coef*K[1], coef*K[2], coef*K[3] };

    for (unsigned char j = 0; j < 4; j++)
    {
        for (unsigned char k = 0; k < j; k++)
        {
            Q(j, k) = Q(k, j) = K[j] * K2[k];
        }
        Q(j, j) = K[j] * K2[j];
    }

    return true;
}

// ----------------------------------------------------------------------------------------------
// returns the position correction delta p.
bool PositionBasedDynamics::solve_IsometricBendingConstraint(
        const Vector3r &p0, Real invMass0,
        const Vector3r &p1, Real invMass1,
        const Vector3r &p2, Real invMass2,
        const Vector3r &p3, Real invMass3,
        const Matrix4r &Q,
        const Real stiffness,
        Vector3r &corr0, Vector3r &corr1, Vector3r &corr2, Vector3r &corr3)
{
    const Vector3r *x[4] = { &p2, &p3, &p0, &p1 };
    Real invMass[4] = { invMass2, invMass3, invMass0, invMass1 };

    Real energy = 0.0;
    for (unsigned char k = 0; k < 4; k++)
        for (unsigned char j = 0; j < 4; j++)
            energy += Q(j, k)*(x[k]->dot(*x[j]));
    energy *= 0.5;

    Vector3r gradC[4];
    gradC[0].setZero();
    gradC[1].setZero();
    gradC[2].setZero();
    gradC[3].setZero();
    for (unsigned char k = 0; k < 4; k++)
        for (unsigned char j = 0; j < 4; j++)
            gradC[j] += Q(j,k) * *x[k];


    Real sum_normGradC = 0.0;
    for (unsigned int j = 0; j < 4; j++)
    {
        // compute sum of squared gradient norms
        if (invMass[j] != 0.0)
            sum_normGradC += invMass[j] * gradC[j].squaredNorm();
    }

    // exit early if required
    if (fabs(sum_normGradC) > eps)
    {
        // compute impulse-based scaling factor
        const Real s = energy / sum_normGradC;

        corr0 = -stiffness * (s*invMass[2]) * gradC[2];
        corr1 = -stiffness * (s*invMass[3]) * gradC[3];
        corr2 = -stiffness * (s*invMass[0]) * gradC[0];
        corr3 = -stiffness * (s*invMass[1]) * gradC[1];

        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------------------------
bool PositionBasedDynamics::solve_DistanceConstraint(
        const Vector3r &p0, Real invMass0,
        const Vector3r &p1, Real invMass1,
        const Real restLength,
        const Real stiffness,
        Vector3r &corr0, Vector3r &corr1)
{
    Real wSum = invMass0 + invMass1;
    if (wSum == 0.0)
        return false;

    Vector3r n = p1 - p0;
    Real d = n.norm();
    n.normalize();

    Vector3r corr;
    corr = stiffness * n * (d - restLength) / wSum;

    corr0 =  invMass0 * corr;
    corr1 = -invMass1 * corr;
    return true;
}

bool PositionBasedDynamics::solve_CollisionConstraint(
        const Vector3r &p0, Real invMass0,
        const Vector3r &q1,
        const Vector3r &qn,
        Vector3r &corr0, double coll_EPS
){
    Vector3r qn_n= qn.normalized();
    Real c_p = (p0-q1).dot(qn_n);

    if(c_p < coll_EPS){
        corr0 = -c_p * qn_n ;
    }else{
        // no correction
        corr0(0) = 0;
        corr0(1) = 0;
        corr0(2) = 0;
    }
    return true;
}