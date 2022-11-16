
#ifndef EXAMPLE_POSITIONBASEDDYNAMICS_H
#define EXAMPLE_POSITIONBASEDDYNAMICS_H
#include <Eigen/Dense>

typedef double Real;

using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
using Matrix4r = Eigen::Matrix<Real, 4, 4, Eigen::DontAlign>;
using Matrix2r = Eigen::Matrix<Real, 2, 2, Eigen::DontAlign>;

class PositionBasedDynamics {
public:
    static bool init_IsometricBendingConstraint(
            const Vector3r &p0,
            const Vector3r &p1,
            const Vector3r &p2,
            const Vector3r &p3,
            Matrix4r &Q
    );

    static bool solve_IsometricBendingConstraint(
            const Vector3r &p0,
            Real invMass0,        // angle on (p2, p3) between triangles (p0, p2, p3) and (p1, p3, p2)
            const Vector3r &p1, Real invMass1,
            const Vector3r &p2, Real invMass2,
            const Vector3r &p3, Real invMass3,
            const Matrix4r &Q,
            const Real stiffness,
            Vector3r &corr0, Vector3r &corr1, Vector3r &corr2, Vector3r &corr3);

    /** Determine the position corrections for a distance constraint between two particles:\n\n
		* \f$C(\mathbf{p}_0, \mathbf{p}_1) = \| \mathbf{p}_0 - \mathbf{p}_1\| - l_0 = 0\f$\n\n
		* More information can be found in the following papers: \cite Mueller07, \cite BMOT2013, \cite BMOTM2014, \cite BMM2015,
		*
		* @param p0 position of first particle
		* @param invMass0 inverse mass of first particle
		* @param p1 position of second particle
		* @param invMass1 inverse mass of second particle
		* @param restLength rest length of distance constraint
		* @param stiffness stiffness coefficient
		* @param corr0 position correction of first particle
		* @param corr1 position correction of second particle
		*/
    static bool solve_DistanceConstraint(
            const Vector3r &p0, Real invMass0,
            const Vector3r &p1, Real invMass1,
            const Real restLength,
            const Real stiffness,
            Vector3r &corr0, Vector3r &corr1);

    bool init_UVStretch(const Real su, const Vector3r& perFaceU,const Vector3r& perFaceV,
                                               const Eigen::MatrixXd& patternCoords, const Eigen::Matrix3d& targetPositions,
                                               const Real stiffness,
                        Vector3r &tarUV0, Vector3r &tarUV1, Vector3r &tarUV2, int uORv);
    // test approach for the new functions
    static bool solve_UVStretch(
            const Real su, const Eigen::Matrix3d& targetPositions,int uORv, const Vector3r &tarUV0, const Vector3r &tarUV1, const Vector3r &tarUV2,
            const Real stiffness,
            Vector3r &corr0, Vector3r &corr1, Vector3r &corr2);

    static bool solve_RigidEnergy(
            const double& rigidEnergy, const double& rigidEPS, const double& rigidStiffness,
            const Vector3r& patternPos0,const Vector3r& patternPos1,const Vector3r& patternPos2,
            const Vector3r& p0, const Vector3r& p1, const Vector3r& p2,
            Vector3r& delta0, Vector3r& delta1, Vector3r& delta2);


    static bool solve_CollisionConstraint(
            const Vector3r &p0,
            const Vector3r &q1,
            const Vector3r &qn,
            Vector3r &corr0, double coll_EPS
            ,const Vector3r & vel1
            );


};
#endif //EXAMPLE_POSITIONBASEDDYNAMICS_H
