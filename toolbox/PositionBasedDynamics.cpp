
#include "PositionBasedDynamics.h"
#include <Eigen/Dense>
#include "MathFunctions.h"
#include "constraint_utils.h"
#include <iostream>
#include "igl/barycentric_interpolation.h"
#include <igl/barycentric_coordinates.h>

const Real eps = static_cast<Real>(1e-6);

typedef double Real;
using namespace std;
using namespace Eigen;
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
    /* trial  no elongation, just shortening */
//    double clampElongation = max(-3.,(d - restLength) );
//    corr = stiffness * n * (clampElongation) / wSum;
    corr = stiffness * n * (d - restLength) / wSum;

    corr0 =  corr;
    corr1 = - corr;
    return true;
}
// ----------------------------------------------------------------------------------------------
bool PositionBasedDynamics::solveUVSimple(MatrixXi& Fg_pattern_curr, MatrixXi& mapFromFg, MatrixXd& p_adaption, MatrixXd& mapFromVg, double stiffness){
    for(int i=0; i<Fg_pattern_curr.rows(); i++){
        VectorXd from0 = mapFromVg.row(mapFromFg(i, 0));
        VectorXd from1 = mapFromVg.row(mapFromFg(i, 1));
        VectorXd from2 = mapFromVg.row(mapFromFg(i, 2));

        VectorXd now0 = p_adaption.row(Fg_pattern_curr(i, 0));
        VectorXd now1 = p_adaption.row(Fg_pattern_curr(i, 1));
        VectorXd now2 = p_adaption.row(Fg_pattern_curr(i, 2)).transpose();

        VectorXd fromBary = (from0 + from1 + from2)/3;
        VectorXd nowBary = (now0 + now1 + now2)/3;

        VectorXd u = fromBary; u(0)+= 1;
        VectorXd v = fromBary; v(1)+= 1;

        VectorXd uBary, vBary;
        igl::barycentric_coordinates(u.transpose(), from0.transpose(), from1.transpose(), from2.transpose(), uBary);
        igl::barycentric_coordinates(v.transpose(), from0.transpose(), from1.transpose(), from2.transpose(), vBary);

        VectorXd uNow = uBary(0) * now0 + uBary(1) * now1 + uBary(2) * now2;
        VectorXd vNow = vBary(0) * now0 + vBary(1) * now1 + vBary(2) * now2;

        double uDiff = (uNow-nowBary).norm();
        double vDiff = (vNow-nowBary).norm();
        if(i==0) cout<<now0.transpose() <<endl<<now1.transpose() <<endl<<now2.transpose()<<endl<<from2.transpose()<<" vert "<<Fg_pattern_curr(i, 2)<<endl ;
        if(i==0)cout<<"uDiff: "<<uDiff<<endl;
        double diff0 = abs(now0(0) - now1(0)) * (1 / uDiff - 1) ;
        double diff1 = abs(now1(0) - now2(0)) * (1 / uDiff - 1);
        double diff2 = abs(now2(0) - now0(0)) * (1 / uDiff - 1);

        if(i==0)cout<<diff0<<" "<<diff1<<" "<<diff2<<" the diffs" <<endl;
        double diffMax = max(max(diff0, diff1),diff2);
        if(diffMax == diff0){
            if (now0(0) < now1(0)){
                now0(0) -= diff0 / 2;
                now1(0) += diff0 / 2;

            } else{
                now1(0) -= diff0 / 2;
                now0(0) += diff0 / 2;
            }
        }else if (diffMax == diff1){
            if (now2(0) < now1(0)){
                now2(0) -= diff1 / 2;
                now1(0) += diff1 / 2;

            } else{
                now1(0) -= diff1 / 2;
                now0(0) += diff1 / 2;
            }
        }else{
            if (now0(0) < now2(0)){
                now0(0) -= diff2 / 2;
                now2(0) += diff2 / 2;

            } else{
                now2(0) -= diff2 / 2;
                now0(0) += diff2 / 2;
            }
        }





        if(i==0)cout<<now0<<" "<<now1<<" "<<now2<<" the nows" <<endl;
        if(i==0) cout<< "update vec "<<now0.transpose() - p_adaption.row(Fg_pattern_curr(i, 0))<<endl;

        p_adaption.row(Fg_pattern_curr(i, 0)) += stiffness * (now0.transpose() - p_adaption.row(Fg_pattern_curr(i, 0)));
        p_adaption.row(Fg_pattern_curr(i, 1)) += stiffness * (now1.transpose() - p_adaption.row(Fg_pattern_curr(i, 1)));
        p_adaption.row(Fg_pattern_curr(i, 2)) += stiffness * (now2.transpose()- p_adaption.row(Fg_pattern_curr(i, 2)));
        if(i==0) cout<< p_adaption.row(Fg_pattern_curr(i, 0))<<" is vert 0 "<<endl;
        if(i==0)cout<< p_adaption.row(Fg_pattern_curr(i, 1))<<" is vert 1 "<<endl;
        if(i==0)cout<< p_adaption.row(Fg_pattern_curr(i, 2))<<" is vert 2 "<<endl;

    }
}

bool PositionBasedDynamics::init_UVStretchPattern( const Vector2r& perFaceU, const Vector2r& perFaceV,
                                            const Eigen::MatrixXd& patternCoords,const Eigen::MatrixXd& targetPositions,
                                            Vector2r &tarUV0, Vector2r &tarUV1, Vector2r &tarUV2, int uORv, double DiagStiffness ){
    Eigen::MatrixXd Jnorm(2, 2);
    Jnorm.col(0)= perFaceU; Jnorm.col(1)= perFaceV;
    Vector2d Jn_0 =  Jnorm.col(0).normalized();
    Vector2d Jn_1 = Jnorm.col(1).normalized();
// J ist für original auf current, dh auch pattern coords muss aus original sein

    if(uORv==1|| uORv==11){// we want to relax u stretch, hence go in direction where u is not stretched i.e it is normalized
        Jnorm.col(0) = Jn_0;
        Jnorm.col(1) = Jn_1;
    }
    if(uORv ==11){
        cout<<Jnorm<< " jnorm" <<endl;
    }

    double angle = acos((Jn_0).dot(Jn_1));
    double deg = angle*180/M_PI;
    double delta = abs(90-deg)/2;

    Eigen::Matrix2d newRot= Eigen::MatrixXd::Identity(2, 2);
    newRot(0, 0)= cos(DiagStiffness * delta);
    newRot(1, 1) = newRot(0, 0);
    newRot(0, 1) = - sin (DiagStiffness * delta);
    newRot(1, 0) =  sin ( DiagStiffness * delta);

    if(uORv ==11){
//        cout<<" deg"<<deg<<endl;
    }
    if(deg<=90){
        Jnorm.col(0) = newRot.transpose() * Jnorm.col(0);
        Jnorm.col(1) = newRot * Jnorm.col(1);
    }else{
        Jnorm.col(1) = newRot.transpose() * Jnorm.col(1);
        Jnorm.col(0) = newRot * Jnorm.col(0);
    }

    Eigen::MatrixXd jacobiStretchedPattern = Jnorm * patternCoords;
    if(uORv ==11) cout<<endl<<patternCoords<<" pattern coords"<<endl ;
//    if(uORv ==11) cout<<endl<<jacobiStretchedPattern<<" jacobi coords"<<endl ;
//    if(uORv ==11) cout<<endl<<Jnorm<<" jacobian norm coords after angle spread"<<endl ;

    // compute rotation and translation of that matrix to best fit the original
    Eigen::MatrixXd R_est;
    Eigen::VectorXd T_est;

    procrustesWORef(jacobiStretchedPattern.transpose(), targetPositions.transpose(), R_est, T_est);

    Eigen::MatrixXd rotTargetPos =   R_est* jacobiStretchedPattern ;
    Eigen::MatrixXd refTargetPos = rotTargetPos.colwise() + T_est;
    if(uORv ==11) cout<<endl<<refTargetPos<<" rot trans coords"<<endl ;
    if( uORv ==11){
//        refTargetPos =  targetPositions;
//        tarUV1 =  targetPositions.col(1);
//        tarUV2 =  targetPositions.col(2);
    }
    tarUV0 = refTargetPos.col(0);
    tarUV1 = refTargetPos.col(1);
    tarUV2 = refTargetPos.col(2);


    return true;

}


bool PositionBasedDynamics::init_UVStretch( const Vector3r& perFaceU, const Vector3r& perFaceV,
                                           const Eigen::MatrixXd& patternCoords, const Eigen::Matrix3d& targetPositions,
                                           Vector3r &tarUV0, Vector3r &tarUV1, Vector3r &tarUV2, int uORv, double DiagStiffness ){

    Eigen::MatrixXd Jnorm(3, 2);
    Jnorm.col(0)= perFaceU; Jnorm.col(1)= perFaceV;
    Eigen::MatrixXd Jn_n = Jnorm;
    Vector3d Jn_0 =  Jnorm.col(0).normalized();
    Vector3d Jn_1 = Jnorm.col(1).normalized();
    if(uORv==1){// we want ro relax u stretch, hence go in direction where u is not stretched i.e it is normalized
        Jnorm.col(0) = Jnorm.col(0).normalized();// TODO NOT SO SURE IF THIS CAN ACTUALLY BE DONE
    }else {
        Jnorm.col(1) = Jnorm.col(1).normalized();
    }

    // test rotation of jacobian vectors
    double angle = acos((Jn_0).dot(Jn_1));
    Vector3d crossVec = Jn_0.cross(Jn_1);
    //TODO THIS IS A HEURISTIC AND DOES NOT WORK WITH NEGATIVE SIGNS
//    if(oldNormalVec.dot(crossVec)<0){
//        newAngle = -newAngle;
//    }
    double deg = angle*180/M_PI;
    double delta = abs(90-deg)/2;

    Eigen::Matrix4d newrotMat= Eigen::MatrixXd::Identity(4, 4);
    MathFunctions::setUpRotationMatrix(DiagStiffness * delta, crossVec, newrotMat);
    Matrix3d newnewRot = newrotMat.block(0,0,3,3);

    // change if negative , fixed 8.12.2022
    if(deg<=90){
        Jnorm.col(0) = newnewRot.transpose() * Jnorm.col(0);
        Jnorm.col(1) = newnewRot * Jnorm.col(1);
    }else{
        Jnorm.col(1) = newnewRot.transpose() * Jnorm.col(1);
        Jnorm.col(0) = newnewRot * Jnorm.col(0);
    }

//    Jnorm.col(0) = newnewRot.transpose() * Jnorm.col(0);
//    Jnorm.col(1) = newnewRot * Jnorm.col(1);
    Eigen::Matrix3d jacobiStretchedPattern = Jnorm * patternCoords;

    // compute rotation and translation of that matrix to best fit the original
    Eigen::MatrixXd R_est;
    Eigen::VectorXd T_est;

    // changed to  WO REF
    procrustesWORef(jacobiStretchedPattern.transpose(), targetPositions.transpose(), R_est, T_est);

    Eigen::Matrix3d rotTargetPos =   R_est* jacobiStretchedPattern ;
    Eigen::Matrix3d refTargetPos = rotTargetPos.colwise() + T_est;

    tarUV0 = refTargetPos.col(0);
    tarUV1 = refTargetPos.col(1);
    tarUV2 = refTargetPos.col(2);
    return true;

}


bool PositionBasedDynamics::init_Diag_Stretch( const Vector3r& perFaceU,const Vector3r& perFaceV,
                                           const Eigen::MatrixXd& patternCoords, const Eigen::Matrix3d& targetPositions,
                                           Vector3r &tarUV0, Vector3r &tarUV1, Vector3r &tarUV2){

    Eigen::MatrixXd Jnorm(3, 2);
    Jnorm.col(0)= perFaceU.normalized(); Jnorm.col(1)= perFaceV.normalized();

    Eigen::Matrix3d jacobiStretchedPattern = Jnorm * patternCoords;

    // compute rotation and translation of that matrix to best fit the original
    Eigen::MatrixXd R_est;
    Eigen::VectorXd T_est;
    procrustes(jacobiStretchedPattern.transpose(), targetPositions.transpose(), R_est, T_est);

    Eigen::Matrix3d rotTargetPos =   R_est* jacobiStretchedPattern ;
    Eigen::Matrix3d refTargetPos = rotTargetPos.colwise() + T_est;

    tarUV0 = refTargetPos.col(0);// cout<<tarUV0.rows()<< " rows and cols "<<tarUV0.cols()<<endl;
    tarUV1 = refTargetPos.col(1);
    tarUV2 = refTargetPos.col(2);
    return true;

}
bool PositionBasedDynamics::solve_Stretch(
        const Eigen::Matrix3d& targetPositions, const Vector3r &tarUV0, const Vector3r &tarUV1, const Vector3r &tarUV2,
        const Real stiffness,
        Vector3r &corr0, Vector3r &corr1, Vector3r &corr2 )
{

    Vector3r dir0 = tarUV0 - targetPositions.col(0) ;
    Vector3r dir1 = tarUV1 - targetPositions.col(1) ;
    Vector3r dir2 = tarUV2 - targetPositions.col(2) ;
//TODO WEIGHTING

    corr0 = stiffness * dir0;//.normalized(); //* abs(su-1) // or normalize dir and multiply by su
    corr1 = stiffness * dir1;//.normalized() ;
    corr2 = stiffness * dir2;//.normalized();
    return true;
}
bool PositionBasedDynamics::solve_Diag_Stretch(
       const Eigen::Matrix3d& targetPositions, const Vector3r &tarUV0, const Vector3r &tarUV1, const Vector3r &tarUV2,
        const Real stiffness,
        Vector3r &corr0, Vector3r &corr1, Vector3r &corr2 )
{

//    Vector3r dir0 = tarUV0 - targetPositions.col(0) ;
//    Vector3r dir1 = tarUV1 - targetPositions.col(1) ;
//    Vector3r dir2 = tarUV2 - targetPositions.col(2) ;
////TODO WEIGHTING
//
//    corr0 = stiffness * dir0;//.normalized();
//    corr1 = stiffness * dir1;//.normalized() ;
//    corr2 = stiffness * dir2;//.normalized();
    return true;
}

// ----------------------------------------------------------------------------------------------
bool PositionBasedDynamics::solve_RigidEnergy(const double& rigidEnergy, const double& rigidEPS, const double& rigidStiffness,
                                              const Vector3r& patternPos0,const Vector3r& patternPos1,const Vector3r& patternPos2,
                                              const Vector3r& p0, const Vector3r& p1, const Vector3r& p2,
                                               Vector3r& delta0, Vector3r& delta1, Vector3r& delta2){
    delta0 = Vector3r ::Zero(3);
    delta1 = Vector3r ::Zero(3);
    delta2 = Vector3r ::Zero(3);

    if(rigidEnergy > rigidEPS){
        // we are not rigid enough, thus take a step in direction of the best fit mapping
        Vector3r dir0 = patternPos0 - p0 ; //cout<<dir0.norm()<<endl;
        Vector3r dir1 = patternPos1 - p1 ;
        Vector3r dir2 = patternPos2 - p2 ;

        delta0 = rigidEnergy * rigidStiffness * dir0.normalized(); // is the rigid energy term even needed since we have it implicitly in the direction vector length (when not normalized)
        delta1 = rigidEnergy * rigidStiffness * dir1.normalized();
        delta2 = rigidEnergy * rigidStiffness * dir2.normalized();
    }
    return true;
}


// ----------------------------------------------------------------------------------------------
bool PositionBasedDynamics::solve_CollisionConstraint(
        const Vector3r &p0,
        const Vector3r &q1,
        const Vector3r &qn,
        Vector3r &corr0, double coll_EPS,const Vector3r & vel1
){
    Vector3r qn_n= qn.normalized();// that normalization does not make a difference
    Vector3r newq1 = q1 + coll_EPS*qn_n;
    Real c_p = (p0-newq1).dot(qn_n);
//    Real c_p = (p0-q1).dot(qn_n);
//    Real c_p = (p0-q1).normalized().dot(qn_n);

//&&  vel1.dot(qn)<0  makes it worse or at least no better
    if(c_p < coll_EPS ){
        corr0 = -c_p * qn_n ;
    }else{
        // no correction
        corr0(0) = 0;
        corr0(1) = 0;
        corr0(2) = 0;
    }
    return true;
}