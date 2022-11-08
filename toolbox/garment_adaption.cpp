//
// Created by Anna Maria Eggler on 31.10.22.
//

#include "garment_adaption.h"
#include <Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;

garment_adaption::garment_adaption(Eigen::MatrixXd &Vg, Eigen::MatrixXi& Fg, Eigen::MatrixXd & V_pattern_orig, Eigen::MatrixXi& Fg_pattern_orig ) {
    numFace= Fg.rows();
    numVert = Vg.rows();
    this->Fg = Fg;
    jacobians.resize(numFace); //    std::vector<Eigen::MatrixXd > jacobians; 3x3  from 2 to 3 for 2 x 1 positions, thus col wise
    inv_jacobians.resize(numFace);// 3 x 3 each
    V_init = Vg;
    V_pattern = V_pattern_orig;
    Fg_pattern = Fg_pattern_orig;


}
void Barycentric(VectorXd& p, VectorXd a, VectorXd b, VectorXd c, VectorXd& baryP)
{
    Vector3d v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = v0.dot( v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot( v1);
    float d20 = v2.dot( v0);
    float d21 = v2.dot( v1);
    float denom = d00 * d11 - d01 * d01;
    baryP(0) = (d11 * d20 - d01 * d21) / denom;
    baryP(1) = (d00 * d21 - d01 * d20) / denom;
    baryP(2) = 1.0f - baryP(0) - baryP(1);
}
void garment_adaption::computeJacobian(){
    for(int j = 0; j<numFace; j++){
        Eigen::MatrixXd jac2to3 (3, 2);
        Eigen::MatrixXd jacobian(3, 3);
        Eigen::MatrixXd inv_jac(2, 3);

        int id0 = Fg_pattern(j, 0);
        int id1 = Fg_pattern(j, 1);
        int id2 = Fg_pattern(j, 2);

        Vector2d u0, u1h, u2h;
        u0(0) = V_pattern(id0, 0);
        u0(1) = V_pattern(id0, 1);

        u1h( 0) = V_pattern(id1, 0);
        u1h(1) = V_pattern(id1, 1);

        u2h( 0) = V_pattern(id2, 0);
        u2h(1) = V_pattern(id2, 1);

        u1h -= u0;
        u2h -= u0;

        double det = u1h( 0)*u2h(1)- (u2h(0)*u1h(1));
        u1h/=det;
        u2h/=det;

        // attention reuse
         id0 = Fg(j, 0);
         id1 = Fg(j, 1);
         id2 = Fg(j, 2);

        Eigen::Vector3d faceAvg = (V_init.row(id0) + V_init.row(id1) + V_init.row(id2)) / 3;

        Vector3d p0 = V_init.row(id0);
        Vector3d p1 = V_init.row(id1);
        Vector3d p2 = V_init.row(id2);

        p1 -= p0;
        p2 -= p0;

        // from Nicos code
        jac2to3.col(0) = p1 * u2h(1) - p2 * u1h(1);
        jac2to3.col(1) = p2 * u1h( 0) - p1 * u2h( 0);

        Vector3d normalVec = p1.cross(p2);
        normalVec = normalVec.normalized();//normalVec /= 200;//

        jacobian.col(0) = jac2to3.col(0);
        jacobian.col(1) = jac2to3.col(1);
        jacobian.col(2) = normalVec;


       // cout<<"ortho test "<< jacobian.col(0).dot(jacobian.col(2))<<" and other "<< jacobian.col(1).dot(jacobian.col(2))<<endl;

        jacobians[j] = jacobian;

        inv_jacobians[j] = jacobian.inverse();

    }
}

void garment_adaption::setUpRotationMatrix(double angle,Vector3d& axis, Matrix4d& rotationMatrix)
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
}

void garment_adaption::performJacobianUpdateAndMerge(Eigen::MatrixXd & V_curr, int iterations){
    V= V_pattern;
   // this is ithe skrinked model, it has too much stress, we want to enlarge it
    // j-1 gives us the pattern such that the stress is just the same as for the original
    // but j-1 gives a different vertex position per face, thus we average it
    // so for each face we apply the inverse and get two new edges. these edges start from v_curr_0

    std::vector<Eigen::Matrix3d> jacobi_adapted_Edges (numFace);

    for (int j = 0; j < numFace; j++) {
        int id0 = Fg(j, 0); int idp0 = Fg_pattern(j, 0);
        int id1 = Fg(j, 1); int idp1 = Fg_pattern(j, 1);
        int id2 = Fg(j, 2); int idp2 = Fg_pattern(j, 2);
        RowVector3d barycenter = V_curr.row(id1) + V_curr.row(id0) + V_curr.row(id2);
        barycenter/=3;

        // from center to all adjacent vertices
        Eigen::MatrixXd positions(3, 3);
        positions.col(0) = V_curr.row(id0) - barycenter;
        positions.col(1) = V_curr.row(id1) - barycenter;
        positions.col(2) = V_curr.row(id2) - barycenter;

        /*test: we take the normal of the current and of the original triangle and align the original to the current
            * , the same rotation is then also applied to the jacobian*/
        Vector3d p2 = V_curr.row(id2) - V_curr.row(id0);
        Vector3d p1 =V_curr.row(id1) - V_curr.row(id0);
        Vector3d normalVec = p1.cross(p2);
        normalVec= normalVec.normalized(); // the new normal vector, how do we get from jacobians[j].col(2) to this?

        Vector3d oldNormalVec =  jacobians[j].col(2);
        // rotate old to normalVec
        if(oldNormalVec!= normalVec){
            //R = I + [v]_x + [v]_x^2 * (1-c)/s^2
            //https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/180436#180436
            Eigen::VectorXd v = oldNormalVec.cross(normalVec);
            double c = oldNormalVec.dot(normalVec);

            // I am sure this can be made faster/simpler but this is a logical derivation
            Eigen::MatrixXd G= Eigen::MatrixXd::Zero(3, 3);
            G(0,0) = c; G(1, 1)= c; G(2, 2)= 1;
            G(0, 1) = -v.norm();
            G(1, 0) = v.norm();

            Eigen::MatrixXd F= Eigen::MatrixXd::Zero(3, 3);
            F.col(0)= oldNormalVec;
            F.col(1) = (normalVec- c*oldNormalVec).normalized();
            F.col(2) = normalVec.cross(oldNormalVec);
            MatrixXd R = F* G*F.inverse();

            // trial - normalization or what is missing? rotate both and then search edge ?
            Vector3d newp2 = ((V_curr.row(id2) - V_curr.row(id0)));

            // we have aligned the normal already (the 3rd col) so it should be a 2D rotation
            // https://math.stackexchange.com/questions/3563901/how-to-find-2d-rotation-matrix-that-rotates-vector-mathbfa-to-mathbfb
            VectorXd zeroRotated = R*(V_init.row(id0).transpose());
            VectorXd twoRotated = R*(V_init.row(id2).transpose());
            VectorXd oldp2Rotated = (twoRotated- zeroRotated);
        // what is the angle between p2 and oldp2Rotated? theyy are normal aligned so it should be a simple rotation around normal axis
            auto angle = acos((p2.normalized()).dot(oldp2Rotated.normalized()));
           if(j==50) cout<<angle<<" angle"<<endl;
            //
            auto degree = angle*180/M_PI;
            Eigen::Matrix4d rotMat= Eigen::MatrixXd::Identity(4, 4);
            setUpRotationMatrix(degree,normalVec, rotMat);
            if(j==50) cout<<rotMat<<" rot mat "<<degree<<endl;



//            // how much we want to rotate about the axis, but we need ot rotate there first and back after
//            Eigen::MatrixXd D= MatrixXd::Identity(4,4);
//            Eigen::MatrixXd D_inv= MatrixXd::Identity(4,4);
//
//            D(0, 3)= -V_curr(id0, 0);
//            D(1, 3)= -V_curr(id0, 1);//-zeroRotated(1);
//            D(2, 3)= -V_curr(id0, 2);//-zeroRotated(2);
//            D_inv = D;
//            D_inv(0, 3)= V_curr(id0, 0);//zeroRotated(0);
//            D_inv(1, 3)= V_curr(id0, 1);//zeroRotated(1);
//            D_inv(2, 3)= V_curr(id0, 2);//zeroRotated(2);
//
//            Eigen::MatrixXd Rx= MatrixXd::Identity(4,4);
//            Eigen::MatrixXd Rx_inv= MatrixXd::Identity(4,4);
//
//            double vval = sqrt(normalVec(1)*normalVec(1)+normalVec(2)*normalVec(2));
//            Rx(1, 1)= normalVec(2)/vval;
//            Rx(2, 2)= Rx(1, 1);
//            Rx(1, 2) = -normalVec(1)/vval;
//            Rx(2, 1) = -Rx(1, 2);
//            Rx_inv = Rx;
//            Rx_inv(2, 1) = -Rx(2, 1);
//            Rx_inv(1, 2) = -Rx(1, 2);
//
//            Eigen::MatrixXd Ry= MatrixXd::Identity(4,4);
//            Eigen::MatrixXd Ry_inv= MatrixXd::Identity(4,4);
//
//            double lnorm = normalVec.norm();
//            Ry(0, 0) = vval/lnorm;
//            Ry(2, 2) = Ry(0,0);
//            Ry(0, 2) = -normalVec(0)/lnorm;
//            Ry(2, 0) = -R(0, 2);
//            Ry_inv = Ry;
//            Ry_inv(0, 2)= -Ry(0, 2);
//            Ry_inv(2, 0)= -Ry(2, 0);


//            MatrixXd rot2d   = MatrixXd::Identity(4, 4);
//            rot2d(0, 0) = p2(0)*oldp2(0) + p2(1)*oldp2(1);
//            rot2d(1, 1) = rot2d(0,0);
//            rot2d(0,1) = p2(0)* oldp2(1) - oldp2(0) * p2(1);
//            rot2d (1, 0) = oldp2(0)*p2(1) - p2(0)* oldp2(1);
//
//            Eigen::MatrixXd T = D_inv * Rx_inv * Ry_inv * rot2d * Ry * Rx *D;

            Eigen::MatrixXd jacobianAdapted =  R*jacobians[j];
            Vector4d j1; j1<<jacobianAdapted.col(0), 1;
            Vector4d j2; j2<<jacobianAdapted.col(1), 1;
            VectorXd outputj1=VectorXd::Zero(4);
            VectorXd outputj2=VectorXd::Zero(4);
            Vector4d j3; j3<<jacobianAdapted.col(2), 1;
            VectorXd outputj3=VectorXd::Zero(4);
            for(int ii=0; ii<4; ii++){
                for(int k =0; k<4; k++){
                    outputj1(ii) += rotMat(ii, k)* j1(k);
                }
            }
            for(int ii=0; ii<4; ii++){
                for(int k =0; k<4; k++){
                    outputj2(ii) += rotMat(ii, k)* j2(k);
                }
            }
            for(int ii=0; ii<4; ii++){
                for(int k =0; k<4; k++){
                    outputj3(ii) += rotMat(ii, k)* j3(k);
                }
            }


            Eigen::MatrixXd adapted(3,3); adapted.col(2)= jacobianAdapted.col(2);

            if(j==50) cout<<outputj3<<" out j3 "<<endl<<endl<<j3<<" and not rot norm " <<endl;
            adapted(0, 0)= outputj1(0);
            adapted(1, 0)= outputj1(1);
            adapted(2, 0)= outputj1(2);

            adapted(0, 1)= outputj2(0);
            adapted(1, 1)= outputj2(1);
            adapted(2, 1)= outputj2(2);

            inv_jacobians[j] = adapted.inverse();
//            if(j==50 ) {
//                cout << jacobianAdapted << " the normal aligned" << endl;
//                cout << endl << rot2d << " rotation around normal" << endl << rot2d * jacobianAdapted << endl << endl
//                     << adapted << endl;
//            }
//            Eigen::MatrixXd jacobianAdapted = R*jacobians[j];
//            inv_jacobians[j] = jacobianAdapted.inverse();

        }

        MatrixXd jacobi_adapted_Edge = inv_jacobians[j] * positions;
        jacobi_adapted_Edges[j]= jacobi_adapted_Edge;


    }

    for (int numIt=0; numIt <iterations; numIt++) {
        std::vector<std::vector<std::pair<Eigen::Vector3d, int>>> perVertexPositions(V_pattern.rows());

        for(int j=0; j<numFace; j++){
            // now the reference is the barycenter of the 2D patter of the unshrinked model
            int idp0 = Fg_pattern(j, 0);
            int idp1 = Fg_pattern(j, 1);
            int idp2 = Fg_pattern(j, 2);

            Eigen::Vector3d ref = (V.row(idp0)+ V.row(idp1)+ V.row(idp2))/3 ;

            perVertexPositions[idp0].push_back(std::make_pair(ref + jacobi_adapted_Edges[j].col(0), j));
            perVertexPositions[idp1].push_back(std::make_pair(ref + jacobi_adapted_Edges[j].col(1), j));
            perVertexPositions[idp2].push_back(std::make_pair(ref + jacobi_adapted_Edges[j].col(2), j));
        }

        // this iteration makes no sense, we average back to the original.
        // instead we should fix one vertex and from there on fix all the others.
        for (int i = 0; i < numVert; i++) {

            Eigen::Vector3d avg = Eigen::VectorXd::Zero(3);

            for (int j = 0; j < perVertexPositions[i].size(); j++) {
                Eigen::Vector3d curr = get<0>(perVertexPositions[i][j]);
                avg += curr;
            }

            avg /= perVertexPositions[i].size();

            V.row(i) = avg.transpose();
        }

    }
    V_curr = V;
}




