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
        //cout<<"curr at j="<<j<<endl;
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

void garment_adaption::performJacobianUpdateAndMerge(Eigen::MatrixXd & V_curr){
    V= V_pattern;
   // this is ithe skrinked model, it has too much stress, we want to enlarge it
    // j-1 gives us the pattern such that the stress is just the same as for the original
    // but j-1 gives a different vertex position per face, thus we average it
    // so for each face we apply the inverse and get two new edges. these edges start from v_curr_0

    for (int numIt=0; numIt <1000; numIt++) {
        std::vector<std::vector<std::pair<Eigen::Vector3d, int>>> perVertexPositions(V_pattern.rows());
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


                //if(j==0 )cout<<oldNormalVec<<" old, "<<endl<< R*oldNormalVec<<" R*old "<<endl<<endl<<normalVec<<" should be=new normal"<<endl<<R<<" the rotation"<<endl;

                Eigen::MatrixXd jacobianAdapted = R*jacobians[j];
                inv_jacobians[j] = jacobianAdapted.inverse();
            }
            
            MatrixXd jacobi_adapted_Edge = inv_jacobians[j] * positions;
//            jacobi_adapted_Edge.row(2)= VectorXd::Zero(3);
            if(jacobi_adapted_Edge.row(2)(0)!=0 || jacobi_adapted_Edge.row(2)(1)!=0 || jacobi_adapted_Edge.row(2)(2)!=0 ){
//                cout<<jacobi_adapted_Edge<<" the jacobi "<<endl<<endl;
//                cout<<inv_jacobians[j]<<" inv jacobi "<<endl<<endl;

            }

            // now the reference is the barycenter of the 2D patter of the unshrinked model
            Eigen::Vector3d ref = (V.row(idp0)+ V.row(idp1)+ V.row(idp2))/3 ;

            perVertexPositions[idp0].push_back(std::make_pair(ref + jacobi_adapted_Edge.col(0), j));
            perVertexPositions[idp1].push_back(std::make_pair(ref + jacobi_adapted_Edge.col(1), j));
            perVertexPositions[idp2].push_back(std::make_pair(ref + jacobi_adapted_Edge.col(2), j));
//            if((ref + jacobi_adapted_Edge.col(2))(2)!= 200){
//                cout<<(ref + jacobi_adapted_Edge.col(2))(2)<<" ref "<<ref.transpose()<<endl<<jacobi_adapted_Edge.col(2).transpose()<<endl<<endl;
//            }

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



