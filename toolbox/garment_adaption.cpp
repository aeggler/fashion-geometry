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
//    perVertexPositions.resize(numVert); // per vert num size can differ, use push back since it should not be much
    V_init = Vg;
    V_pattern = V_pattern_orig;
    Fg_pattern = Fg_pattern_orig;


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
        normalVec= normalVec.normalized();

        jacobian.col(0) = jac2to3.col(0);
        jacobian.col(1) = jac2to3.col(1);
        jacobian.col(2) = normalVec;

        jacobians[j] = jacobian;

        inv_jacobians[j] = jacobian.inverse();
        if(j==50){
            int currid0 = Fg_pattern(j, 0);
            int currid1 = Fg_pattern(j, 1);
            auto edge1= V_pattern.row(currid1).transpose()-V_pattern.row(currid0).transpose();
            cout<<edge1 << " ... orig in uv and p1"<<p1<<endl <<endl<<endl;
            cout<<" now after jacobian"<<endl;
            cout<<jacobians[j]*edge1<<"...is the result"<<endl;
            cout<<"for comparison, the orig position in 3D "<<endl<<p1<<endl;
            cout<<endl<<endl;

            cout<<"now the inverse test"<<endl;
            cout<<inv_jacobians[j]*p1<<"...is the result"<<endl;
            cout<<u0+inv_jacobians[j]*p1 <<" should be "<<  V_pattern.row(currid1)<<endl;
        }

    }
cout<<"fin"<<endl;
};

void garment_adaption::performJacobianUpdateAndMerge(Eigen::MatrixXd & V_curr){
    V= V_pattern;
   // this is ithe skrinked model, it has too much stress, we want to enlarge it
    // j-1 gives us the pattern such that the stress is just the same as for the original
    // but j-1 gives a different vertex position per face, thus we average it
    // so for each face we apply the inverse and get two new edges. these edges start from v_curr_0

    for (int numIt=0; numIt <1; numIt++) {
        std::vector<std::vector<Eigen::Vector3d>> perVertexPositions(V_pattern.rows());
        for (int j = 0; j < numFace; j++) {
            int id0 = Fg(j, 0); int idp0 = Fg_pattern(j, 0);
            int id1 = Fg(j, 1); int idp1 = Fg_pattern(j, 1);
            int id2 = Fg(j, 2); int idp2 = Fg_pattern(j, 2);
            Eigen::MatrixXd positions(3, 2);
            positions.col(0) = V_curr.row(id1) - V_curr.row(id0);
            positions.col(1) = V_curr.row(id2) - V_curr.row(id0);
//            cout<<j<<" j 0 "<<endl;

            MatrixXd jacobi_adapted_Edge = inv_jacobians[j] * positions;
            Eigen::Vector3d ref = V.row(idp0);// no reason why it should start here, initial guess
//            cout<<j<<"=j "<<numFace<<endl;

            perVertexPositions[idp0].push_back(ref);
//            cout<<j<<" j a "<<endl;
            perVertexPositions[idp1].push_back(ref + jacobi_adapted_Edge.col(0));
//            cout<<j<<" j b "<<idp2<<" per vert pos"<<perVertexPositions.size()<<endl;

            perVertexPositions[idp2].push_back(ref + jacobi_adapted_Edge.col(1));
//            cout<<j<<" j c "<<endl;

        }
cout<<"test"<<endl;
        for (int i = 0; i < numVert; i++) {

            Eigen::Vector3d avg = Eigen::VectorXd::Zero(3);

            for (int j = 0; j < perVertexPositions[i].size(); j++) {
                avg += perVertexPositions[i][j];
                if (i == 50) cout << perVertexPositions[i][j] << endl << endl;
            }
            avg /= perVertexPositions[i].size();// not surre if size is the right one
            if (i == 50) cout << avg << " average position" << endl;

            V.row(i) = avg.transpose(); //(V_curr.row(i)+ avg.transpose())/2; // average curr position and new computed
            //V.row(i) /= 2;
            if (i == 50) cout << V.row(i) << endl;
        }
        // set the curr position as the positions we just computed
    }
    V_curr = V;
}



