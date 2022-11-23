//
// Created by Anna Maria Eggler on 31.10.22.
//

#include "garment_adaption.h"
#include <Eigen/Dense>
#include <iostream>
#include <igl/doublearea.h>
#include <cmath>
#include "adjacency.h"
#include<Eigen/SparseCholesky>

using namespace std;
using namespace Eigen;

typedef Eigen::SparseMatrix<double, RowMajor> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

garment_adaption::garment_adaption(Eigen::MatrixXd &Vg, Eigen::MatrixXi& Fg, Eigen::MatrixXd & V_pattern_orig, Eigen::MatrixXi& Fg_pattern_orig ) {
    numFace= Fg.rows();
    numVertGarment = Vg.rows();
    numVertPattern = V_pattern_orig.rows();
    this->Fg = Fg;
    jacobians.resize(numFace); //    std::vector<Eigen::MatrixXd > jacobians; 3x3  from 2 to 3 for 2 x 1 positions, thus col wise
    inv_jacobians.resize(numFace);// 3 x 3 each
    V_init = Vg;
    V_pattern = V_pattern_orig;
    Fg_pattern = Fg_pattern_orig;
    createVertexFaceAdjacencyList(Fg_pattern, vfAdj);

    int rowCount = 0;
    std::vector<T> tripletList;

    for(int i =0; i<numVertPattern; i++){
        vector<int> neigh = vfAdj[i]; // number of faces of this vertex ,attention remove -1 faces
        //sort(neigh.begin(), neigh.end());
        for (int j = 0; j<neigh.size(); j++){
            if(neigh[j]== -1){
                continue;
            }
            tripletList.push_back(T(3 * rowCount, 3 * i, 1));
            tripletList.push_back(T(3 * rowCount + 1, 3 * i + 1, 1));
            tripletList.push_back(T(3 * rowCount + 2, 3 * i + 2, 1));
            rowCount++;

        }
    }
    A.resize(3*rowCount, 3*numVertPattern);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    SpMat LHS = (A.transpose()*A);

    cholSolver.analyzePattern(LHS);
    cholSolver.factorize(LHS);


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
    perFaceTargetNorm.resize(numFace);
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

        jacobians[j] = jacobian;
//        cout<<jacobians[j]<<" face j "<<endl<<endl;

        inv_jacobians[j] = jacobian.inverse();

        //they should have stretch 1, hence add deviation from 1 to measure
        perFaceTargetNorm(j) = (jacobian.col(0).norm()-1) * (jacobian.col(0).norm()-1);
        perFaceTargetNorm(j) += (jacobian.col(1).norm()-1) * (jacobian.col(1).norm()-1);
        // add some kind of angle measure
        // they should be orthogonal, hence add dot squared as norm
        double dot = jacobian.col(0).normalized().dot(jacobian.col(1).normalized());
        perFaceTargetNorm(j) += (dot * dot);

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

void garment_adaption::performJacobianUpdateAndMerge(Eigen::MatrixXd & V_curr, int iterations, const MatrixXd& baryCoords1, const MatrixXd& baryCoords2,Eigen::MatrixXd & V_newPattern ){
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
        MatrixXd R;
        if(oldNormalVec!= normalVec) {
            //R = I + [v]_x + [v]_x^2 * (1-c)/s^2
            //https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/180436#180436
            Eigen::VectorXd v = oldNormalVec.cross(normalVec);
            double c = oldNormalVec.dot(normalVec);

            // I am sure this can be made faster/simpler but this is a logical derivation
            Eigen::MatrixXd G = Eigen::MatrixXd::Zero(3, 3);
            G(0, 0) = c;
            G(1, 1) = c;
            G(2, 2) = 1;
            G(0, 1) = -v.norm();
            G(1, 0) = v.norm();

            Eigen::MatrixXd F = Eigen::MatrixXd::Zero(3, 3);
            F.col(0) = oldNormalVec;
            F.col(1) = (normalVec - c * oldNormalVec).normalized();
            F.col(2) = normalVec.cross(oldNormalVec);
            R = F * G * F.inverse();
        }else {
            R= MatrixXd::Identity(3, 3);
        }
            // rotate normalVec to old
            MatrixXd Rinv = R.inverse();
            // Update: we do not rotate the jacobian but actually the positions, then from these rotated positions
            // we align the u or v axis using the bary coordinates and difference between reference 3D of the original and our new 3D
            VectorXd zeroInv = Rinv * V_curr.row(id0).transpose();
            VectorXd oneInv = Rinv * V_curr.row(id1).transpose();
            VectorXd twoInv = Rinv * V_curr.row(id2).transpose();

            VectorXd currU = baryCoords1(j, 0)*zeroInv + baryCoords1(j, 1)* oneInv +baryCoords1(j, 2)*twoInv ;
            VectorXd currV = baryCoords2(j, 0)*zeroInv + baryCoords2(j, 1)* oneInv +baryCoords2(j, 2)*twoInv ;

            VectorXd oldU = baryCoords1(j, 0)*V_init.row(id0).transpose() + baryCoords1(j, 1)* V_init.row(id1).transpose() +baryCoords1(j, 2)*V_init.row(id2).transpose() ;
            VectorXd oldV = baryCoords2(j, 0)*V_init.row(id0).transpose() + baryCoords2(j, 1)* V_init.row(id1).transpose() +baryCoords2(j, 2)*V_init.row(id2).transpose() ;
            currU -=((zeroInv+oneInv+twoInv)/3.);
            oldU -= (V_init.row(id0)+V_init.row(id1)+V_init.row(id2))/3;

            // if u aligned , else other, can later be made user input
            Vector3d alignFrom = currU.normalized();
            Vector3d alignTo = oldU.normalized();

            // they are normal aligned so it should be a simple rotation around normal axis
            //https://www.euclideanspace.com/maths/algebra/vectors/angleBetween/ and https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
            // praise stackoverflow
            double newAngle = acos((alignFrom).dot(alignTo));
            auto crossVec = alignFrom.cross(alignTo); // or other way round?
            if(oldNormalVec.dot(crossVec)<0){
                newAngle = -newAngle;
            }

            double newdegree = newAngle*180/M_PI;

        // https://math.stackexchange.com/questions/3563901/how-to-find-2d-rotation-matrix-that-rotates-vector-mathbfa-to-mathbfb
            Eigen::Matrix4d newrotMat= Eigen::MatrixXd::Identity(4, 4);
            setUpRotationMatrix(newdegree,oldNormalVec, newrotMat);

            Eigen::MatrixXd jacobianAdapted =   jacobians[j];
            Matrix3d newnewRot = newrotMat.block(0,0,3,3);

        // when applied to the positions we rotate them back to align the normals,then we apply the jacobian inverse
        //then we apply the rotation around the normal
        // and finally the jinverse
        MatrixXd jacobi_adapted_Edge = inv_jacobians[j] * newnewRot * Rinv * positions;
        jacobi_adapted_Edges[j]= jacobi_adapted_Edge;

    }

    Eigen::VectorXd v_asVec (3*numVertPattern);
    Eigen::VectorXd v_asVecOld = VectorXd::Zero(3*numVertPattern);

    for(int i = 0; i<numVertPattern; i++){
        v_asVec(3 * i) = V(i, 0);
        v_asVec(3*i+1) = V(i, 1);
        v_asVec(3*i+2) = V(i, 2);
    }

    Eigen::VectorXd b (A.rows());

    for (int numIt=0; numIt< iterations ; numIt++) {
        std::vector<std::vector<std::pair<Eigen::Vector3d, int>>> perVertexPositions(numVertPattern);

        // the local step
        for(int j=0; j<numFace; j++){
            // now the reference is the barycenter of the 2D patter of the unshrinked model
            int idp0 = Fg_pattern(j, 0);
            int idp1 = Fg_pattern(j, 1);
            int idp2 = Fg_pattern(j, 2);

            Eigen::Vector3d ref = (v_asVec.block(3*idp0, 0, 3, 1)+
                    v_asVec.block(3*idp1, 0, 3, 1)+
                    v_asVec.block(3*idp2, 0, 3, 1))/3 ;

            perVertexPositions[idp0].push_back(std::make_pair(ref + jacobi_adapted_Edges[j].col(0), j));
            perVertexPositions[idp1].push_back(std::make_pair(ref + jacobi_adapted_Edges[j].col(1), j));
            perVertexPositions[idp2].push_back(std::make_pair(ref + jacobi_adapted_Edges[j].col(2), j));
        }

       // the global step
        int counter = 0;
        for(int j =0; j<numVertPattern; j++) {
            vector<std::pair<Eigen::Vector3d, int>> vi_pos = perVertexPositions[j];

            for (int curr = 0; curr < vi_pos.size(); curr++) {
                Eigen::Vector3d currPos = get<0>(vi_pos[curr]);
                b(3 * counter) = currPos(0);
                b(3 * counter + 1) = currPos(1);
                b(3 * counter + 2) = currPos(2);
                counter++;
            }
        }
        MatrixXd RHS = A.transpose() * b;
        v_asVecOld= v_asVec;
        v_asVec = cholSolver.solve(RHS);


        // not claiming this was very efficient ,just for debug purposes
        auto diff = (v_asVecOld-v_asVec);
        MatrixXd norms(numVertPattern, 3);
        for(int i=0; i<numVertPattern; i++){
            norms.row(i) = diff.block(3*i, 0, 3, 1).transpose();
        }
      //needs a more principled approach
        double maxPatternNorm = V_pattern.rowwise().norm().maxCoeff();
        double maxChange = norms.rowwise().norm().maxCoeff();
        cout<<"change in percent  "<<(maxChange/maxPatternNorm)<<endl;
        // TODO SET THIS FINISH PARAMETER
        if(maxChange/maxPatternNorm < 0.0001){
            cout<<"fin in iteration " <<numIt<<endl;
            break;
        }
    }

    V_newPattern.resize(numVertPattern, 3);

    for(int i=0; i<numVertPattern; i++){
        V_newPattern.row(i) = v_asVec.block(3*i, 0, 3, 1).transpose();
    }

}




