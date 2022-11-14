////
//// Created by Anna Maria Eggler on 26.10.22.
////

#include "body_interpolation.h"

#include <igl/polar_dec.h>
#include <igl/svd3x3.h>
#include <igl/AABB.h>
#include <igl/centroid.h>
#include <igl/sparse_cached.h>

#include <iostream>


BodyInterpolator::BodyInterpolator(MatrixXd& curr_V_from, MatrixXd& curr_V_to, MatrixXi& curr_F ){
    V_from = curr_V_from;
    V_to = curr_V_to;
    F = curr_F;
    cout<<"in setup"<<endl;

    int numFaces = F.rows();
    int numVert = V_from.rows();
    if(V_to.rows()!= numVert){
        // we run into a huge problem
        //TODO
    }
    U.resize(numFaces);
    U2.resize(numFaces);
    K.resize(numFaces);
    q.resize(numFaces);

    // just like katja we find the closest vertex
    double smallest_dist = numeric_limits<double>::max();
    fixed_id = -1;

    for(int i=0; i<numVert; i++){
        double dist = (V_from.row(i)- V_to.row(i)).norm();
        if(dist<smallest_dist){
            smallest_dist = dist;
            fixed_id = i;
        }
    }
    // set S and Q
    for(int f =0; f < numFaces; f++){
        Vector3d v1 = V_from.row(F(f, 0));
        Vector3d v2 = V_from.row(F(f, 1));
        Vector3d v3 = V_from.row(F(f, 2));
        // calculate the 4th vertex for this face
        Vector3d a = (v2 - v1).cross(v3 - v1);
        Vector3d v4 = v1 + a / sqrt(a.norm());

        Vector3d v1_def = V_to.row(F(f, 0));
        Vector3d v2_def = V_to.row(F(f, 1));
        Vector3d v3_def = V_to.row(F(f, 2));
        Vector3d a_def = (v2_def - v1_def).cross(v3_def - v1_def);
        Vector3d v4_def = v1_def + a_def / sqrt(a_def.norm());

        // calculate matrices V, V_def and Q
        Matrix3d V, V_def, Q;
        V.col(0) = v2 - v1;
        V.col(1) = v3 - v1;
        V.col(2) = v4 - v1;
        V_def.col(0) = v2_def - v1_def;
        V_def.col(1) = v3_def - v1_def;
        V_def.col(2) = v4_def - v1_def;
        Q = V_def * V.inverse();

        // polar decomposition
        Matrix3d R, S;
        igl::polar_dec(Q, R, S);

        // S= U *K*U2
        Matrix3d Uf, U2f;
        Vector3d Kf;
        igl::svd3x3(S, Uf, Kf, U2f);

        U[f] = Uf;
        U2[f] = U2f;
        K[f] = Kf;

        Quaterniond qf(R);
        q[f] = qf;
    }

    cout<<" after svd"<<endl;
    // Precompute LHS
    // since A is a sparse matrix, we set it from triplets

    vector< Triplet<double>> tri;
    tri.reserve(9 * numFaces);
    A_remcol = VectorXd::Zero(3 * numFaces);

    for(int f = 0; f<numFaces; f++){
        MatrixXd W(3, 2);
        W.col(0) = V_from.row(F(f, 1)) - V_from.row(F(f, 0));
        W.col(1) = V_from.row(F(f, 2)) - V_from.row(F(f, 0));

        HouseholderQR<MatrixXd> qr(W);
        Matrix3d Q_qr = qr.householderQ();
        MatrixXd R_qr = Q_qr.transpose() * W;
        MatrixXd T = R_qr.block<2, 2>(0, 0).inverse() * Q_qr.block<3, 2>(0, 0).transpose();

        // fill A - refer to page 65 of Robert Sumners Thesis
        // check each time if we are in the removed column instead
        int i1 = F(f, 0);
        int i2 = F(f, 1);
        int i3 = F(f, 2);

        if (i2 != fixed_id) {
            if (i2 > fixed_id) i2 -= 1;
            for (int i = 0; i < 3; i++)
                tri.push_back(Triplet<double>(3 * f + i, i2, T(0, i)));
        }
        else
            for (int i = 0; i < 3; i++)
                A_remcol(3 * f + i) = T(0, i);
        if (i3 != fixed_id) {
            if (i3 > fixed_id) i3 -= 1;
            for (int i = 0; i < 3; i++)
                tri.push_back(Triplet<double>(3 * f + i, i3, T(1, i)));
        }
        else
            for (int i = 0; i < 3; i++)
                A_remcol(3 * f + i) = T(1, i);
        if (i1 != fixed_id) {
            if (i1 > fixed_id) i1 -= 1;
            for (int i = 0; i < 3; i++)
                tri.push_back(Triplet<double>(3 * f + i, i1, -T(0, i) - T(1, i)));
        }
        else
            for (int i = 0; i < 3; i++)
                A_remcol(3 * f + i) = -T(0, i) - T(1, i);

    }
    cout<<"setting the A mesh now"<<endl;
    //A.setFromTriplets(tri.begin(), tri.end());
    if (A_data.rows() == 0) {
        A = SparseMatrix<double>(3 * numFaces, numVert - 1);
        igl::sparse_cached_precompute(tri, A_data, A);
    }
    else
        igl::sparse_cached(tri, A_data, A);
    ATA = A.transpose() * A;
    cholSolver.analyzePattern(ATA);
    cholSolver.factorize(ATA);
    cout<<"finished A setup"<<endl;
}

void BodyInterpolator::interpolateMesh(double p, Eigen::MatrixXd &V_updated) {

    double eps= 0.00001;
    if(p>= 1-eps)  { V_updated = V_to; return; }
    if(p<eps) { V_updated = V_from; return; }

    int numVert = V_from.rows();
    int numFace = F.rows();
    // the real interpolation
    V_updated= MatrixXd::Zero(numVert, 3);
    MatrixXd RHS_Q = MatrixXd::Zero(3 * numFace, 3);

    for (int f =0; f < numFace; f++){
        // interpolate S
        Vector3d K_inter = (1.0 - p) * Vector3d::Ones() + p * K[f];
        Matrix3d S_inter = U[f] * K_inter.asDiagonal() * U2[f].transpose();
        // interpolate R
        Quaterniond qI(Matrix3d::Identity());
        Quaterniond q_inter = qI.slerp(p, q[f]);
        MatrixXd R_inter = q_inter.toRotationMatrix();

        Matrix3d Q_inter = R_inter * S_inter;
        RHS_Q.block<3, 3>(3 * f, 0) = Q_inter.transpose();
    }
    // subtract the missing vertex and row from the right hand side
    // refer to page 68 of Robert Sumners thesis
    MatrixXd RHS = RHS_Q;
    Vector3d v_fixed = (1.0 - p) * V_from.row(fixed_id) + p * V_to.row(fixed_id);
    RHS.col(0) -= v_fixed(0) * A_remcol;
    RHS.col(1) -= v_fixed(1) * A_remcol;
    RHS.col(2) -= v_fixed(2) * A_remcol;

    // Build the linear system:
    RHS = A.transpose() * RHS;

    MatrixXd V_solved = cholSolver.solve(RHS);


    // plug in the one fixed vertex into the solution
    V_updated.block(0, 0, fixed_id, 3) = V_solved.block(0, 0, fixed_id, 3);
    V_updated.row(fixed_id) = v_fixed;
    V_updated.block(fixed_id + 1, 0, V_from.rows() - fixed_id - 1, 3) =
            V_solved.block(fixed_id, 0, V_from.rows() - fixed_id - 1, 3);
}