//
// Created by Anna Maria Eggler on 21.10.22.
//

#include "constraint_utils.h"

#include <Eigen/Dense>
#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/LU>

using namespace Eigen;
using namespace std;
// see https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
//  ATTENTION THIS INCLUDES A REFLECTION
void procrustes(const Eigen::MatrixXd& points1,    // from
                const Eigen::MatrixXd& points2,    //to
                Eigen::MatrixXd& R_est,
                Eigen::VectorXd& T_est
                ){
    Eigen::MatrixXd points1t = points1.transpose();
    Eigen::MatrixXd points2t = points2.transpose();

    Eigen::VectorXd pb = points1t.rowwise().mean();
    Eigen::VectorXd qb = points2t.rowwise().mean();

    Eigen::MatrixXd X = (points1t.colwise() - pb);
    Eigen::MatrixXd Y = (points2t.colwise() - qb);
    Eigen::MatrixXd S = X * Y.transpose();

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::MatrixXd sigma = Eigen::MatrixXd::Identity(svd.matrixV().cols(), svd.matrixU().cols());
    sigma(sigma.rows() - 1, sigma.cols() - 1) = -(svd.matrixV() * svd.matrixU().transpose()).determinant();

    R_est = svd.matrixV() * sigma * svd.matrixU().transpose();
    //R_est = MatrixXd::Identity(3, 3);
    T_est = qb - R_est * pb;
}

void initProcrustesPatternTo3D( const Eigen::MatrixXd& Vg_pattern,const Eigen::MatrixXi& Fg_pattern,
                               const Eigen::MatrixXi& Fg_orig, const Eigen::MatrixXd& p, Eigen::MatrixXd& procrustesPatternIn3D ){

    int numFace = Fg_orig.rows();
    procrustesPatternIn3D.resize(numFace*3, 3);

    for(int j=0; j < numFace; j++ ){
        // we map each triangle individually

        Matrix3d fromMat;
        fromMat.row(0) = Vg_pattern.row(Fg_pattern(j, 0));
        fromMat.row(1) = Vg_pattern.row(Fg_pattern(j, 1));
        fromMat.row(2) = Vg_pattern.row(Fg_pattern(j, 2));

        int id0 = Fg_orig(j, 0);
        int id1 = Fg_orig(j, 1);
        int id2 = Fg_orig(j, 2);

        Matrix3d toMat;
        toMat.row(0)= p.row(id0);
        toMat.row(1)= p.row(id1);
        toMat.row(2)= p.row(id2);

        Eigen::MatrixXd R_est;
        Eigen::VectorXd T_est;
        procrustes( fromMat , toMat,  R_est, T_est);

        Eigen::MatrixXd rotated = (R_est * fromMat.transpose());

        Eigen::MatrixXd rotandTransT  = rotated.colwise()+ T_est;
        Eigen::MatrixXd rotandTrans = rotandTransT.transpose();

        procrustesPatternIn3D.row(3*j+0 ) = rotandTrans.row(0);
        procrustesPatternIn3D.row(3*j+1 ) = rotandTrans.row(1);
        procrustesPatternIn3D.row(3*j+2 ) = rotandTrans.row(2);

    }

}