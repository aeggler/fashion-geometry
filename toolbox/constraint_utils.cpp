//
// Created by Anna Maria Eggler on 21.10.22.
//

#include "constraint_utils.h"

#include <Eigen/Dense>
#include <iostream>
#include <Eigen/SVD>

using namespace Eigen;
using namespace std;
// see https://igl.ethz.ch/projects/ARAP/svd_rot.pdf

void procrustes(const Eigen::MatrixXd& points1,    // from
                const Eigen::MatrixXd& points2,    //to
                Eigen::MatrixXd& Rot_est,
                Eigen::VectorXd& Trans_est
                ){
    MatrixXd points1t= points1.transpose();
    MatrixXd points2t= points2.transpose();

    VectorXd pcentroid = points1t.rowwise().mean();
    VectorXd qcentroid = points2t.rowwise().mean();

    MatrixXd X = (points1t.colwise() - pcentroid);
    MatrixXd Y = (points2t.colwise() - qcentroid);
    MatrixXd S = X * Y.transpose();

    JacobiSVD<MatrixXd> svd (S, ComputeThinU | Eigen::ComputeThinV);
    MatrixXd sigma = MatrixXd::Identity(svd.matrixU().cols(), svd.matrixV().cols());
    int srows = sigma.rows();
    int scols = sigma.cols();

    sigma(srows-1, scols-1) = - -(svd.matrixV() * svd.matrixU().transpose()).determinant();
    Rot_est = svd.matrixV() * sigma * svd.matrixU().transpose();
    Trans_est = qcentroid - Rot_est * pcentroid;
}