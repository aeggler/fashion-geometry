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
    sigma(sigma.rows() - 1, sigma.cols() - 1) = (svd.matrixV() * svd.matrixU().transpose()).determinant();

    R_est = svd.matrixV() * sigma * svd.matrixU().transpose();
    //R_est = MatrixXd::Identity(3, 3);
    T_est = qb - R_est * pb;
}