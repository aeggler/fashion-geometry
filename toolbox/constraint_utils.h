//
// Created by Anna Maria Eggler on 21.10.22.
//

#ifndef EXAMPLE_CONSTRAINT_UTILS_H
#define EXAMPLE_CONSTRAINT_UTILS_H
#include <Eigen/Dense>
#endif //EXAMPLE_CONSTRAINT_UTILS_H

/**
 * @brief Computes an ideal *reflection* between points1 and points2 using
 * procrustean analysis.
 * Usage:
 * FROM 2 TO 1
 * points3t = points2.transpose().colwise() - T_est;
 * points3t = (R_est.transpose() * points3);
 * FROM 1 TO 2
 * points3t = (R_est * points1.transpose());
 * points3t = points3t.colwise() + T_est;
 * points3 = points3t.transpose();
 */
 //https://igl.ethz.ch/projects/ARAP/svd_rot.pdf

    void procrustes(const Eigen::MatrixXd &points1,
                    const Eigen::MatrixXd &points2,
                    Eigen::MatrixXd &Rot_est,
                    Eigen::VectorXd &Trans_est
    );
void procrustesWORef(const Eigen::MatrixXd &points1,
                const Eigen::MatrixXd &points2,
                Eigen::MatrixXd &Rot_est,
                Eigen::VectorXd &Trans_est
);

void initProcrustesPatternTo3D(const Eigen::MatrixXd& Vg_pattern,const Eigen::MatrixXi& Fg_pattern,
                           const Eigen::MatrixXi& Fg_orig, const Eigen::MatrixXd& p, Eigen::MatrixXd& procrustesPatternIn3D);

