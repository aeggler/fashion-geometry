//
// Created by Anna Maria Eggler on 31.10.22.
//

#ifndef EXAMPLE_GARMENT_ADAPTION_H
#define EXAMPLE_GARMENT_ADAPTION_H
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Eigen/Geometry>
using namespace std;
using namespace Eigen;


class garment_adaption {
private:
    int numFace;
    int numVert;
    std::vector<Eigen::MatrixXd > jacobians;
    std::vector<Eigen::MatrixXd > inv_jacobians;
    Eigen::MatrixXd V_init;
    Eigen::MatrixXd V;
    Eigen::MatrixXi Fg;
    Eigen::MatrixXd V_pattern;
    Eigen::MatrixXi Fg_pattern;


public:
    garment_adaption(Eigen::MatrixXd& V, Eigen::MatrixXi& Fg, Eigen::MatrixXd & V_pattern, Eigen::MatrixXi& Fg_pattern_orig);
    void computeJacobian();
    void setUpRotationMatrix(double angle,Vector3d& axis, Matrix4d& rotationMatrix);
    void performJacobianUpdateAndMerge(Eigen::MatrixXd & V_curr, int iteratitons);

};


#endif //EXAMPLE_GARMENT_ADAPTION_H
