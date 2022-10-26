//
// Created by Anna Maria Eggler on 26.10.22.
//
// interpolation according to http://people.csail.mit.edu/sumner/thesis/Sumner2005MMU.pdf
// Deformation Transfer for Triangle Meshes
//Robert W. Sumner Jovan Popovic ́

//Poisson shape interpolation
//Dong Xu, Hongxin Zhang, Qing Wang, Hujun Bao*

// for polar decomp see https://research.cs.wisc.edu/graphics/Courses/838-s2002/Papers/polar-decomp.pdf
// based on custom fit garments by Katja Wolff

#ifndef EXAMPLE_BODY_INTERPOLATION_H
#define EXAMPLE_BODY_INTERPOLATION_H

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace Eigen; using namespace std;

class BodyInterpolator{
private:
    MatrixXd V_from, V_to;
    MatrixXi F;

    // for the pre-computation of R, S , R*S= Q, a 3x3 matrix for each face
    std::vector< Eigen::Matrix3d > U;	 //  S can be factored into diagonal form, S = U*K*U2.transpose()
    std::vector< Eigen::Matrix3d > U2;
    std::vector< Eigen::Vector3d > K;
    std::vector< Eigen::Quaterniond > q; // quaternion of R

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholSolver;
//TODO REWRITE DOC
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXi A_data;			 // store the non-zero index data for A
    Eigen::SparseMatrix<double> ATA; // we need to store this here, since cg only keeps a reference to this matrix

    int fixed_id;				// id of the one vertex that gets fixed in space
    Eigen::VectorXd A_remcol;	// removed column of A for that fixed vertex


public:
    BodyInterpolator(MatrixXd& V_from, MatrixXd& V_to, MatrixXi& F );
    void interpolateMesh(double p, MatrixXd& V_updated);
};
#endif //EXAMPLE_BODY_INTERPOLATION_H