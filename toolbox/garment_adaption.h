//
// Created by Anna Maria Eggler on 31.10.22.
//

#ifndef EXAMPLE_GARMENT_ADAPTION_H
#define EXAMPLE_GARMENT_ADAPTION_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include "seam.h"
#include <Eigen/Geometry>
using namespace std;
using namespace Eigen;


class garment_adaption {
private:
    int numFace;
    int numVertGarment;
    int numVertPattern;

    std::vector<Eigen::MatrixXd > inv_jacobians;
    Eigen::MatrixXd V_init;
    Eigen::MatrixXd V;
    Eigen::MatrixXi Fg;
    Eigen::MatrixXd V_pattern;
    Eigen::MatrixXi Fg_pattern;
    std::vector< std::vector<int> > faceFaceAdjecencyList_3D;
    vector<vector<int> > vfAdj;
    Eigen::SparseMatrix<double, RowMajor> A;
    Eigen::SparseMatrix<double, RowMajor> W;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholSolver;


public:
    /* Initialiize the class instance and set up matrix A for the local global computation */
    garment_adaption(Eigen::MatrixXd& V, Eigen::MatrixXi& Fg, Eigen::MatrixXd & V_pattern, Eigen::MatrixXi& Fg_pattern_orig,
                     vector<seam*>& seamsList, std::vector<std::vector<int>>& boundaryL
    );
    /* A area weighted average over the neighbors to smooth the target Jacobian*/
    void smoothJacobian();
    /* Computes per face a jacobian matrix and saves the target norm for visualization. The jacobian is computed based on initial vertex positions
     * of the 3D garment and the corresponding 2D pattern. */
    void computeJacobian();
    /* Even if the norms of the triangles are aligned, we have to rotate along the normal to align u (or v) direction before we apply the inverse jacobian
     * Align EITHER u OR v direction. By default it is v. */
    void setUpRotationMatrix(double angle,Vector3d& axis, Matrix4d& rotationMatrix);
    /* Use local global to compute new positions of the pattern vertices, given the current 3D model. Input iterations defines how many local global steps are
     * performed. Bary coords are needed to align u or v direction. Output: v_newPattern  */
    void performJacobianUpdateAndMerge(Eigen::MatrixXd & V_curr, int iteratitons, const MatrixXd& baryCoords1, const MatrixXd& baryCoords2,
                                       Eigen::MatrixXd & V_newPattern, vector<seam*>& seamsList, std::vector<std::vector<int>>& boundaryL);
    /* additional feature to slim or widen a garment in a specific area byy a specific amount */
    void changeFitViaJacobian(bool geoDistU,bool geoDistV,double geoDistChange, const Eigen::VectorXd& affectedFaces );

    std::vector<std::pair<double, double>> perFaceTargetNorm;
    std::vector<Eigen::MatrixXd > jacobians;
};


#endif //EXAMPLE_GARMENT_ADAPTION_H
