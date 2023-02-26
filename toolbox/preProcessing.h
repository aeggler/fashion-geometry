//
// Created by Anna Maria Eggler on 30.01.23.
//

#ifndef EXAMPLE_PREPROCESSING_H
#define EXAMPLE_PREPROCESSING_H

#endif //EXAMPLE_PREPROCESSING_H
#include <Eigen/Dense>
#include <iostream>
#include <map>


using namespace std;
using namespace Eigen;

void initCollMeshCall( MatrixXd& Vm_left, MatrixXi& Fm_left,
                       MatrixXd& Vm_right, MatrixXi& Fm_right);

void setupCollisionConstraintsCall(Eigen::MatrixXi& collisionVert, vector<int> & pureCollVert,
                                   MatrixXd& testMorph_V1left,
                                   MatrixXi& testMorph_F1left,
                                   MatrixXd& testMorph_V1right,
                                   MatrixXi& testMorph_F1right,
                                   MatrixXd& p, int& numVert, double coll_EPS,
                               std::map<int, int> & leftHalfToFullFaceMap, std::map<int, int> & rightHalfToFullFaceMap, vector<VectorXd> & CleftRight , vector<VectorXd>& NleftRight, VectorXi& closestFaceId,
                               MatrixXd& Vm, MatrixXi& Fm, MatrixXi& Fg);

void createHalfSewingPattern(MatrixXd& Vg, MatrixXi& Fg, MatrixXd& Vg_pattern, MatrixXi& Fg_pattern, MatrixXd& Vg_pattern_half, MatrixXi& Fg_pattern_half,
                             map<int, int>& halfPatternFaceToFullPatternFace, map<int, int>& fullPatternFaceToHalfPatternFace, map<int, int>& halfPatternVertToFullPatternVert ,
                             map<int, int>& fullPatternVertToHalfPatternVert, map<int, int>& insertedIdxToPatternVert, VectorXi& isLeftVertPattern , MatrixXd& R_sym,
                             VectorXd& T_sym, MatrixXd& rightVert);
//void oneShotLengthSolve(MatrixXd& p_adaption, MatrixXi& Fg_pattern_curr, MatrixXd& baryCoordsUPattern, MatrixXd& baryCoordsVPattern, MatrixXd& mapFromVg, MatrixXi& mapFromFg);
